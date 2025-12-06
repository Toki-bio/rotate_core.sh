#!/bin/bash
# rotate_core.sh - v22_safe - Alignment-based universal rotation with sequence integrity checks
# 
# Algorithm:
# 1. MAFFT align all sequences (unrotated)
# 2. ssearch36 anchor vs alignment → find best hit
# 3. Best hit's subject position (ss) defines the universal rotation column
# 4. For each sequence: find ungapped position at that column → rotate there
# 5. Re-align rotated sequences
#
# Usage: rotate_core.sh <input.fa> <anchor_name>
# Output: <input>.rotated.fa, <input>.rotated.al

set -uo pipefail

[[ $# -ge 2 ]] || { echo "Usage: $0 <input.fa> <anchor_name>"; exit 1; }

INFA="$1"
ANCHOR_NAME="$2"
OUTFA="${INFA%.fa}.rotated.fa"
[[ "$OUTFA" == "$INFA" ]] && OUTFA="${INFA}.rotated.fa"
OUTAL="${OUTFA%.fa}.al"

command -v ssearch36 >/dev/null || { echo "Error: ssearch36 not found"; exit 1; }
command -v seqkit >/dev/null || { echo "Error: seqkit not found"; exit 1; }
command -v mafft >/dev/null || { echo "Error: mafft not found"; exit 1; }

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR" 2>/dev/null || true' EXIT

echo "=== rotate_core.sh v22_safe ==="
echo "Input  : $INFA"
echo "Anchor : $ANCHOR_NAME"
echo "Output : $OUTFA"
echo

# --- Clean input ---
echo "=== Loading and validating input ==="
tr -d '\r' < "$INFA" | sed '/^$/d' > "$TMPDIR/in.fa"
seqkit seq -g -w 0 "$TMPDIR/in.fa" > "$TMPDIR/clean.fa"

TOTAL=$(grep -c "^>" "$TMPDIR/clean.fa")
echo "Sequences: $TOTAL"

# --- Validate all sequences contain only valid DNA characters ---
echo "Validating sequence characters..."
INVALID_SEQS=$(awk '
    /^>/ { next }
    {
        # Remove valid DNA characters (including ambiguous bases)
        line = toupper($0)
        gsub(/[ACGTURYSWKMBDHVN.-]/, "", line)
        if (line != "") {
            print "Invalid characters in sequence: " line
            print "Full line: " $0
        }
    }
' "$TMPDIR/clean.fa")

if [[ -n "$INVALID_SEQS" ]]; then
    echo "ERROR: Sequences contain invalid DNA characters:"
    echo "$INVALID_SEQS"
    exit 1
fi

# --- Extract anchor sequence ---
ANCHOR_SEQ=$(awk -v name="$ANCHOR_NAME" '
    /^>/ { 
        hdr=substr($0,2); 
        split(hdr,a," "); 
        current=a[1]; 
        next 
    }
    current == name { print; exit }
' "$TMPDIR/clean.fa")

[[ -z "$ANCHOR_SEQ" ]] && { echo "ERROR: Anchor '$ANCHOR_NAME' not found!"; exit 1; }
ANCHOR_LEN=${#ANCHOR_SEQ}
echo "Anchor: $ANCHOR_NAME ($ANCHOR_LEN bp)"
echo

# --- Save original sequences for verification ---
echo "=== Saving original sequences for verification ==="
awk 'BEGIN {OFS="\t"}
    /^>/ { 
        if (seq != "") print short, full, seq, length(seq)
        full = substr($0, 2)
        split(full, a, " ")
        short = a[1]
        seq = ""
        next
    }
    { seq = seq $0 }
    END { if (seq != "") print short, full, seq, length(seq) }
' "$TMPDIR/clean.fa" > "$TMPDIR/original_seqs.tsv"

echo "Saved original sequences: $TMPDIR/original_seqs.tsv"

# --- Create header mapping (short ID -> full header) ---
awk 'BEGIN {OFS="\t"}
    /^>/ { 
        if (full != "") print short, full
        full = substr($0, 2)
        split(full, a, " ")
        short = a[1]
        next
    }
    END { if (full != "") print short, full }
' "$TMPDIR/clean.fa" > "$TMPDIR/header_map.tsv"

# --- STEP 1: Initial MAFFT alignment ---
echo "=== STEP 1: Initial MAFFT alignment ==="
mafft --thread $(nproc) --threadtb $(nproc) --threadit $(nproc) \
    --localpair --maxiterate 1000 --ep 0.123 \
    --nuc --reorder --preservecase --quiet "$TMPDIR/clean.fa" > "$TMPDIR/initial.al"
echo "Done: $TMPDIR/initial.al"
echo

# --- STEP 2: Find best hit with ssearch36 ---
echo "=== STEP 2: Find universal rotation point ==="

# Extract anchor from alignment
printf ">%s\n%s\n" "$ANCHOR_NAME" "$ANCHOR_SEQ" > "$TMPDIR/anchor.fa"

# Create doubled sequences (EXCLUDE anchor) - use short IDs
awk -v anchor="$ANCHOR_NAME" '
    /^>/ {
        if (seq != "") {
            # Process previous sequence
            gsub(/-/, "", seq)  # remove gaps
            if (short_id != "" && short_id != anchor) {
                # Use short ID only for ssearch36
                print ">" short_id
                print seq seq  # double it
            }
        }
        # Parse current header - extract short ID
        full_header = substr($0, 2)
        split(full_header, a, " ")
        short_id = a[1]
        seq = ""
        next
    }
    { seq = seq $0 }
    END { 
        if (seq != "") {
            gsub(/-/, "", seq)
            if (short_id != "" && short_id != anchor) {
                print ">" short_id
                print seq seq
            }
        }
    }
' "$TMPDIR/initial.al" > "$TMPDIR/all_doubled.fa"

# ssearch36 - will output short IDs
ssearch36 -Q -n -z 11 -E 2 -m 8C "$TMPDIR/anchor.fa" "$TMPDIR/all_doubled.fa" 3 2>/dev/null \
    | grep -v "^#" | sort -k12,12nr | head -1 > "$TMPDIR/best_hit.m8"

if [[ ! -s "$TMPDIR/best_hit.m8" ]]; then
    echo "ERROR: No ssearch36 hits found!"
    exit 1
fi

read -r qid sid pident alen mm gap qs qe ss se eval score < "$TMPDIR/best_hit.m8"
echo "Best hit: $sid"
echo "  anchor[$qs:$qe] aligns to $sid[$ss:$se]"
echo "  score=$score, identity=$pident%"

# The rotation point: where does anchor position 1 map to?
ROTATION_POINT=$((ss - qs + 1))
echo "  Universal rotation point (in best hit): position $ROTATION_POINT"
echo

# --- STEP 3: Find alignment column for rotation point ---
echo "=== STEP 3: Map rotation point to alignment column ==="

# Find the best-hit sequence in alignment by matching the beginning of headers
BEST_SEQ_ALIGNED=$(awk -v pattern="$sid" '
    BEGIN { found = 0; seq = "" }
    /^>/ { 
        if (found) { 
            print seq
            exit
        }
        current_full = substr($0, 2)
        # Check if pattern matches the beginning of the header
        if (index(current_full, pattern) == 1) {
            found = 1
        }
        seq = ""
        next 
    }
    found { seq = seq $0 }
    END { 
        if (found) print seq
        else print ""  # not found
    }
' "$TMPDIR/initial.al")

if [[ -z "$BEST_SEQ_ALIGNED" ]]; then
    # If direct pattern matching fails, try to find by short ID (first word)
    SHORT_SID=$(echo "$sid" | awk '{print $1}')
    echo "Note: Matching truncated header '$sid', using short ID: '$SHORT_SID'"
    
    BEST_SEQ_ALIGNED=$(awk -v short_id="$SHORT_SID" '
        BEGIN { found = 0; seq = "" }
        /^>/ { 
            if (found) { 
                print seq
                exit
            }
            current_full = substr($0, 2)
            split(current_full, a, " ")
            current_short = a[1]
            if (current_short == short_id) {
                found = 1
            }
            seq = ""
            next 
        }
        found { seq = seq $0 }
        END { 
            if (found) print seq
            else print ""  # not found
        }
    ' "$TMPDIR/initial.al")
    
    if [[ -z "$BEST_SEQ_ALIGNED" ]]; then
        echo "ERROR: Could not find sequence matching '$sid' in alignment"
        exit 1
    fi
fi

# Find alignment column corresponding to ungapped position ROTATION_POINT
ROTATION_COLUMN=$(awk -v pos="$ROTATION_POINT" '
{
    ungapped = 0
    for (i = 1; i <= length($0); i++) {
        c = substr($0, i, 1)
        if (c != "-") {
            ungapped++
            if (ungapped == pos) {
                print i
                exit
            }
        }
    }
    # If we get here, the position is beyond the ungapped length
    print 1  # default to first position
}' <<< "$BEST_SEQ_ALIGNED")

echo "Found best hit sequence in alignment"
echo "Rotation point (ungapped pos $ROTATION_POINT) = alignment column $ROTATION_COLUMN"
echo

# --- STEP 4: Rotate each sequence based on alignment column ---
echo "=== STEP 4: Rotating all sequences with integrity checks ==="

> "$OUTFA"
> "$TMPDIR/rotation_log.tsv"
echo -e "Sequence\tOriginal_Length\tRotation_Position\tRotated_Length\tSame_Characters\tValid_Rotation" > "$TMPDIR/rotation_log.tsv"

while IFS=$'\t' read -r SHORT FULL SEQ LEN; do
    # Get this sequence's aligned version
    SEQ_ALIGNED=$(awk -v name="$SHORT" '
    BEGIN { found = 0; seq = "" }
    /^>/ { 
        if (found) { 
            print seq
            exit
        }
        current_full = substr($0, 2)
        split(current_full, a, " ")
        current_short = a[1]
        if (current_short == name) {
            found = 1
        }
        seq = ""
        next 
    }
    found { seq = seq $0 }
    END { 
        if (found) print seq
        else print ""  # not found
    }
' "$TMPDIR/initial.al" | tr -d '\n')
    
    if [[ -z "$SEQ_ALIGNED" ]]; then
        echo "  $SHORT: not found in alignment, keeping as-is"
        # Save original sequence
        printf ">%s\n%s\n" "$FULL" "$SEQ" >> "$OUTFA"
        echo -e "$SHORT\t$LEN\t0\t$LEN\tYES\tNO_ALIGNMENT" >> "$TMPDIR/rotation_log.tsv"
        continue
    fi
    
    # Find ungapped position at ROTATION_COLUMN
    ROT_POS=$(awk -v col="$ROTATION_COLUMN" '
    {
        ungapped = 0
        for (i = 1; i <= col && i <= length($0); i++) {
            c = substr($0, i, 1)
            if (c != "-") ungapped++
        }
        print (ungapped > 0) ? ungapped : 1
        exit
    }' <<< "$SEQ_ALIGNED")
    
    # Trim whitespace
    ROT_POS=$(echo "$ROT_POS" | tr -d '[:space:]')
    
    # Validate ROT_POS is a number
    if ! [[ "$ROT_POS" =~ ^[0-9]+$ ]]; then
        echo "  WARNING: $SHORT: invalid rotation position '$ROT_POS', keeping as-is"
        printf ">%s\n%s\n" "$FULL" "$SEQ" >> "$OUTFA"
        echo -e "$SHORT\t$LEN\tERROR\t$LEN\tYES\tINVALID_POS" >> "$TMPDIR/rotation_log.tsv"
        continue
    fi
    
    # Convert to 0-indexed
    START=$((ROT_POS - 1))
    
    # Handle edge cases
    if (( START < 0 )); then START=0; fi
    if (( START >= LEN )); then START=$((START % LEN)); fi
    
    # Rotate
    ROTATED="${SEQ:$START}${SEQ:0:$START}"
    ROTATED_LEN=${#ROTATED}
    
    # --- INTEGRITY CHECKS ---
    # Check 1: Length must be preserved
    if [[ "$LEN" -ne "$ROTATED_LEN" ]]; then
        echo "  ERROR: $SHORT: Length changed from $LEN to $ROTATED_LEN!"
        echo "  Keeping original sequence"
        printf ">%s\n%s\n" "$FULL" "$SEQ" >> "$OUTFA"
        echo -e "$SHORT\t$LEN\t$ROT_POS\t$ROTATED_LEN\tERROR\tLENGTH_CHANGE" >> "$TMPDIR/rotation_log.tsv"
        continue
    fi
    
    # Check 2: Character set must be identical (case-insensitive)
    ORIG_SORTED=$(echo "$SEQ" | tr '[:lower:]' '[:upper:]' | grep -o . | sort | tr -d '\n')
    ROT_SORTED=$(echo "$ROTATED" | tr '[:lower:]' '[:upper:]' | grep -o . | sort | tr -d '\n')
    
    if [[ "$ORIG_SORTED" != "$ROT_SORTED" ]]; then
        echo "  ERROR: $SHORT: Character composition changed!"
        echo "  Keeping original sequence"
        printf ">%s\n%s\n" "$FULL" "$SEQ" >> "$OUTFA"
        echo -e "$SHORT\t$LEN\t$ROT_POS\t$ROTATED_LEN\tNO\tCHARACTER_CHANGE" >> "$TMPDIR/rotation_log.tsv"
        continue
    fi
    
    # Check 3: Verify rotation mathematically (rotate back should give original)
    BACK_START=$((LEN - START))
    if (( BACK_START >= LEN )); then BACK_START=$((BACK_START % LEN)); fi
    ROTATED_BACK="${ROTATED:$BACK_START}${ROTATED:0:$BACK_START}"
    
    if [[ "$ROTATED_BACK" != "$SEQ" ]]; then
        echo "  ERROR: $SHORT: Rotation validation failed!"
        echo "  Keeping original sequence"
        printf ">%s\n%s\n" "$FULL" "$SEQ" >> "$OUTFA"
        echo -e "$SHORT\t$LEN\t$ROT_POS\t$ROTATED_LEN\tYES\tROTATION_VALIDATION_FAILED" >> "$TMPDIR/rotation_log.tsv"
        continue
    fi
    
    # Check 4: Ensure no gaps or invalid characters introduced
    if echo "$ROTATED" | grep -q '[^ACGTURYSWKMBDHVNacgturyswkmbdhvn]'; then
        echo "  ERROR: $SHORT: Invalid characters in rotated sequence!"
        echo "  Keeping original sequence"
        printf ">%s\n%s\n" "$FULL" "$SEQ" >> "$OUTFA"
        echo -e "$SHORT\t$LEN\t$ROT_POS\t$ROTATED_LEN\tNO\tINVALID_CHARACTERS" >> "$TMPDIR/rotation_log.tsv"
        continue
    fi
    
    # All checks passed
    echo "  $SHORT: rotate to position $ROT_POS (0-indexed: $START) [VALIDATED]"
    printf ">%s\n%s\n" "$FULL" "$ROTATED" >> "$OUTFA"
    echo -e "$SHORT\t$LEN\t$ROT_POS\t$ROTATED_LEN\tYES\tOK" >> "$TMPDIR/rotation_log.tsv"
    
done < "$TMPDIR/original_seqs.tsv"

echo
echo "Rotated sequences: $OUTFA"
echo "Rotation log: $TMPDIR/rotation_log.tsv"
echo

# --- Summary of rotation results ---
echo "=== Rotation Summary ==="
TOTAL_ROTATED=$(tail -n +2 "$TMPDIR/rotation_log.tsv" | wc -l)
SUCCESSFUL=$(tail -n +2 "$TMPDIR/rotation_log.tsv" | awk -F'\t' '$6 == "OK" {count++} END {print count+0}')
FAILED=$((TOTAL_ROTATED - SUCCESSFUL))

echo "Total sequences processed: $TOTAL_ROTATED"
echo "Successfully rotated: $SUCCESSFUL"
echo "Failed/kept original: $FAILED"

if [[ "$FAILED" -gt 0 ]]; then
    echo ""
    echo "Failed rotations:"
    tail -n +2 "$TMPDIR/rotation_log.tsv" | awk -F'\t' '$6 != "OK" {print $1 ": " $6}' | head -10
    if [[ "$FAILED" -gt 10 ]]; then
        echo "... and $((FAILED - 10)) more"
    fi
fi
echo

# --- Verify final output matches input sequences (rotated or original) ---
echo "=== Final verification ==="
awk 'BEGIN {OFS="\t"}
    /^>/ { 
        if (seq != "") print short, full, seq, length(seq)
        full = substr($0, 2)
        split(full, a, " ")
        short = a[1]
        seq = ""
        next
    }
    { seq = seq $0 }
    END { if (seq != "") print short, full, seq, length(seq) }
' "$OUTFA" > "$TMPDIR/final_seqs.tsv"

# Compare with originals
MISMATCHES=0
while IFS=$'\t' read -r SHORT FULL SEQ LEN; do
    # Get original sequence
    ORIG_SEQ=$(awk -F'\t' -v short="$SHORT" '$1 == short {print $3; exit}' "$TMPDIR/original_seqs.tsv")
    ORIG_LEN=$(awk -F'\t' -v short="$SHORT" '$1 == short {print $4; exit}' "$TMPDIR/original_seqs.tsv")
    
    if [[ -z "$ORIG_SEQ" ]]; then
        echo "  WARNING: $SHORT not found in original sequences"
        continue
    fi
    
    # Check if it's a valid rotation
    if [[ "$LEN" -eq "$ORIG_LEN" ]] && [[ "$ORIG_SEQ$ORIG_SEQ" == *"$SEQ"* ]]; then
        # Valid rotation
        :
    else
        echo "  ERROR: $SHORT final sequence is not a valid rotation of original!"
        echo "    Original: $ORIG_SEQ"
        echo "    Final:    $SEQ"
        MISMATCHES=$((MISMATCHES + 1))
    fi
done < "$TMPDIR/final_seqs.tsv"

if [[ "$MISMATCHES" -eq 0 ]]; then
    echo "All sequences verified: rotations are valid circular permutations"
else
    echo "ERROR: $MISMATCHES sequences failed final verification!"
    exit 1
fi
echo

# --- STEP 5: Final MAFFT alignment ---
if [[ -s "$OUTFA" ]]; then
    echo "=== STEP 5: Final MAFFT alignment ==="
    mafft --thread $(nproc) --threadtb $(nproc) --threadit $(nproc) \
        --localpair --maxiterate 1000 --ep 0.123 \
        --nuc --reorder --preservecase --quiet "$OUTFA" > "$TMPDIR/final.al"
    
    # Restore full headers
    awk '/^>/{print substr($0,2)}' "$OUTFA" | while read -r hdr; do
        short=$(echo "$hdr" | cut -d' ' -f1)
        printf "%s\t%s\n" "$short" "$hdr"
    done > "$TMPDIR/hdr_map_final.tsv"
    
    awk -F'\t' 'NR==FNR {map[$1]=$2; next} 
        /^>/ {
            short=substr($0,2)
            if (short in map) print ">"map[short]
            else print
            next
        }
        {print}
    ' "$TMPDIR/hdr_map_final.tsv" "$TMPDIR/final.al" > "$OUTAL"
    
    echo "Final alignment: $OUTAL"
    echo
else
    echo "ERROR: No sequences in output FASTA!"
    exit 1
fi

echo "=== Done ==="
echo "Rotated FASTA: $OUTFA ($(grep -c '^>' "$OUTFA") sequences)"
echo "Aligned: $OUTAL"
echo
echo "=== Safety checks performed ==="
echo "1. Sequence character validation"
echo "2. Length preservation check"
echo "3. Character composition check"
echo "4. Rotation mathematical validation"
echo "5. Final rotation verification"
echo
echo "All sequences are guaranteed to be circular permutations of originals."
