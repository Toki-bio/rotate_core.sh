#!/bin/bash
# rotate_core.sh - v16 - Alignment-based universal rotation
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

echo "=== rotate_core.sh v16 ==="
echo "Input  : $INFA"
echo "Anchor : $ANCHOR_NAME"
echo "Output : $OUTFA"
echo

# --- Clean input ---
tr -d '\r' < "$INFA" | sed '/^$/d' > "$TMPDIR/in.fa"
seqkit seq -g -w 0 "$TMPDIR/in.fa" > "$TMPDIR/clean.fa"

TOTAL=$(grep -c "^>" "$TMPDIR/clean.fa")
echo "Sequences: $TOTAL"

# --- Extract anchor sequence ---
ANCHOR_SEQ=$(awk -v name="$ANCHOR_NAME" '
    /^>/ { hdr=substr($0,2); split(hdr,a," "); current=a[1]; next }
    current == name { print; exit }
' "$TMPDIR/clean.fa")

[[ -z "$ANCHOR_SEQ" ]] && { echo "ERROR: Anchor '$ANCHOR_NAME' not found!"; exit 1; }
ANCHOR_LEN=${#ANCHOR_SEQ}
echo "Anchor: $ANCHOR_NAME ($ANCHOR_LEN bp)"
echo

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

# Search anchor against all aligned sequences (doubled for circular)
# EXCLUDE anchor itself from search
> "$TMPDIR/all_doubled.fa"
awk -v anchor="$ANCHOR_NAME" '
    /^>/ { 
        if (seq != "" && name != anchor) {
            gsub(/-/, "", seq)  # remove gaps
            print ">" name
            print seq seq  # double it
        }
        name = substr($0, 2)
        split(name, a, " ")
        name = a[1]
        seq = ""
        next
    }
    { seq = seq $0 }
    END { 
        if (seq != "" && name != anchor) {
            gsub(/-/, "", seq)
            print ">" name
            print seq seq
        }
    }
' "$TMPDIR/initial.al" > "$TMPDIR/all_doubled.fa"

# ssearch36 to find best overall hit (excluding anchor)
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
# anchor[qs] ~ subject[ss], so anchor[1] ~ subject[ss - qs + 1]
ROTATION_POINT=$((ss - qs + 1))
echo "  Universal rotation point (in best hit): position $ROTATION_POINT"
echo

# --- STEP 3: Find alignment column for rotation point ---
echo "=== STEP 3: Map rotation point to alignment column ==="

# Get the best-hit sequence from alignment (with gaps)
BEST_SEQ_ALIGNED=$(awk -v name="$sid" '
    BEGIN { found = 0; seq = "" }
    /^>/ { 
        if (found) { print seq; exit }
        hdr = substr($0, 2)
        split(hdr, a, " ")
        if (a[1] == name) found = 1
        next 
    }
    found { seq = seq $0 }
    END { if (found) print seq }
' "$TMPDIR/initial.al")

if [[ -z "$BEST_SEQ_ALIGNED" ]]; then
    echo "ERROR: Could not find $sid in alignment"
    exit 1
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
}' <<< "$BEST_SEQ_ALIGNED")

echo "Best hit sequence: $sid"
echo "Rotation point (ungapped pos $ROTATION_POINT) = alignment column $ROTATION_COLUMN"
echo

# --- STEP 4: Rotate each sequence based on alignment column ---
echo "=== STEP 4: Rotating all sequences ==="

# Build sequence list with full headers
awk 'BEGIN {OFS="\t"}
    /^>/ { 
        if (seq != "") print short, full, seq
        full = substr($0, 2)
        split(full, a, " ")
        short = a[1]
        seq = ""
        next
    }
    { seq = seq $0 }
    END { if (seq != "") print short, full, seq }
' "$TMPDIR/clean.fa" > "$TMPDIR/seqlist.tsv"

# For each sequence, find its position at ROTATION_COLUMN in alignment
> "$OUTFA"

while IFS=$'\t' read -r SHORT FULL SEQ; do
    LEN=${#SEQ}
    
    # Get this sequence's aligned version (single line)
    SEQ_ALIGNED=$(awk -v name="$SHORT" '
        BEGIN { found = 0; seq = "" }
        /^>/ { 
            if (found) { print seq; exit }
            hdr = substr($0, 2)
            split(hdr, a, " ")
            if (a[1] == name) found = 1
            next 
        }
        found { seq = seq $0 }
        END { if (found) print seq }
    ' "$TMPDIR/initial.al" | tr -d '\n')
    
    if [[ -z "$SEQ_ALIGNED" ]]; then
        echo "  $SHORT: not found in alignment, keeping as-is"
        printf ">%s\n%s\n" "$FULL" "$SEQ" >> "$OUTFA"
        continue
    fi
    
    # Find ungapped position at ROTATION_COLUMN (single line input)
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
        echo "  $SHORT: invalid rotation position '$ROT_POS', keeping as-is"
        printf ">%s\n%s\n" "$FULL" "$SEQ" >> "$OUTFA"
        continue
    fi
    
    # Convert to 0-indexed
    START=$((ROT_POS - 1))
    
    # Handle edge cases
    if (( START < 0 )); then START=0; fi
    if (( START >= LEN )); then START=$((START % LEN)); fi
    
    # Rotate
    ROTATED="${SEQ:$START}${SEQ:0:$START}"
    
    echo "  $SHORT: rotate to position $ROT_POS (0-indexed: $START)"
    printf ">%s\n%s\n" "$FULL" "$ROTATED" >> "$OUTFA"
    
done < "$TMPDIR/seqlist.tsv"

echo
echo "Rotated sequences: $OUTFA"
echo

# --- STEP 5: Final MAFFT alignment ---
echo "=== STEP 5: Final MAFFT alignment ==="
mafft --thread $(nproc) --threadtb $(nproc) --threadit $(nproc) \
    --localpair --maxiterate 1000 --ep 0.123 \
    --nuc --reorder --preservecase --quiet "$OUTFA" > "$TMPDIR/final.al"

# Restore full headers
awk '/^>/{print substr($0,2)}' "$OUTFA" | while read -r hdr; do
    short=$(echo "$hdr" | cut -d' ' -f1)
    printf "%s\t%s\n" "$short" "$hdr"
done > "$TMPDIR/hdr_map.tsv"

awk -F'\t' 'NR==FNR {map[$1]=$2; next} 
    /^>/ {
        short=substr($0,2)
        if (short in map) print ">"map[short]
        else print
        next
    }
    {print}
' "$TMPDIR/hdr_map.tsv" "$TMPDIR/final.al" > "$OUTAL"

echo "Final alignment: $OUTAL"
echo

echo "=== Done ==="
echo "Rotated FASTA: $OUTFA ($(grep -c '^>' "$OUTFA") sequences)"
echo "Aligned: $OUTAL"
