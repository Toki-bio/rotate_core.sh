#!/bin/bash
# rotate_core.sh - v13 FINAL - Core-based rotation using ssearch36
#
# Rotates satellite DNA sequences to align at conserved core region
# Uses ssearch36 to find best alignment, handles both strands
#
# Usage: rotate_core.sh <input.fa> <anchor_name>
# Output: <input>.rotated.fa

set -uo pipefail

[[ $# -ge 2 ]] || { echo "Usage: $0 <input.fa> <anchor_name>"; exit 1; }

INFA="$1"
ANCHOR_NAME="$2"
OUTFA="${INFA%.fa}.rotated.fa"
[[ "$OUTFA" == "$INFA" ]] && OUTFA="${INFA}.rotated"

command -v ssearch36 >/dev/null || { echo "Error: ssearch36 not found"; exit 1; }
command -v seqkit >/dev/null || { echo "Error: seqkit not found"; exit 1; }

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR" 2>/dev/null || true' EXIT

echo "=== rotate_core.sh v13 ==="
echo "Input  : $INFA"
echo "Anchor : $ANCHOR_NAME"
echo "Output : $OUTFA"
echo

# --- Preprocessing ---
tr -d '\r' < "$INFA" | sed '/^$/d' > "$TMPDIR/in.fa"
seqkit seq -g -w 0 "$TMPDIR/in.fa" > "$TMPDIR/clean.fa"

# --- Extract anchor ---
ANCHOR_SEQ=$(awk -v name="$ANCHOR_NAME" '
    /^>/ { hdr=substr($0,2); split(hdr,a," "); current=a[1]; next }
    current == name { print; exit }
' "$TMPDIR/clean.fa")

[[ -z "$ANCHOR_SEQ" ]] && { echo "ERROR: Anchor '$ANCHOR_NAME' not found!"; exit 1; }

ANCHOR_LEN=${#ANCHOR_SEQ}
printf ">%s\n%s\n" "$ANCHOR_NAME" "$ANCHOR_SEQ" > "$TMPDIR/anchor.fa"
echo "Anchor: $ANCHOR_NAME ($ANCHOR_LEN bp)"

TOTAL=$(grep -c "^>" "$TMPDIR/clean.fa")
echo "Sequences: $TOTAL"
echo

# --- Create sequence list (critical: file-based to avoid subshell issues) ---
awk '/^>/ { 
    if (seq != "") print name "\t" seq
    name = substr($0, 2); split(name, a, " "); name = a[1]; seq = ""
    next
}
{ seq = seq $0 }
END { if (seq != "") print name "\t" seq }
' "$TMPDIR/clean.fa" > "$TMPDIR/seqlist.tsv"

# --- PASS 1: Align anchor vs doubled sequences ---
echo "=== PASS 1: Aligning ===" 

> "$TMPDIR/hits.tsv"
SEQ_NUM=0

while IFS=$'\t' read -r NAME SEQ; do
    SEQ_NUM=$((SEQ_NUM + 1))
    LEN=${#SEQ}
    
    [[ "$NAME" == "$ANCHOR_NAME" ]] && { echo "[$SEQ_NUM/$TOTAL] $NAME -> anchor"; continue; }
    
    printf ">d\n%s%s\n" "$SEQ" "$SEQ" > "$TMPDIR/dbl.fa"
    
    ssearch36 -Q -n -z 11 -E 2 -m 8C "$TMPDIR/anchor.fa" "$TMPDIR/dbl.fa" 3 2>/dev/null \
        | grep -v "^#" > "$TMPDIR/hit.m8" || true
    
    if [[ -s "$TMPDIR/hit.m8" ]]; then
        sort -k12,12nr "$TMPDIR/hit.m8" | head -1 > "$TMPDIR/best.tmp"
        read -r _ _ _ _ _ _ qs qe ss se _ score _ < "$TMPDIR/best.tmp" || true
        strand="+"; (( qs > qe )) && strand="-"
        printf "%d\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n" \
            "$SEQ_NUM" "$NAME" "$LEN" "$qs" "$qe" "$ss" "$se" "$strand" "$score" >> "$TMPDIR/hits.tsv"
        echo "[$SEQ_NUM/$TOTAL] $NAME -> $strand $score"
    else
        printf "%d\t%s\t%d\t0\t0\t0\t0\t?\t0\n" "$SEQ_NUM" "$NAME" "$LEN" >> "$TMPDIR/hits.tsv"
        echo "[$SEQ_NUM/$TOTAL] $NAME -> no hit"
    fi
done < "$TMPDIR/seqlist.tsv"

# --- PASS 2: Determine anchor offset ---
echo
echo "=== PASS 2: Anchor offset ==="

awk -F'\t' '$8 != "?" {print $9 "\t" $0}' "$TMPDIR/hits.tsv" | sort -k1,1nr | head -1 | cut -f2- > "$TMPDIR/best_hit.tsv"
BEST=$(cat "$TMPDIR/best_hit.tsv")

if [[ -z "$BEST" ]]; then
    ANCHOR_OFFSET=0
    echo "No valid hits -> offset=0"
else
    qs=$(echo "$BEST" | cut -f4)
    qe=$(echo "$BEST" | cut -f5)
    strand=$(echo "$BEST" | cut -f8)
    score=$(echo "$BEST" | cut -f9)
    
    [[ "$strand" == "-" ]] && eff_q=$qe || eff_q=$qs
    ANCHOR_OFFSET=$((eff_q - 1))
    echo "Best: score=$score strand=$strand -> offset=$ANCHOR_OFFSET"
fi

# Rotate anchor
if (( ANCHOR_OFFSET > 0 && ANCHOR_OFFSET < ANCHOR_LEN )); then
    ANCHOR_ROTATED="${ANCHOR_SEQ:$ANCHOR_OFFSET}${ANCHOR_SEQ:0:$ANCHOR_OFFSET}"
else
    ANCHOR_ROTATED="$ANCHOR_SEQ"
fi
printf ">%s\n%s\n" "$ANCHOR_NAME" "$ANCHOR_ROTATED" > "$OUTFA"

# --- PASS 3: Rotate all sequences ---
echo
echo "=== PASS 3: Rotating ==="

SEQ_NUM=0
while IFS=$'\t' read -r NAME SEQ; do
    SEQ_NUM=$((SEQ_NUM + 1))
    LEN=${#SEQ}
    
    [[ "$NAME" == "$ANCHOR_NAME" ]] && continue
    
    HIT=$(awk -F'\t' -v n="$SEQ_NUM" '$1 == n {print; exit}' "$TMPDIR/hits.tsv") || true
    
    if [[ -z "$HIT" ]]; then
        printf ">%s\n%s\n" "$NAME" "$SEQ" >> "$OUTFA"
        continue
    fi
    
    qs=$(echo "$HIT" | cut -f4)
    qe=$(echo "$HIT" | cut -f5)
    ss=$(echo "$HIT" | cut -f6)
    se=$(echo "$HIT" | cut -f7)
    strand=$(echo "$HIT" | cut -f8)
    
    if [[ "$strand" == "?" || -z "$strand" ]]; then
        printf ">%s\n%s\n" "$NAME" "$SEQ" >> "$OUTFA"
        continue
    fi
    
    if [[ "$strand" == "-" ]]; then
        SEQ_USE=$(printf ">_\n%s\n" "$SEQ" | seqkit seq -rp -t dna 2>/dev/null | seqkit seq -s -w 0 2>/dev/null) || {
            printf ">%s\n%s\n" "$NAME" "$SEQ" >> "$OUTFA"; continue
        }
        pos=$(( (se - 1) % LEN ))
        start=$(( LEN - 1 - pos + ANCHOR_OFFSET ))
    else
        SEQ_USE="$SEQ"
        start=$(( ss - qs + ANCHOR_OFFSET ))
    fi
    
    while (( start < 0 )); do start=$((start + LEN)); done
    start=$((start % LEN))
    
    ROTATED="${SEQ_USE:$start}${SEQ_USE:0:$start}"
    printf ">%s\n%s\n" "$NAME" "$ROTATED" >> "$OUTFA"
    
done < "$TMPDIR/seqlist.tsv"

echo
echo "=== Done ==="
echo "Output: $OUTFA ($(grep -c '^>' "$OUTFA") sequences)"
echo "Next: mafft --localpair --maxiterate 1000 $OUTFA > ${OUTFA%.fa}.al"
