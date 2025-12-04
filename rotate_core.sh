#!/bin/bash
# rotate_core.sh - v15 - FIXED rotation logic
# Fixed: correct rotation using (ss - 1 + ANCHOR_OFFSET) % LEN
# No longer uses incorrect (ss - qs) formula

set -uo pipefail

[[ $# -ge 2 ]] || { echo "Usage: $0 <input.fa> <anchor_name>"; exit 1; }

INFA="$1"
ANCHOR_NAME="$2"
OUTFA="${INFA%.fa}.rotated.fa"
[[ "$OUTFA" == "$INFA" ]] && OUTFA="${INFA}.rotated.fa"
OUTAL="${OUTFA%.fa}.al"

command -v ssearch36 >/dev/null || { echo "Error: ssearch36 not found"; exit 1; }
command -v seqkit >/dev/null || { echo "Error: seqkit not found"; exit 1; }

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR" 2>/dev/null || true' EXIT

echo "=== rotate_core.sh v15 ==="
echo "Input  : $INFA"
echo "Anchor : $ANCHOR_NAME"
echo "Output : $OUTFA"
echo

# --- Clean input ---
tr -d '\r' < "$INFA" | sed '/^$/d' > "$TMPDIR/in.fa"
seqkit seq -g -w 0 "$TMPDIR/in.fa" > "$TMPDIR/clean.fa"

# --- Extract anchor sequence ---
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

# --- Build sequence list: shortname \t fullheader \t sequence ---
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

# --- PASS 1: Align anchor vs doubled sequences ---
echo "=== PASS 1: Aligning ===" 

> "$TMPDIR/hits.tsv"
SEQ_NUM=0

while IFS=$'\t' read -r SHORTNAME FULLHDR SEQ; do
    SEQ_NUM=$((SEQ_NUM + 1))
    LEN=${#SEQ}
    
    [[ "$SHORTNAME" == "$ANCHOR_NAME" ]] && { echo "[$SEQ_NUM/$TOTAL] $SHORTNAME -> anchor"; continue; }
    
    printf ">d\n%s%s\n" "$SEQ" "$SEQ" > "$TMPDIR/dbl.fa"
    
    ssearch36 -Q -n -z 11 -E 2 -m 8C "$TMPDIR/anchor.fa" "$TMPDIR/dbl.fa" 3 2>/dev/null \
        | grep -v "^#" > "$TMPDIR/hit.m8" || true
    
    if [[ -s "$TMPDIR/hit.m8" ]]; then
        sort -k12,12nr "$TMPDIR/hit.m8" | head -1 > "$TMPDIR/best.tmp"
        read -r _ _ _ _ _ _ qs qe ss se _ score _ < "$TMPDIR/best.tmp" || true
        strand="+"; (( qs > qe )) && strand="-"
        printf "%d\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n" \
            "$SEQ_NUM" "$SHORTNAME" "$LEN" "$qs" "$qe" "$ss" "$se" "$strand" "$score" >> "$TMPDIR/hits.tsv"
        echo "[$SEQ_NUM/$TOTAL] $SHORTNAME -> $strand $score"
    else
        printf "%d\t%s\t%d\t0\t0\t0\t0\t?\t0\n" "$SEQ_NUM" "$SHORTNAME" "$LEN" >> "$TMPDIR/hits.tsv"
        echo "[$SEQ_NUM/$TOTAL] $SHORTNAME -> no hit"
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
while IFS=$'\t' read -r SHORTNAME FULLHDR SEQ; do
    SEQ_NUM=$((SEQ_NUM + 1))
    LEN=${#SEQ}
    
    [[ "$SHORTNAME" == "$ANCHOR_NAME" ]] && continue
    
    HIT=$(awk -F'\t' -v n="$SEQ_NUM" '$1 == n {print; exit}' "$TMPDIR/hits.tsv") || true
    
    if [[ -z "$HIT" ]]; then
        printf ">%s\n%s\n" "$FULLHDR" "$SEQ" >> "$OUTFA"
        continue
    fi
    
    qs=$(echo "$HIT" | cut -f4)
    qe=$(echo "$HIT" | cut -f5)
    ss=$(echo "$HIT" | cut -f6)
    se=$(echo "$HIT" | cut -f7)
    strand=$(echo "$HIT" | cut -f8)
    
    if [[ "$strand" == "?" || -z "$strand" ]]; then
        printf ">%s\n%s\n" "$FULLHDR" "$SEQ" >> "$OUTFA"
        continue
    fi
    
    if [[ "$strand" == "-" ]]; then
        SEQ_USE=$(printf ">_\n%s\n" "$SEQ" | seqkit seq -rp -t dna 2>/dev/null | seqkit seq -s -w 0 2>/dev/null) || {
            printf ">%s\n%s\n" "$FULLHDR" "$SEQ" >> "$OUTFA"; continue
        }
        # FIXED: use (se - 1 + ANCHOR_OFFSET) % LEN
        start=$(( (se - 1 + ANCHOR_OFFSET) % LEN ))
    else
        SEQ_USE="$SEQ"
        # FIXED: use (ss - 1 + ANCHOR_OFFSET) % LEN
        start=$(( (ss - 1 + ANCHOR_OFFSET) % LEN ))
    fi
    
    while (( start < 0 )); do start=$((start + LEN)); done
    start=$((start % LEN))
    
    ROTATED="${SEQ_USE:$start}${SEQ_USE:0:$start}"
    printf ">%s\n%s\n" "$FULLHDR" "$ROTATED" >> "$OUTFA"
    
done < "$TMPDIR/seqlist.tsv"

# --- PASS 4: MAFFT alignment ---
if command -v mafft >/dev/null; then
    echo
    echo "=== PASS 4: MAFFT alignment ==="
    mafft --thread $(nproc) --threadtb $(nproc) --threadit $(nproc) \
        --localpair --maxiterate 1000 --ep 0.123 \
        --nuc --reorder --preservecase --quiet "$OUTFA" > "$TMPDIR/mafft.al"

    # Restore full headers (MAFFT truncates at first space)
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
    ' "$TMPDIR/hdr_map.tsv" "$TMPDIR/mafft.al" > "$OUTAL"
    
    echo "Aligned: $OUTAL"
else
    echo "mafft not found - skipping alignment"
fi

echo
echo "=== Done ==="
echo "Rotated: $OUTFA ($(grep -c '^>' "$OUTFA") sequences)"
