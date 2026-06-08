#!/bin/bash
# rotate_to_anchor_v2.sh - Rotate sequences to match anchor using doubled-sequence alignment
# Fix: degap all sequences before processing (input may be aligned FASTA with '-' gaps)
# Usage: rotate_to_anchor_v2.sh <cluster.fa> <anchor_name>
# Output: <cluster.fa>.rotated (ungapped rotated sequences)

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <cluster.fa> <anchor_name>"
    echo "Output: <cluster.fa>.rotated"
    exit 1
fi

INFA="$1"
ANCHOR_NAME="$2"
OUTFA="${INFA}.rotated"

command -v water >/dev/null 2>&1 || { echo "Error: EMBOSS water not found"; exit 1; }
command -v seqkit >/dev/null 2>&1 || { echo "Error: seqkit not found"; exit 1; }

TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Extract anchor and DEGAP it
seqkit grep -p "$ANCHOR_NAME" "$INFA" | seqkit seq -w 0 > "$TMPDIR/anchor_raw.fa"

if [[ ! -s "$TMPDIR/anchor_raw.fa" ]]; then
    echo "Error: Anchor '$ANCHOR_NAME' not found in $INFA"
    exit 1
fi

# Degap anchor sequence
ANCHOR_SEQ=$(awk 'NR==2' "$TMPDIR/anchor_raw.fa" | tr -d '-' | tr -d '.')
ANCHOR_LEN=${#ANCHOR_SEQ}
printf ">%s\n%s\n" "$ANCHOR_NAME" "$ANCHOR_SEQ" > "$TMPDIR/anchor.fa"
echo "Anchor: $ANCHOR_NAME ($ANCHOR_LEN bp ungapped)"

# Initialize output with degapped anchor
cat "$TMPDIR/anchor.fa" > "$OUTFA"

# Get all non-anchor sequences
seqkit grep -v -p "$ANCHOR_NAME" "$INFA" | seqkit seq -w 0 > "$TMPDIR/others.fa"

TOTAL=$(grep -c "^>" "$TMPDIR/others.fa" || echo 0)
echo "Sequences to rotate: $TOTAL"

if [[ $TOTAL -eq 0 ]]; then
    echo "No other sequences to rotate"
    echo "Output: $OUTFA"
    exit 0
fi

COUNT=0
while IFS= read -r header; do
    IFS= read -r seq_raw

    COUNT=$((COUNT + 1))
    NAME=$(echo "$header" | sed 's/^>//' | awk '{print $1}')

    # DEGAP: strip gap characters before all processing
    seq=$(echo "$seq_raw" | tr -d '-' | tr -d '.')
    LEN=${#seq}

    echo -n "[$COUNT/$TOTAL] $NAME ($LEN bp ungapped): "

    # Double the sequence (forward)
    DOUBLED_FWD="${seq}${seq}"
    printf ">doubled_fwd\n%s\n" "$DOUBLED_FWD" > "$TMPDIR/doubled_fwd.fa"

    # Reverse complement and double it
    SEQ_RC=$(echo "$seq" | tr 'ACGTacgtNn' 'TGCAtgcaNn' | rev)
    DOUBLED_RC="${SEQ_RC}${SEQ_RC}"
    printf ">doubled_rc\n%s\n" "$DOUBLED_RC" > "$TMPDIR/doubled_rc.fa"

    # Align anchor vs forward doubled
    water -asequence "$TMPDIR/anchor.fa" -bsequence "$TMPDIR/doubled_fwd.fa" \
          -gapopen 10 -gapextend 0.5 -stdout -auto 2>/dev/null > "$TMPDIR/water_fwd.out"

    # Align anchor vs RC doubled
    water -asequence "$TMPDIR/anchor.fa" -bsequence "$TMPDIR/doubled_rc.fa" \
          -gapopen 10 -gapextend 0.5 -stdout -auto 2>/dev/null > "$TMPDIR/water_rc.out"

    SCORE_FWD=$(awk '/^# Score:/ {print $3}' "$TMPDIR/water_fwd.out" || echo "0")
    SCORE_RC=$(awk  '/^# Score:/ {print $3}' "$TMPDIR/water_rc.out"  || echo "0")
    [[ -z "$SCORE_FWD" ]] && SCORE_FWD="0"
    [[ -z "$SCORE_RC"  ]] && SCORE_RC="0"

    SCORE_FWD_INT=$(echo "$SCORE_FWD" | awk '{printf "%.0f", $1 * 10}')
    SCORE_RC_INT=$( echo "$SCORE_RC"  | awk '{printf "%.0f", $1 * 10}')

    if [[ $SCORE_RC_INT -gt $SCORE_FWD_INT ]]; then
        WATER_OUT="$TMPDIR/water_rc.out"
        DOUBLED="$DOUBLED_RC"
        ORIENTATION="RC"
        BEST_SCORE="$SCORE_RC"
    else
        WATER_OUT="$TMPDIR/water_fwd.out"
        DOUBLED="$DOUBLED_FWD"
        ORIENTATION="FWD"
        BEST_SCORE="$SCORE_FWD"
    fi

    # Parse subject start position from water output
    START=$(awk '/^doubled/ {for(i=1;i<=NF;i++) if($i~/^[0-9]+$/) {print $i; exit}}' "$WATER_OUT")
    if [[ -z "$START" || "$START" == "0" ]]; then
        START=$(grep -A1 "^doubled" "$WATER_OUT" | head -1 | awk '{print $2}')
    fi
    if [[ -z "$START" || "$START" == "0" ]]; then
        echo "warning: no alignment found, keeping original (ungapped)"
        echo "$header"      >> "$OUTFA"
        echo "$seq"         >> "$OUTFA"
        continue
    fi

    START_IDX=$((START - 1))
    (( START_IDX >= LEN )) && START_IDX=$((START_IDX - LEN))

    # Rotate from the degapped doubled sequence
    ROTATED="${DOUBLED:$START_IDX:$LEN}"

    echo "orient: $ORIENTATION, score: $BEST_SCORE, rotation: $START_IDX"

    echo "$header"  >> "$OUTFA"
    echo "$ROTATED" >> "$OUTFA"

done < "$TMPDIR/others.fa"

echo ""
echo "Output: $OUTFA  (all sequences ungapped and rotated)"
echo "Run 'mafft --localpair --maxiterate 1000 $OUTFA > ${OUTFA}.al' to align"
