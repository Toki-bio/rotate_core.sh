# rotate_core.sh - Documentation

## Purpose
Rotate circular satellite DNA sequences to align at a conserved core region 
using a reference anchor sequence.

## Problem
Satellite DNA monomers are circular - same sequence can start at any position:
```
Original:  ABCDEFGHIJ  (starts at A)
Rotation:  FGHIJABCDE  (same sequence, starts at F)
```
When comparing sequences from different genomic locations, they may be 
in different rotational phases, making alignment difficult.

## Algorithm

### Pseudocode
```
INPUT: sequences.fa, anchor_name

1. PREPROCESS
   - Remove gaps, linearize FASTA
   - Extract anchor sequence

2. PASS 1: COLLECT ALIGNMENT DATA
   For each sequence (except anchor):
       doubled = seq + seq                    # Handle circular wrap
       hit = ssearch36(anchor, doubled)       # Find best local alignment
       Store: qstart, qend, sstart, send, strand, score
       
       # Strand detection: if qstart > qend → reverse complement needed

3. PASS 2: DETERMINE ANCHOR OFFSET  
   best_hit = hit with highest bitscore
   if best_hit.strand == "-":
       anchor_offset = best_hit.qend - 1      # Use qend for reverse
   else:
       anchor_offset = best_hit.qstart - 1    # Use qstart for forward
   
   Rotate anchor by anchor_offset

4. PASS 3: ROTATE ALL SEQUENCES
   For each sequence:
       if strand == "-":
           seq = reverse_complement(seq)
           start = LEN - 1 - ((send - 1) % LEN) + anchor_offset
       else:
           start = sstart - qstart + anchor_offset
       
       start = start % LEN                    # Normalize
       rotated = seq[start:] + seq[:start]    # Cut and rejoin

OUTPUT: sequences.rotated.fa
```

## Toy Example

### Input
```
>Anchor
AAACTTGTTGGTGTGTTTT    (19bp, core at position 1)

>Seq1  
TGTTTTAAACTTGTTGGTG    (same seq, rotated by 13)

>Seq2_revcomp
AAACACACCCAACAAGTTT    (reverse complement, needs RC + rotation)
```

### Step-by-step

**PASS 1: Align Anchor vs Doubled sequences**

Seq1 doubled: `TGTTTTAAACTTGTTGGTGTGTTTTAAACTTGTTGGTG` (38bp)

ssearch36 finds: Anchor[1-19] matches Doubled[7-25]
```
qstart=1, qend=19, sstart=7, send=25, strand=+
```

**PASS 2: Anchor offset**
```
Best hit: qstart=1, strand=+
anchor_offset = 1 - 1 = 0  (no anchor rotation)
```

**PASS 3: Rotate Seq1**
```
strand = "+"
start = sstart - qstart + anchor_offset
start = 7 - 1 + 0 = 6

Seq1 = TGTTTTAAACTTGTTGGTG
       0123456789...
              ^cut here (position 6)

Rotated = AAACTTGTTGGTGTGTTTT  ✓ (matches anchor)
```

## Key Concepts

### Why double the sequence?
```
Original (10bp): FGHIJABCDE
Looking for:     ABCDE (at position 6)

Problem: Can't extract ABCDE...FGH in one substring

Solution - Double it: FGHIJABCDEFGHIJABCDE
Now ABCDEFGHIJ is contiguous at position 6-15
```

### Strand detection (ssearch36 m8C format)
```
Forward hit:  qstart=1,  qend=146  (qstart < qend)
Reverse hit:  qstart=146, qend=1   (qstart > qend)
```

### Reverse complement rotation
When strand="-", the sequence matched in reverse orientation:
1. Reverse complement the sequence first
2. Calculate rotation point based on `send` (not sstart)
3. Adjust for coordinate transformation

## Limitations

1. **Divergent sequences**: If sequence is too different from anchor,
   ssearch36 may find suboptimal alignment region

2. **Core vs flanking**: Script aligns at ssearch36's "best" match,
   which may not be the biologically meaningful core region

3. **Short anchors**: Very short anchors may have multiple equally 
   good matches, leading to inconsistent rotation

## Usage
```bash
./rotate_core.sh input.fa anchor_name
# Output: input.rotated.fa

# Then align:
mafft --localpair --maxiterate 1000 input.rotated.fa > aligned.fa
```

## Files
- rotate_core_v13.sh - Clean production version
- rotate_core_debug.sh - Verbose debug version (keeps temp files)
