#!/bin/bash
set -e  # zatrzymaj przy pierwszym błędzie

INPUT_DIR="$1"
OUTPUT_DIR="$2"
THREADS="$3"  # liczba wątków dla MAFFT
LOGFILE="$4"

exec > "$LOGFILE" 2>&1

mkdir -p "$OUTPUT_DIR"

echo "==== MAFFT pipeline started ===="
START_TOTAL=$(date +%s)

# iteracja po wszystkich plikach FASTA
for f in "$INPUT_DIR"/*.faa; do
    base=$(basename "$f" .faa)
    echo "Aligning $base..."
    START=$(date +%s)

    mafft --auto --anysymbol --thread "$THREADS" "$f" > "$OUTPUT_DIR/${base}_aln.faa"

    END=$(date +%s)
    echo "Alignment for $base took $((END-START)) seconds"
done

END_TOTAL=$(date +%s)
echo "==== MAFFT pipeline finished ===="
echo "Total time: $((END_TOTAL-START_TOTAL)) seconds"
