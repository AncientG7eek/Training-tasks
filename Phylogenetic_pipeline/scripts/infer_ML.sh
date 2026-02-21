#!/bin/bash
set -e

ALN_DIR="$1"
TREE_DIR="$2"
THREADS="$3"
LOGFILE="$4"
BOOTSTRAP="$5"

mkdir -p "$TREE_DIR"

exec > "$LOGFILE" 2>&1

echo "==== IQ-TREE parallel pipeline started ===="
START_TOTAL=$(date +%s)

export TREE_DIR

parallel -j "$THREADS" '
    aln={}
    base=$(basename "$aln" _aln.faa)
    echo "Inferring tree for $base"

    iqtree2 \
        -s "$aln" \
        -m MFP \
        -bb "$BOOTSTRAP" \
        -alrt 1000 \
        -nt 1 \
        -pre "$TREE_DIR/$base"
' ::: "$ALN_DIR"/*_aln.faa

END_TOTAL=$(date +%s)
echo "==== IQ-TREE pipeline finished ===="
echo "Total time: $((END_TOTAL-START_TOTAL)) seconds"

