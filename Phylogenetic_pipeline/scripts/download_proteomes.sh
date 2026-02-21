#!/usr/bin/env bash
set -euo pipefail

INPUT="configs/test_species_list.txt" #"$1"
THREADS=4 #"$2"
OUTDIR=intermediate_results/test/proteomes #"$3"
ALL_PROTEINS=intermediate_results/test/all_proteins.faa #"$4"


mkdir -p "$OUTDIR"

while read -r GENOME; do
    [[ -z "$GENOME" ]] && continue
    echo ">>> Downloading: $GENOME"

    SAFE_NAME=$(echo "$GENOME" | tr ' /' '__')
    WORKDIR="$OUTDIR/$SAFE_NAME"

    mkdir -p "$WORKDIR"
    cd "$WORKDIR"

    datasets download genome taxon "$GENOME" \
        --include protein \
        --reference \
        --filename dataset.zip

    unzip -q dataset.zip
    rm dataset.zip

    # Collect proteome
    find ncbi_dataset -name "*.faa" -exec cat {} \; > proteome.faa || true

    if [[ ! -s proteome.faa ]]; then
        echo "WARNING: No proteome found for $GENOME"
    else
        echo "Proteome saved: $WORKDIR/proteome.faa"
    fi

    cd - >/dev/null

done < "$INPUT"

find "$OUTDIR" -name "*.faa" -exec cat {} + > "$ALL_PROTEINS"


echo "All downloads finished"
