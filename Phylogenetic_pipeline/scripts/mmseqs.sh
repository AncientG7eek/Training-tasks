#!/bin/bash
set -e  # stop script if any command fails

# Measure total runtime
START_TOTAL=$(date +%s)

mkdir -p mmseqs 
cd mmseqs

# Variables
FASTA="$1"
DB="$2"
RESULT="$3"
TMP="$4"
CLUSTERS="$5"
THREADS="$6"
CONVERTALIS="$7"
TSV="$8"

echo "==== MMseqs2 pipeline started ===="

# Create MMseqs database
echo "Creating database..."
START=$(date +%s)
mmseqs createdb "$FASTA" "$DB"
END=$(date +%s)
echo "Database creation took $((END-START)) seconds"

# All-vs-all search
echo "Running all-vs-all search..."
START=$(date +%s)
mmseqs search "$DB" "$DB" "$RESULT" "$TMP" --threads $THREADS -s 7.5 --min-seq-id 0.3 -c 0.5
END=$(date +%s)
echo "All-vs-all search took $((END-START)) seconds"

mmseqs convertalis proteinsDB proteinsDB result "$CONVERTALIS"

# Clustering
echo "Clustering sequences..."
START=$(date +%s)
mmseqs cluster "$DB" "$CLUSTERS" "$TMP" --min-seq-id 0.3 -c 0.5 --cov-mode 0 --threads $THREADS
END=$(date +%s)
echo "Clustering took $((END-START)) seconds"

# Export clusters to TSV
echo "Exporting clusters to TSV..."
START=$(date +%s)
mmseqs createtsv "$DB" "$DB" "$CLUSTERS" "$TSV"
END=$(date +%s)
echo "Export to TSV took $((END-START)) seconds"

# Total runtime
END_TOTAL=$(date +%s)
echo "==== MMseqs2 pipeline finished ===="
echo "Total time: $((END_TOTAL-START_TOTAL)) seconds"
