import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict
import os
from tqdm import tqdm

# ====== CONFIG ======
CLUSTERS_TSV = snakemake.input.clusters_tsv
SEARCH_RESULTS_TSV = snakemake.input.search_result
PROTEINS_FASTA = snakemake.inputs.all_protein_fasta
OUTPUT_DIR = snakemake.output.clusters_fasta
CHUNK_SIZE = snakemake.params.chunk_size  # dla dużego pliku TSV
# ====================

RELATIONS_PATH = "mmseqs/clusters.tsv"

relations = pd.read_csv(RELATIONS_PATH, sep='\t', header=None, names=["centroid", "cluster_member"])
relations["genomes"] = relations["cluster_member"].str.split("|").str[0]

clusters = relations.groupby("centroid")

def has_paralogs(df):
    return df["genomes"].nunique() < len(df)


clusters_with_paralogs = []
small_clusters_without_paralogs = []
clusters_1to1 = []
all_genomes_clusters = 0
for centroid, group in clusters:
    if group["genomes"].nunique() == 10:
        all_genomes_clusters += 1
    if has_paralogs(group):
        clusters_with_paralogs.append(group)
    elif len(group) == 10:
        clusters_1to1.append(group)
    else:
        small_clusters_without_paralogs.append(group)

def keep_best_paralog():

    # Wczytanie mapping {protein : centroid}
    print("Loading clusters mapping...")
    clusters_df = pd.read_csv(CLUSTERS_TSV, sep="\t", header=None, names=["centroid","protein"])
    protein_to_centroid = dict(zip(clusters_df["protein"], clusters_df["centroid"]))

    # Przygotowanie struktury do przechowywania najlepszego paraloga
    # best_hits[(centroid, genome)] = (protein, mean_bits)
    best_bitscores = defaultdict(list)

    # Chunkowe wczytywanie search_results.tsv
    print("Processing search results in chunks...")

    cols = [
        "query","target","fident","alnlen","mismatches","gaps",
        "qstart","qend","tstart","tend","evalue","bits" 
    ]
    import subprocess

    total_lines = int(
        subprocess.check_output(["wc", "-l", SEARCH_RESULTS_TSV]).split()[0]
    )
    total_chunks = total_lines // CHUNK_SIZE + 1

    for chunk in tqdm(pd.read_csv(SEARCH_RESULTS_TSV, sep="\t", header=None, names=cols, chunksize=CHUNK_SIZE), total= total_chunks):
        # wyciągnij genome z nagłówka
        chunk["query_genome"] = chunk["query"].str.split("|").str[0]
        chunk["target_genome"] = chunk["target"].str.split("|").str[0]



        # centroidy
        chunk["query_centroid"]  = chunk["query"].map(protein_to_centroid)
        chunk["target_centroid"] = chunk["target"].map(protein_to_centroid)

        chunk = chunk[
            (chunk["query_genome"] != chunk["target_genome"]) &
            (chunk["query_centroid"].notna()) &
            (chunk["query_centroid"] == chunk["target_centroid"])
        ]

        edges = chunk[["query_centroid", "query_genome", "query", "bits"]].rename(columns={"query_genome":"genome", "query_centroid":"centroid", "query":"protein"})
            




        # agregacja: paralog → średni bits
        grouped = (
            edges
            .groupby(["centroid","genome","protein"])["bits"]
            .mean()
            .reset_index()
        )

        # zapamiętaj
        for _, row in grouped.iterrows():
            key = (row["centroid"], row["genome"])
            best_bitscores[key].append((row["protein"], row["bits"]))
    # Wybór najlepszego paraloga dla każdej pary (centroid, genome)
    print("Selecting best paralogs per genome per cluster...")
    best_hits = {}
    for key, hits in best_bitscores.items():
        # sortowanie malejąco po bits i wybór największego
        best_protein = max(hits, key=lambda x: x[1])[0]
        best_hits[key] = best_protein

    # select 1:1 clusters
    all_genomes = set(
        clusters_df["protein"].str.split("|").str[0]
    )
    N_GENOMES = len(all_genomes)

    cluster_to_genomes = defaultdict(set)

    for (centroid, genome), protein in best_hits.items():
        cluster_to_genomes[centroid].add(genome)

    # from collections import Counter
    # print(f"Distribution of genomes number in clusters {Counter(len(g) for g in cluster_to_genomes.values())}", )
   
    one_to_one_clusters = {
        c for c, genomes in cluster_to_genomes.items()
        if len(genomes) == N_GENOMES
    }
    print(f"Found {len(one_to_one_clusters)} 1to1 clusters")

    # proteomy do pamięci (dict)
    print("Loading protein sequences...")
    protein_seqs = SeqIO.to_dict(
        SeqIO.parse(PROTEINS_FASTA, "fasta")
    )

    for centroid in one_to_one_clusters:
        records = []
        for genome in all_genomes:
            prot = best_hits[(centroid, genome)]
            records.append(
                SeqRecord(
                    protein_seqs[prot].seq,
                    id=genome,          # ← bardzo ważne!
                    description=""
                )
            )

        out_fasta = os.path.join(OUTPUT_DIR, f"{centroid}.faa")
        SeqIO.write(records, out_fasta, "fasta")

print("All done! FASTA files saved in:", OUTPUT_DIR)


print(f"Number of 1 to 1 clusters: {len(clusters_1to1)}")
print(f"Number of all genomes clusters in total (also these with paralogs): {all_genomes_clusters}")
print(f"Number of clusters without paralogs but not covering all genomes: {len(small_clusters_without_paralogs)}")
print(f"Number of clusters with paralogs: {len(clusters_with_paralogs)}")

keep_best_paralog()
