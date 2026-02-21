import re
import glob
import os

INPUT = snakemake.input.trees
OUTPUT = snakemake.output.trees_filtered
UFBOOT_THRESHOLD = snakemake.params.ufboot_thr
SHARLT_THRESHOLD = snakemake.params.sharlt_thr

def support_values(treefile):
    with open(treefile) as f:
        nwk = f.read()

    # match )X/Y:Z
    matches = re.findall(r"\)([\d.]+)/([\d.]+):", nwk)
    return [(float(shalrt), float(ufboot)) for shalrt, ufboot in matches]

def keep_tree(treefile):
    supports = support_values(treefile)
    if not supports:
        return False
    shalrt, ufboot = zip(*supports)
    return min(shalrt) >= SHARLT_THRESHOLD and min(ufboot) >= UFBOOT_THRESHOLD

treefiles = glob.glob(os.path.join(INPUT, "*.treefile"))


with open(OUTPUT, "w") as out:
    for tree in treefiles:
        if keep_tree(tree):
            with open(tree) as f:

                clean = re.sub("[0-9.:/]", "", f.read())
                out.write(clean.strip() + "\n")
