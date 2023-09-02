import sys

# logging
sys.stderr = open(snakemake.log[0], "w")
sys.stdout = open(snakemake.log[0], "w")

import pandas as pd
import time

start = time.time()

print(snakemake.config["samples"])

# Read annotations
annotation_df = pd.read_csv(snakemake.config["samples"], sep='\t', index_col=0)

print("Load anotations: {}".format(time.time() - start))

# Read transcription
tpm_df = pd.read_csv(snakemake.input.normcounts, sep='\t',  index_col=0)

print("Load tpm: {}".format(time.time() - start))

print(list(snakemake.config["diffexp"]["variables_of_interest"].keys()))

labels = annotation_df[list(snakemake.config["diffexp"]["variables_of_interest"].keys())]

tpm_df = pd.merge(tpm_df, labels, left_index=True, right_index=True)

print("Merge with Labels: {}".format(time.time() - start))


# Save file
tpm_df.to_csv(snakemake.output.parsedcounts, sep="\t")
