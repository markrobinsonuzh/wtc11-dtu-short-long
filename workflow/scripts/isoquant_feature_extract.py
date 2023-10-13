import numpy as np
import pandas as pd
from pyranges import read_gtf
import dask.dataframe as dd
import time
import os
import warnings

warnings.filterwarnings("ignore")

# from dask.distributed import Client
# client = Client()

# os.chdir(r"/home/asamant/pacbio_benchmark/isoquant/unfiltered/all_samples_25%_anno/OUT")

def perc_strand(x):
    perc_pos = sum(x == "+")/len(x)
    if perc_pos > 0.5:
        return perc_pos
    else:
        return 1 - perc_pos


extended_anno_file =  snakemake.input["extended_anno"]
corrected_reads_file =  snakemake.input["corrected_reads"]
gene_counts_file = snakemake.input["gene_counts"]
transcript_counts_file = snakemake.input["transcript_counts"]
model_reads_file = snakemake.input["model_reads"]
converted_fasta =  snakemake.input["converted_fasta"]

feature_table = dd.from_pandas(pd.DataFrame(), npartitions=1).reset_index(drop=True)

extended_anno = dd.from_pandas(read_gtf(extended_anno_file).df, npartitions = 1).reset_index(drop=True)
extended_anno = extended_anno[extended_anno["Feature"] == "transcript"].reset_index(drop=True)

transcript_counts = dd.read_csv(transcript_counts_file, delimiter="\t").rename(columns={"#feature_id":"transcript_id"})
    
linked_transcript_counts = extended_anno.merge(transcript_counts, on="transcript_id")

gene_counts = dd.from_pandas(pd.DataFrame(), npartitions=1).reset_index(drop=True)

gene_counts = linked_transcript_counts.groupby("gene_id").sum().reset_index()

linked_transcript_counts = linked_transcript_counts.merge(gene_counts, on="gene_id",suffixes=("_transcript_counts","_gene_counts"))

linked_transcript_counts["gene_prop"] = linked_transcript_counts["transcript_counts"]/linked_transcript_counts["gene_counts"]
linked_transcript_counts["gene_prop"] = linked_transcript_counts["gene_prop"].fillna(0)

feature_table = linked_transcript_counts[["transcript_id","gene_prop","tpm","Strand"]]

del(linked_transcript_counts, gene_counts, transcript_counts)

model_reads = dd.read_csv(model_reads_file,delimiter='\t', header=None)
model_reads.columns = ["name","transcript_id"]

corrected_reads = dd.read_csv(corrected_reads_file, delimiter = '\t')

merged = dd.merge(model_reads, corrected_reads[["chromStart","chromEnd","strand","name"]], on='name')

merged["5_end_sd"] = merged.map_partitions(lambda df: df.apply(lambda x: x["chromStart"] if x["strand"] == "+" else x["chromEnd"], axis = 1))
merged["3_end_sd"] = merged.map_partitions(lambda df: df.apply(lambda x: x["chromEnd"] if x["strand"] == "+" else x["chromStart"], axis = 1))

per_transcript_sd = merged[["transcript_id","5_end_sd","3_end_sd","strand"]].groupby("transcript_id").apply(lambda x: pd.Series({"5_end_sd":x["5_end_sd"].std(),
                                                                                                                        "3_end_sd":x["3_end_sd"].std(),
                                                                                                                        "perc_strand":perc_strand(x["strand"])})).reset_index()
per_transcript_sd = per_transcript_sd.fillna(0)
feature_table = feature_table.merge(per_transcript_sd, on="transcript_id")

converted_fasta = dd.read_csv(converted_fasta, delimiter='\t', header=None)
converted_fasta.columns = ["transcript_id","sequence"]

feature_table = feature_table.merge(converted_fasta, on="transcript_id")

print(feature_table.columns)

feature_table["freq_A_start_20"] = feature_table.map_partitions(lambda df: df.apply(lambda x: x["sequence"][:20].count("A") if x["Strand"] == "+" else x["sequence"][-20:].count("T"), axis = 1))
feature_table["freq_T_start_20"] = feature_table.map_partitions(lambda df: df.apply(lambda x: x["sequence"][:20].count("T") if x["Strand"] == "+" else x["sequence"][-20:].count("A"), axis = 1))
feature_table["freq_A_end_20"] = feature_table.map_partitions(lambda df: df.apply(lambda x: x["sequence"][-20:].count("A") if x["Strand"] == "+" else x["sequence"][:20].count("T"), axis = 1))
feature_table["freq_T_end_20"] = feature_table.map_partitions(lambda df: df.apply(lambda x: x["sequence"][-20:].count("T") if x["Strand"] == "+" else x["sequence"][:20].count("A"), axis = 1))


feature_table = feature_table.drop(["sequence","Strand"], axis = 1)
feature_table["label"] = feature_table.apply(lambda x: 1 if x["transcript_id"].startswith("ENST") or x["transcript_id"].startswith("SIRV") else 0, axis = 1)

feature_table = feature_table.compute()
feature_table.to_csv(snakemake.output["feature_table"],sep='\t',index=False)






