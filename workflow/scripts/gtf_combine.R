
library(GenomicRanges, quietly = TRUE)

relevant_cols <- c("seqnames","start","end","width","strand","source","type","score","phase","gene_id",
                   "transcript_id")

gencode_gtf_filepath <- snakemake@input[["GENCODE_annotation"]]

gencode_gtf <- rtracklayer::import(gencode_gtf_filepath)
gencode_gtf <- as.data.frame(gencode_gtf)
gencode_gtf <- gencode_gtf[,relevant_cols]

if (snakemake@wildcards[["perc"]] != "0%"){

    SIRV_gtf_filepath <- snakemake@input[["SIRV_annotation"]]

    SIRV_gtf <- rtracklayer::import(SIRV_gtf_filepath)
    SIRV_gtf <- as.data.frame(SIRV_gtf)
    SIRV_gtf <- SIRV_gtf[,relevant_cols]

    merged_gtf <- rbind(gencode_gtf, SIRV_gtf)

    merged_gtf_granges <- makeGRangesFromDataFrame(merged_gtf, keep.extra.columns = TRUE)
    rtracklayer::export(merged_gtf_granges, snakemake@output[["out"]])
} else {
    merged_gtf_granges <- makeGRangesFromDataFrame(gencode_gtf, keep.extra.columns = TRUE)
    rtracklayer::export(merged_gtf_granges, snakemake@output[["out"]])
}