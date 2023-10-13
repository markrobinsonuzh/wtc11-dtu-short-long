library(GenomicRanges, quietly = TRUE)

relevant_cols <- c("seqnames","start","end","width","strand","source","type","score","phase","gene_id",
                   "transcript_id","exon_assignment")

gencode_gtf_filepath <- snakemake@input[["GENCODE_annotation"]]

gencode_gtf <- rtracklayer::import(gencode_gtf_filepath)
gencode_gtf <- as.data.frame(gencode_gtf)
gencode_gtf <- subset(gencode_gtf, select= -c(exon_assignment))

merged_gtf_granges <- makeGRangesFromDataFrame(gencode_gtf, keep.extra.columns = TRUE)
rtracklayer::export(merged_gtf_granges, snakemake@output[["out"]])
