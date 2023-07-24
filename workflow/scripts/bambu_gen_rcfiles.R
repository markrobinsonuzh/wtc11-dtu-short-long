log <- file(snakemake@log[[1]], open="wt")
sink(log)

library(bambu)

bam_paths <- snakemake@input[["bam_files"]]
fa_file <- snakemake@input[["reference"]]
gtf_file <- snakemake@input[["annotation"]]
out_dir <- snakemake@output[["out"]]

if (!dir.exists(out_dir)){
  dir.create(out_dir,recursive=TRUE)
}

if (grepl("_0%",gtf_file)){
  se <- bambu(reads = bam_paths, rcOutDir = out_dir, genome = fa_file, discovery = FALSE, quant = FALSE, ncore = 4)
} else {
  bambuAnnotations <- prepareAnnotations(gtf_file)
  se <- bambu(reads = bam_paths, rcOutDir = out_dir, annotations = bambuAnnotations, genome = fa_file, discovery = FALSE, quant = FALSE, ncore = 4)
}