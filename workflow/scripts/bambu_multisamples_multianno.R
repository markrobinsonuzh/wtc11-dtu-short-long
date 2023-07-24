log <- file(snakemake@log[[1]], open="wt")
sink(log)

library(bambu, quietly = TRUE)

library(BiocFileCache)
bfc <- BiocFileCache(snakemake@input[["rc_files"]], ask = FALSE)
info <- bfcinfo(bfc)

gtf_file <- snakemake@input[["annotation"]]
fa_file <- snakemake@input[["reference"]]
  
if (grepl("_0%",gtf_file)){
  se <- bambu(reads = info$rpath, genome = fa_file, NDR = 1)
  writeBambuOutput(se, path = snakemake@output[["out"]])
} else {
    bambuAnnotations <- prepareAnnotations(gtf_file)
    se <- bambu(reads = info$rpath, annotations = bambuAnnotations, genome = fa_file, NDR = 1)
    writeBambuOutput(se, path = snakemake@output[["out"]])
}


