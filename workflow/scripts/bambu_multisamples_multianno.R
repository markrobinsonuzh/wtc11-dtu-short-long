log <- file(snakemake@log[["log"]], open="wt")
sink(log,type="message")

library(bambu, quietly = TRUE)

library(BiocFileCache, quietly = TRUE)
bfc <- BiocFileCache(snakemake@input[["rc_files"]], ask = FALSE)
info <- bfcinfo(bfc)

NDR_thresh <- 0


gtf_file <- snakemake@input[["annotation"]]
fa_file <- snakemake@input[["reference"]]

threshold <- snakemake@wildcards[["threshold"]]
  
if (grepl("_0%",gtf_file)){
  se <- bambu(reads = info$rpath, genome = fa_file, NDR = 1, ncore = snakemake@threads)
  writeBambuOutput(se, path = snakemake@output[["out"]])
} else {
    bambuAnnotations <- prepareAnnotations(gtf_file)
    if (threshold == "default_threshold"){
        se <- bambu(reads = info$rpath, annotations = bambuAnnotations, genome = fa_file, ncore = snakemake@threads)
    } else {
        threshold <- as.integer(strsplit(threshold,"_")[[1]][1])/100
        se <- bambu(reads = info$rpath, annotations = bambuAnnotations, genome = fa_file, NDR = threshold, ncore = snakemake@threads)
    }
    writeBambuOutput(se, path = snakemake@output[["out"]])
}

sink()
unlink(log)



