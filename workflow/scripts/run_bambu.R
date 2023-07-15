run_bambu <- function(reads, annotations, genome, discovery, output_directory, ncore, lowMemory = TRUE) {
    library(bambu)
    bambuAnnotations <- prepareAnnotations(annotations)
    bambu_result <- bambu::bambu(reads = reads, annotations = bambuAnnotations, genome = genome, discovery = discovery,
                                 ncore = ncore, lowMemory = lowMemory
    )
    bambu::writeBambuOutput(bambu_result, path = output_directory)
}

run_bambu(
     reads = unlist(strsplit(snakemake@input[["reads"]], ",")),
     annotations = snakemake@input[["annotation"]],
     genome = snakemake@input[["reference"]],
     discovery = FALSE,
     output_directory = snakemake@output,
     ncore = snakemake@threads,
     lowMemory = FALSE
     )
