run_bambu <- function(reads, annotations, genome, discovery, output_directory, ncore, lowMemory = TRUE) {
    library(bambu)
    bambuAnnotations <- bambu::prepareAnnotations(annotations)
    bambu_result <- bambu::bambu(reads = reads, annotations = bambuAnnotations, genome = genome, discovery = discovery,
                                 ncore = ncore, lowMemory = lowMemory
    )
    bambu::writeBambuOutput(bambu_result, path = output_directory)
}

run_bambu(
     reads = unlist(strsplit(snakemake@input[["reads"]], ",")),
     annotations = snakemake@input[["annotation"]],
     genome = snakemake@input[["reference"]],
     discovery = snakemake@params[["discovery"]],
     output_directory = snakemake@params[["output_dir"]],
     ncore = snakemake@threads,
     lowMemory = TRUE
     )
