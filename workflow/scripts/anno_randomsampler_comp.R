library(GenomicRanges, quietly = TRUE)

full_gtf_filepath <- snakemake@input[["full_annotation"]]
perc_value <- snakemake@wildcards[["perc"]]

full_gtf <- rtracklayer::import(full_gtf_filepath)
full_gtf <- as.data.frame(full_gtf)

temp <- strsplit(full_gtf_filepath,split = ".",fixed = TRUE)[[1]]
full_gtf_filepath <- temp[1]

if (perc_value == "0%"){

    rtracklayer::export(full_gtf, sprintf("%s_sampled_0%%_complement.gtf",full_gtf_filepath))

} else if (perc_value == "100%") {

    dummy_gtf <- data.table::data.table(dummy1=numeric(), dummy2=numeric(), dummy3=numeric())
    data.table::fwrite(dummy_gtf, sprintf("%s_sampled_100%%_complement.gtf",full_gtf_filepath))

} else {

    sampled_gtf <- snakemake@input[["sampled_annotation"]]
    sampled_gtf <- rtracklayer::import(sampled_gtf)
    sampled_gtf <- as.data.frame(sampled_gtf)

    sampled_transcripts <- unique(sampled_gtf[,"transcript_id"])

    comp_gtf <- full_gtf[!(full_gtf[,"transcript_id"] %in% sampled_transcripts), ]
    comp_gtf <- makeGRangesFromDataFrame(comp_gtf, keep.extra.columns = TRUE)

    rtracklayer::export(comp_gtf, sprintf("%s_sampled_%s_complement.gtf",full_gtf_filepath,perc_value))

}