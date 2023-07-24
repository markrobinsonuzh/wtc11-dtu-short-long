# PARAMS:
# args[1] - input file
# args[2] - output folder

library(GenomicRanges, quietly = TRUE)


args <- commandArgs(trailingOnly = TRUE)

temp <- strsplit(args[1],split = ".",fixed = TRUE)[[1]]
gtf_file <- temp[1]

sampling_percs <- seq(1,0,-0.1)

full_gtf <- rtracklayer::import(sprintf("%s.gtf",gtf_file))
full_gtf <- as.data.frame(full_gtf)

sampled_transcripts <- unique(full_gtf[,"transcript_id"])
full_anno_size <- length(sampled_transcripts)
anno_size <- length(sampled_transcripts)

for (perc in sampling_percs){
  
  anno_sample <- sample(1:anno_size, ceiling(perc*full_anno_size), replace = FALSE)
  
  sampled_transcripts <- sampled_transcripts[anno_sample]
  anno_size <- length(sampled_transcripts)

  sampled_gtf <- full_gtf[full_gtf[,"transcript_id"] %in% sampled_transcripts, ]
  sampled_gtf <- makeGRangesFromDataFrame(sampled_gtf, keep.extra.columns = TRUE)
  if (perc > 0){
  rtracklayer::export(sampled_gtf, sprintf("%s_sampled_%d%%.gtf",gtf_file,as.integer(ceiling(perc*100))))
  } else {
    dummy_gtf <- data.table::data.table(va=numeric(), vb=numeric(), vc=numeric())
    data.table::fwrite(dummy_gtf, sprintf("%s_sampled_%d%%.gtf",gtf_file,as.integer(ceiling(perc*100))))
  }
}


