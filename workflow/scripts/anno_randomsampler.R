# PARAMS:
# args[1] - input file
# args[2] - output folder

set.seed(1234)

library(GenomicRanges, quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)

temp <- strsplit(args[1],split = ".",fixed = TRUE)[[1]]
gtf_file <- temp[1]

sampling_percs <- seq(1,0,-0.25)

full_gtf <- rtracklayer::import(sprintf("%s.gtf",gtf_file))
full_gtf <- as.data.frame(full_gtf)

chromosomes <- unique(full_gtf[,"seqnames"])

trans_per_chrom <- list()

for (chrom in chromosomes){
  chrom_trans <- full_gtf[full_gtf[,"seqnames"] == chrom, ]
  trans_per_chrom[[chrom]] <- length(unique(chrom_trans[,"transcript_id"]))
} 

sampled_gtf <- full_gtf

for (perc in sampling_percs){
  
  sampled_transcripts <- c()
  
  for (chrom in chromosomes) {
    
    chrom_trans <- sampled_gtf[sampled_gtf[,"seqnames"] == chrom, ]
    
    chrom_transcripts <- unique(chrom_trans[,"transcript_id"])
    anno_size <- length(chrom_transcripts)
  
    anno_sample <- sample(1:anno_size, ceiling(perc*trans_per_chrom[[chrom]]), replace = FALSE)
    sampled_transcripts <- append(sampled_transcripts, chrom_transcripts[anno_sample])
        
  }
  
  sampled_gtf <- full_gtf[full_gtf[,"transcript_id"] %in% sampled_transcripts, ]
  sampled_gtf_granges <- makeGRangesFromDataFrame(sampled_gtf, keep.extra.columns = TRUE)
  
  if (perc > 0){
  rtracklayer::export(sampled_gtf_granges, sprintf("%s_sampled_%d%%.gtf",gtf_file,as.integer(ceiling(perc*100))))
  } else {
    dummy_gtf <- data.table::data.table(va=numeric(), vb=numeric(), vc=numeric())
    data.table::fwrite(dummy_gtf, sprintf("%s_sampled_%d%%.gtf",gtf_file,as.integer(ceiling(perc*100))))
  }
}


