library(rtracklayer)
library(GenomicRanges)
library(GenomicAlignments)
library(data.table)

setwd("/home/asamant/pacbio_benchmark/isoquant/unfiltered_per_sample/all_samples_75%_anno/day0-rep1-sirvs/OUT/")

extended_anno <- import("OUT.transcript_models.gtf")
extended_anno <- extended_anno[extended_anno$type == "exon"]
extended_anno <- as.data.frame(extended_anno)
extended_anno <- setDT(extended_anno)

extended_anno <- GRangesList(split(extended_anno,extended_anno$transcript_id))
extended_anno <- sort(extended_anno)

bamFile <- Rsamtools::BamFile("/home/Shared/data/seq/sheynkman_pb_mas_wtc11/data-raw/aligned/day0-rep1-sirvs.aligned.sorted.bam")
bf <- open(bamFile)
readGrgList <- list()
counter <- 1
while (isIncomplete(bf)) {
    readGrgList[[counter]] <-
        grglist(readGAlignments(bf,
        use.names = TRUE))
    counter <- counter + 1
}
on.exit(close(bf))

readGrgList <- readGrgList[[1]]
overlaps <- findOverlaps(readGrgList, extended_anno, type="equal", maxgap = 1000)

bed_file <- import("test.bed")
bed_file <- as.data.frame(bed_file)
bed_file <- setDT(bed_file)

isassigned = c()
for (i in 1:length(overlaps)){

    read <- readGrgList[[queryHits(overlaps)[i]]]
    anno <- extended_anno[[subjectHits(overlaps)[i]]]
    if (length(anno) == length(read)){
        x <- pintersect(anno, read)
        percentOverlap <- sum(width(read)) / sum(width(anno))
        isassigned = append(isassigned,percentOverlap > 0.95)
    }
    else {
        percentOverlap <- 1
        isassigned = append(isassigned,FALSE)
    }

}




