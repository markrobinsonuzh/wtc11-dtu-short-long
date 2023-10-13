library(tidyr)
library(dplyr)
library(readr)
library(tibble)
library(edgeR)
library(pheatmap)

assignments <- read.delim(snakemake@input[["assignments"]],sep='\t', header = TRUE)
counts <- read.delim(snakemake@input[["counts"]],sep='\t', header = TRUE)

col_order <- c("TXNAME","GENEID","day0.rep1.aligned", "day0.rep2.aligned", "day0.rep3.aligned",
               "day5.rep1.aligned", "day5.rep2.aligned", "day5.rep3.aligned")

counts <- counts[, col_order]

merged <- merge(x=counts,y=assignments, 
      by.x=c("TXNAME"), 
      by.y=c("transcript_id"))

for (i in 1:nrow(merged)){
    if (merged[i,"label"] == "TP" & startsWith(merged[i,"TXNAME"],"Bambu")) {
        merged[i,"TXNAME"] = merged[i,"ref_transcript"]
        merged[i,"GENEID"] = merged[i,"ref_gene"]
    }
}


trans_counts_table <- merged[,col_order]  %>% as.data.frame

write.table(trans_counts_table,snakemake@output[["modified_counts"]], sep='\t',quote=FALSE,row.names=FALSE)