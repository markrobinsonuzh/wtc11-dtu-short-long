library(GenomicRanges, quietly = TRUE)

get_transcript_id <- function(x) {
    if (x != "-"){
    stringr::str_split(x,stringr::fixed("|"))[[1]][2]
    } else {
        "-"
    }
}

get_gene_id <- function(x) {
    if (x != "-"){
    stringr::str_split(x,stringr::fixed("|"))[[1]][1]
    } else {
        "-"
    }
}

total_annotations = 69
perc_anno = as.numeric(stringr::str_split(snakemake@wildcards[["perc"]],stringr::fixed("%"))[[1]][1])
retained_anno = ceiling(perc_anno * total_annotations / 100)
removed_anno = total_annotations - retained_anno

gffcmp_retained <- read.delim(paste(snakemake@input[["retained_gff"]],"gffcmp.tracking",sep='/'),sep='\t', header = FALSE)
names(gffcmp_retained) <- c("Retained_Query_Transfrag_ID","Retained_Query_Locus_ID","Retained_Reference_Gene_ID","Retained_Class_Code","V5")
gffcmp_retained$transcript_id <- sapply(gffcmp_retained[,"V5"],get_transcript_id)
gffcmp_retained$retained_ref_transcript_id <- sapply(gffcmp_retained[,"Retained_Reference_Gene_ID"],get_transcript_id)
gffcmp_retained$retained_ref_gene_id <- sapply(gffcmp_retained[,"Retained_Reference_Gene_ID"],get_gene_id)

gffcmp_removed <- read.delim(paste(snakemake@input[["removed_gff"]],"gffcmp.tracking",sep='/'),sep='\t', header = FALSE)
names(gffcmp_removed) <- c("Removed_Query_Transfrag_ID","Removed_Query_Locus_ID","Removed_Reference_Gene_ID","Removed_Class_Code")
gffcmp_removed$transcript_id <- sapply(gffcmp_removed[,5],get_transcript_id)
gffcmp_removed$removed_ref_transcript_id <- sapply(gffcmp_removed[,"Removed_Reference_Gene_ID"],get_transcript_id)
gffcmp_removed$removed_ref_gene_id <- sapply(gffcmp_removed[,"Removed_Reference_Gene_ID"],get_gene_id)


gffcmp_full <- read.delim(paste(snakemake@input[["full_gff"]],"gffcmp.tracking",sep='/'),sep='\t', header = FALSE)
names(gffcmp_full) <- c("Full_Query_Transfrag_ID","Full_Query_Locus_ID","Full_Reference_Gene_ID","Full_Class_Code")
gffcmp_full$transcript_id <- sapply(gffcmp_full[,5],get_transcript_id)
gffcmp_full$full_ref_transcript_id <- sapply(gffcmp_full[,"Full_Reference_Gene_ID"],get_transcript_id)
gffcmp_full$full_ref_gene_id <- sapply(gffcmp_full[,"Full_Reference_Gene_ID"],get_gene_id)


relevant_columns <- c("transcript_id","retained_ref_gene_id","retained_ref_transcript_id","Retained_Class_Code","removed_ref_gene_id",
                      "removed_ref_transcript_id","Removed_Class_Code","full_ref_gene_id","full_ref_transcript_id","Full_Class_Code")

merged <- merge(gffcmp_retained, gffcmp_removed)
merged <- merge(merged, gffcmp_full)
merged <- merged[,relevant_columns]

write.table(merged,snakemake@output[["merged_gff"]], sep='\t',quote=FALSE, row.names = FALSE)

TP_retained <- 0
TP_removed <- 0
FP_retained <- 0
FP_removed <- 0
ambiguous <- 0

transcript_ids <- merged$transcript_id
annotations <- c()
labels <- c()
refseq <- c()
refseqgene <- c()

print(head(merged))

for (i in 1:nrow(merged)){
    if (merged[i,"Retained_Class_Code"] == merged[i,"Full_Class_Code"] & merged[i,"Removed_Class_Code"] != merged[i,"Full_Class_Code"]) {
        annotations = append(annotations, "Retained")
        refseq <- append(refseq, merged[i,"retained_ref_transcript_id"])
        refseqgene <- append(refseqgene, merged[i,"retained_ref_gene_id"])
        if (merged[i,"Retained_Class_Code"] == "=") {
            TP_retained <- TP_retained + 1
            labels <- append(labels, "TP")
        } else {
            FP_retained <- FP_retained + 1
            labels <- append(labels, "FP")
        }
    } else if (merged[i,"Retained_Class_Code"] != merged[i,"Full_Class_Code"] & merged[i,"Removed_Class_Code"] == merged[i,"Full_Class_Code"]) {
       annotations = append(annotations, "Removed")
       refseq <- append(refseq, merged[i,"removed_ref_transcript_id"])
       refseqgene <- append(refseqgene, merged[i,"removed_ref_gene_id"])
       if (merged[i,"Removed_Class_Code"] == "=") {
            TP_removed <- TP_removed + 1
            labels <- append(labels, "TP")
        } else {
            FP_removed <- FP_removed + 1
            labels <- append(labels, "FP")
        } 
    } else {
        ambiguous <- ambiguous + 1
        annotations <- append(annotations, "Ambiguous")
        refseq <- append(refseq, paste(merged[i,"removed_ref_transcript_id"],merged[i,"retained_ref_transcript_id"],sep="/"))
        refseqgene <- append(refseqgene, paste(merged[i,"removed_ref_gene_id"],merged[i,"retained_ref_gene_id"],sep="/"))
        labels <- append(labels, "Ambiguous")
    }
}

precision_retained <- TP_retained / (TP_retained + FP_retained)
precision_removed <- TP_removed / (TP_removed + FP_removed)

recall_retained <- TP_retained / retained_anno
recall_removed <- TP_removed / removed_anno


assignment_table <- data.frame(transcript_id = transcript_ids, annotation = annotations, label = labels, ref_transcript = refseq,
                                ref_gene = refseqgene)
stats_table <- data.frame(recall_retained = c(recall_retained), recall_removed = c(recall_removed),
                         precision_retained = c(precision_retained), precision_removed = c(precision_removed),
                         TP_retained = c(TP_retained), FP_retained = c(FP_retained), TP_removed = c(TP_removed),
                        FP_removed = c(FP_removed), Ambiguous = c(ambiguous))

write.table(assignment_table,snakemake@output[["assignments"]], sep='\t',quote=FALSE, row.names = FALSE)
write.table(stats_table,snakemake@output[["stats"]], sep='\t',quote=FALSE, row.names = FALSE)




