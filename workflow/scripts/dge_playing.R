## DGE analysis

## Load packages

library(tidyr)
library(dplyr)
library(readr)
library(tibble)
library(edgeR)
library(pheatmap)
library(GGally)


## Read in the transcript-level data and aggregate it

wd <- getwd()

trans_counts_file <- snakemake@input[["counts"]]

trans_counts <- read_tsv(trans_counts_file) %>% select(-GENEID) %>%
              column_to_rownames("TXNAME") %>% replace(is.na(.), 0) %>% as.data.frame 

if (!dir.exists(snakemake@output[["out"]])){
  dir.create(snakemake@output[["out"]])
}


## Make an `DGEList` object and do some filteri



# create edgeR container for the counts
d <- DGEList(counts=trans_counts)

# filter: CPM > 1 in at least 3 samples
cps <- cpm(d)
rs <- rowSums(cps > 1) >= 3
d <- d[rs,]

# normalization
d <- calcNormFactors(d)


# grouping variable (for design matrix)
grp <- sapply(colnames(d),function(x){strsplit(x,split=".",fixed=TRUE)[[1]][1]})

# design matrix
(mm <- model.matrix(~0+grp))

# book-keeping to include things in the table
d$genes <- data.frame(associated_transcript = rownames(trans_counts)[rs])
cps <- cpm(d)
colnames(cps) <- paste0("s", colnames(cps))
d$genes <- data.frame(d$genes, round(cpm(d), 2))
head(d$genes)

# estimate dispersion
d <- estimateDisp(d,mm)

# quick plot for diagnosis
plotBCV(d)
plotMDS(d, cex = 2)

# fit GLM; test for diff b/w day 5 and day 0
f <- glmQLFit(d,mm)
mc <- makeContrasts(grpday5-grpday0, levels=colnames(mm))
qlf <- glmQLFTest(f,contrast = mc)

# quick look at top genes
topTags(qlf)

# pull out table of statistics
tt <- topTags(qlf, n = Inf, sort.by = "none")

# crude filter on number of DE genes
w <- tt$table$FDR < .05 & abs(tt$table$logFC) > 1
table(w)

# P-value histogram
hist(tt$table$PValue, 200)

# MA plot
plotMD(qlf, ylim=c(-5,5))

# spot check on hypervariable genes
# w <- tt$table$logFC > 4 & d$tagwise.dispersion > 1.5 & tt$table$logCPM > 6 & tt$table$FDR > .2
# table(w)
# cpm(d)[w,][,c(1:3,9:11)]

table(w <- tt$table$FDR < .005 & abs(tt$table$logFC) > 3)

setwd(snakemake@output[["out"]])

pdf("quick_heatmap.pdf", height=20, width=8)
pheatmap(log2(1+cps[w,c(1:3,4:6)]), 
         fontsize_row = 8)
dev.off()

pdf("ggpairs.pdf", width = 8, height = 8)
ggpairs( as.data.frame(log2(cps[,c(1:3,4:6)]+1)) )
dev.off()

png("ggpairs_t.png", width = 800, height = 800)
ggpairs( as.data.frame(log2(1+trans_counts[,c(1:3,4:6)])) )
dev.off()

setwd(wd)
# write statistics to table
write.table(tt$table, file.path(snakemake@output[["out"]],"08062023_DGE_gloria_pacbio.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

saveRDS(tt, file.path(snakemake@output[["out"]],"pacbio_edgeR-DGE-topTags-table.rds"))

cor(log2(cps[,c(1:3,4:6)]+1), method = "spearman")

cor(log2(cps[,c(1:3,4:6)]+1), method = "pearson")

