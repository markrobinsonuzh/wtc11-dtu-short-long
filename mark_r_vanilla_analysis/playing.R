
library(rtracklayer)
library(Biostrings)
library(readr)
library(dplyr)
library(pheatmap)


fafs <- list.files("REF", ".fasta$", 
                   recursive = TRUE,  full.names = TRUE)

nap1 <- readDNAStringSet(fafs[1])
nap4 <- readDNAStringSet(fafs[3])

setd <- setdiff(names(nap4), names(nap1))
length(setd)

x <- import("day0-rep1/07_collapse/GS_230501_Day0-1.lima.0--0.lima.IsoSeqX_bc10_5p--IsoSeqX_3p.collapsed_transcripts.gff")
length(unique(x$gene_id))
length(unique(x$transcript_id))



rls <- read_csv("day0-rep1/02_skera/GS_230501_Day0-1.lima.0--0.skera.read_lengths.csv")

read_lengths <- rls %>% select(zmw, original_read_length) %>% unique

m_read_lengths <- rls %>% group_by(zmw) %>% summarize(mean_m_read = mean(m_read_length),
                                                      n = n())

# x.sirv6 <- x[seqnames(x)=="SIRV6"]
# x.sirv6
# length(unique(x.sirv6$gene_id))



# y <- read_tsv("GS_230501_Day0-1.lima.0--0.lima.IsoSeqX_bc10_5p--IsoSeqX_3p.collapsed_transcripts.sorted.gff", 
#               skip=3, col_names=c("seqnames","source","type","start","end","dot1","dash","dot2","ids"))
# y.sirv6 <- y %>% filter(seqnames=="SIRV6")
# y.sirv6
# y.sirv6$ids
# strsplit(y.sirv6$ids, ";", fixed=TRUE)
# sapply(strsplit(y.sirv6$ids, ";", fixed=TRUE), .subset, 1)
# gids <- sapply(strsplit(y.sirv6$ids, ";", fixed=TRUE), .subset, 1)
# tids <- sapply(strsplit(y.sirv6$ids, ";", fixed=TRUE), .subset, 2)
# length(unique(gids))
# length(unique(tids))



z <- read_tsv("day0-rep1/08_pigeon_spikein/GS_230501_Day0-1.lima.0--0.lima.IsoSeqX_bc10_5p--IsoSeqX_3p.collapsed_transcripts_classification.txt") %>% 
  filter(chrom=="SIRV6")

table(z$structural_category, z$subcategory)
table(z$associated_gene)


# read in reference
gtf <- list.files("REF", "SIRV4.gtf$", 
                  recursive = TRUE,  full.names = TRUE)
ref <- import(gtf)
ref.sirv6 <- ref[seqnames(ref)=="SIRV6"]

refs <- split(ref.sirv6, ref.sirv6$transcript_id)


cls <- read_tsv("day0-rep1/08_pigeon_spikein/GS_230501_Day0-1.lima.0--0.lima.IsoSeqX_bc10_5p--IsoSeqX_3p.collapsed_transcripts_classification.filtered_lite_classification.txt")

cls <- read_tsv("day0-rep1/08_pigeon/GS_230501_Day0-1.lima.0--0.lima.IsoSeqX_bc10_5p--IsoSeqX_3p.collapsed_transcripts_classification.filtered_lite_classification.txt")



# # reads <- read_tsv("day0-rep1/07_collapse/GS_230501_Day0-1.lima.0--0.lima.IsoSeqX_bc10_5p--IsoSeqX_3p.collapsed_transcripts.read_stat.txt")
# reads <- read_csv("day0-rep1/07_collapse/GS_230501_Day0-1.lima.0--0.lima.IsoSeqX_bc10_5p--IsoSeqX_3p.collapsed_transcripts.flnc_count.txt")
# sum(reads$fl_assoc)

# cls <- cls %>% left_join(reads, by = c("isoform" = "id"))

library(ggplot)

qs <- unique(quantile(cls$fl_assoc, p=(0:20)/20))

cls$fl_assoc_cut <- cut(cls$fl_assoc, breaks = qs, include.lowest = TRUE)

table(cls$fl_assoc_cut)

tab <- table(cls$structural_category, cls$fl_assoc_cut)
pr <- prop.table(tab, margin = 2)*100

df <- as.data.frame(pr) %>% dplyr::filter(!(Var1 %in% c("antisense", "fusion", "genic",
                                                      "intergenic", "moreJunctions"))) %>%
  transmute(expression = Var2, category = Var1, freq = Freq)

p <- ggplot(df, aes(y = freq)) +
  geom_bar(aes(x = expression, fill = category), position = "fill", stat = "identity")
p


tab <- table(cls$subcategory, cls$fl_assoc_cut)
pr <- prop.table(tab, margin = 2)*100

df <- as.data.frame(pr) %>% dplyr::filter(!(Var1 %in% c("antisense", "fusion", "genic",
                                                        "intergenic", "moreJunctions"))) %>%
  transmute(expression = Var2, category = Var1, freq = Freq)

q <- ggplot(df, aes(y = freq)) +
  geom_bar(aes(x = expression, fill = category), position = "fill", stat = "identity")
q


ggplot(cls, aes(x = cut(, 
                        breaks = unique(quantile(fl_assoc, p=(0:10/10)))), 
                y=structural_category)) + 
  geom_bar()

geom_bar(aes_string(x = "sample_id", 
                    fill = "cluster_id"), position = "fill", stat = "identity")


sirvs <- import("REF/human_plus_SIRV_Set4/SIRV_Set4_Norm_Sequences_20210507/SIRV_ERCC_longSIRV_multi-fasta_20210507.with_gene_names.gtf")
sirvs <- sirvs[seqnames(sirvs)=="SIRV6"]


x <- import("day0-rep1/08_pigeon/GS_230501_Day0-1.lima.0--0.lima.IsoSeqX_bc10_5p--IsoSeqX_3p.collapsed_transcripts.sorted.filtered_lite.gff")
x.sirv6 <- x[seqnames(x)=="SIRV6"]
x.sirv6
length(unique(x.sirv6$transcript_id))
length(unique(x.sirv6$gene_id))


# fo <- findOverlaps(sirvs, x.sirv6)
# 
# x <- import("day0-rep1/08_pigeon/GS_230501_Day0-1.lima.0--0.lima.IsoSeqX_bc10_5p--IsoSeqX_3p.collapsed_transcripts.sorted.gff")
# 
# length(unique(x$gene_id))
# length(unique(x$transcript_id))
# 
# sirvs_s <- split(ranges(sirvs), sirvs$transcript_id)
# transcripts_s <- split(ranges(x.sirv6), x.sirv6$transcript_id)
# 
# fo <- findOverlaps(sirvs_s, transcripts_s)
# f# olr <- overlapsRanges(sirvs_s, transcripts_s, fo)
# 
# 
# 
# 
# fox.sirv6 <- x[seqnames(x)=="SIRV6"]
# x.sirv6

z <- x.sirv6[x.sirv6$type=="exon"]
zs <- split(z, z$transcript_id)

mx <- max(end(ranges(x.sirv6)))

covs <- sapply(c(zs, refs), function(u) {
  as.integer(coverage(ranges(u), width = mx))
})

colnames(covs) <- c(names(zs), names(refs))
colnames(covs)[grep("PB", colnames(covs))] <- ""
covs <- t(covs)

w <- which(colSums(covs) > 0)[1]
covs <- covs[,w:ncol(covs)]


covs_noisy <- covs + rnorm(prod(dim(cov)), mean = 0, sd = .05)



pdf("this.pdf", width=20, height=10)
pheatmap(covs, cluster_cols = FALSE, 
         clustering_distance_rows = "binary")
dev.off()



sce <- SingleCellExperiment(assays=list(logcounts=t(covs_noisy)))

library(scran)

mgv <- modelGeneVar(sce)

sce <- runPCA(sce)
sce <- runUMAP(sce)

df <- data.frame(reducedDim(sce, "PCA")[,1:2], label = colnames(sce))


library(ggrepel)

ggplot(df, aes(x=PC1, y=PC2)) + geom_point() +
  geom_label_repel(data = df %>% filter(label != ""), aes(label=label))


plotUMAP(sce)


x.sirv6
ranges(x.sirv6)
end(ranges(x.sirv6))
max(end(ranges(x.sirv6)))
x.sirv6
zs
sapply(zs, function(u) as.integer(coverage(u, width=max(end(ranges(x.sirv6)))))
)
lapply(zs, function(u) as.integer(coverage(u, width=max(end(ranges(x.sirv6)))))
       mx
       lapply(zs, function(u) coverage(u, width=mx)) 
       lapply(zs, function(u) coverage(u, width=mx)[["SIRV6"]]) 
       lapply(zs, function(u) coverage(u, width=mx)[["SIRV6"]]) 
       history(100)
       
       
       
