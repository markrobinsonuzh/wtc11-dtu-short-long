library(tidyr)
library(dplyr)
library(readr)
library(tibble)
library(DRIMSeq)
library(pheatmap)

trans_counts_file <- snakemake@input[["counts"]]

trans_counts <- read.table(trans_counts_file,row.names=1) %>% 
  replace(is.na(.), 0) %>%
  as.data.frame %>%
  transmute(gene_id = associated_gene,
         feature_id = isoform,
         day0_1 = `0-1`,
         day0_2 = `0-2`,
         day0_3 = `0-3`,
         day5_1 = `5-1`,
         day5_2 = `5-2`,
         day5_3 = `5-3`)

md <- data.frame(sample_id = paste0("day",c(0,0,0,5,5,5),"_",c(1:3,1:3)),
                 group = rep(c("day0","day5"), each=3))
md

d <- dmDSdata(counts = trans_counts, samples = md)
d

head(counts(d))

plotData(d)


d <- dmFilter(d, min_samps_gene_expr = 3, min_samps_feature_expr = 3,
              min_gene_expr = 100, min_feature_expr = 20)

design_full <- model.matrix(~ group, data = samples(d))
design_full


d <- dmPrecision(d, design = design_full)

d <- dmFit(d, design = design_full, verbose = 1)

d <- dmTest(d, coef = "groupday5", verbose = 1)

res <- results(d)
res %>% arrange(adj_pvalue) %>% head(20)

sum(res$adj_pvalue < .05, na.rm = TRUE)


plotProportions(d, gene_id = "RBM7", group_variable = "group")
trans_counts %>% filter(gene_id=="RBM7")

plotProportions(d, gene_id = "RPF2", group_variable = "group")
trans_counts %>% filter(gene_id=="RPF2")

plotProportions(d, gene_id = "RBM7", group_variable = "group")
trans_counts %>% filter(gene_id=="RBM7")

plotProportions(d, gene_id = "ARHGDIB", group_variable = "group")
trans_counts %>% filter(gene_id=="ARHGDIB")


out <- trans_counts %>% left_join(res)

# write statistics to table
write.table(out, "08062023_DTU_gloria_pacbio.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)