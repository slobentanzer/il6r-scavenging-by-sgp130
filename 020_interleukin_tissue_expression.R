# Analysis of tissue specific gene expression from regulatory circuits
# as described in "Scavenging of Interleukin 6 Receptor by Bioidentical Recombinant gp130 as 
# Intervention in Covid-19 Exacerbation" by Sebastian Lobentanzer (2020) 
# Preprint available at https://osf.io/3gwmp

setwd("~/GitHub/il6r-scavenging-by-sgp130")
rm(list = ls())

pacman::p_load(ggplot2, reshape2)

options(stringsAsFactors = F)

# load
tissues <- readRDS("data/analysed_tissues.rds")
ensg2symbol <- readRDS(file = "data/ensg2symbol_all.rds")

# TRANSCRIPTIONAL ACTIVITY IN IMMUNE TISSUES####
# Tissue transcription factor activities were downloaded from the accompanying website of
# Marbach et al. (2016) Tissue-specific regulatory circuits reveal variable modular perturbations 
# across complex diseases. Nature Methods, 13(4):366–370. doi:10.1038/nmeth.3799.

tissues

# >create expression table of each tissue####
files <- paste0("data/", list.files("data/"))
# gene selection
idx <- which(ensg2symbol$hgnc_symbol %in% c("IL6", "IL6R", "IL6ST", "IL1R1", "IL1RL1", "IL17RE", "IL22RA1", "ACE2"))
exp <- ensg2symbol[idx, ]
for(tis in tissues){
  idx <- which(tissues == tis)
  message(idx)
  file <- files[grep(tis, files)]
  tar <- readRDS(file)
  tar_sum <- aggregate(tar$r.tfa, by = list(g.name = tar$g.name, g.ensg = tar$g.ensg), "sum")
  tar_sum <- tidyr::separate_rows(tar_sum, g.ensg, sep = ", ")
  
  #select genes
  col <- tar_sum$x[match(exp$ensembl_gene_id, tar_sum$g.ensg)]
  col[is.na(col)] <- 0
  exp <- cbind(exp, col)
  colnames(exp)[idx+2] <- tis
}
head(exp)
colnames(exp)[1:2] <- c("ensg", "name")

il6exp <- exp[exp$name %in% c("IL6R", "IL6ST", "IL17RE", "IL22RA1", "IL1R1", "IL1RL1", "ACE2"),]
il6exp_sum <- sort(colSums(il6exp["73172", 3:ncol(il6exp)]))

exp_melt <- melt(il6exp, id.vars = c("ensg", "name"))
exp_melt$variable <- as.character(exp_melt$variable)
exp_melt$variable <- factor(exp_melt$variable, levels = names(il6exp_sum))
exp_melt$name <- factor(exp_melt$name, levels = c("IL6R", "IL6ST", "IL1R1", "IL1RL1", "IL17RE", "IL22RA1", "ACE2"))
ggplot(exp_melt, aes(variable, value, color = name, group = name)) + geom_line(aes(linetype = name)) +
  coord_flip() +
  theme_bw()
ggsave("img/marbach_immune_il_receptor_expression.pdf", width = 10, height = 7)
