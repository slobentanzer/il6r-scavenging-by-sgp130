# Re-analysis of differential gene expression from a cohort of SARS patients
# as described in "Scavenging of Interleukin 6 Receptor by Bioidentical Recombinant gp130 as 
# Intervention in Covid-19 Exacerbation" by Sebastian Lobentanzer (2020) 
# Preprint available at https://osf.io/3gwmp

setwd("~/GitHub/il6r-scavenging-by-sgp130")
rm(list = ls())

pacman::p_load(
  tidyverse,
  GEOquery,
  ggplot2,
  limma,
  reshape2,
  quantro,
  doParallel)

options(stringsAsFactors = F)


# load data and prepare metadata ####
# Download deposited data of Cameron et al. (2007). Interferon-Mediated Immunopathological Events Are
# Associated with Atypical Innate and Adaptive Immune Responses in Patients with Severe Acute
# Respiratory Syndrome. J. Virol. 81(16): 8692–8706. doi: 10.1128/JVI.00527-07
geo0 <- getGEO("GSE5972")
geo <- geo0[[1]]

pData(geo) %>% head()

# meta curation
meta0 <- pData(geo)
head(meta0)
meta <- data.frame(row.names = rownames(meta0), title = meta0$title, ch2 = meta0$`Tissue:ch2`)

meta$id <- unlist(lapply(strsplit(meta$ch2, ", "), "[", 2))
meta$id <- as.factor(gsub("Patient id:", "", meta$id))

meta$sex <- unlist(lapply(strsplit(meta$ch2, ", "), "[", 4))
meta$sex <- as.factor(gsub("Sex:", "", meta$sex))

meta$status <- unlist(lapply(strsplit(meta$ch2, ", "), "[", 6))
meta$status <- as.factor(gsub("Status:", "", meta$status))
levels(meta$status) <- c("conv", "ctl", "post", "pre")
meta$status <- relevel(meta$status, "ctl")

meta$age <- unlist(lapply(strsplit(meta$ch2, ", "), "[", 3))
meta$age <- as.numeric(gsub("Age:", "", meta$age))

meta <- meta[, c("id", "sex", "age", "status", "title")]

pData(geo) <- meta
features <- fData(geo)

geo

# quantro ####
# A test for when to use global normalization methods, such as quantile normalization.
ex <- exprs(geo)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
exsq <- 2^ex
qxsq <- as.numeric(quantile(exsq, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))

head(exsq)
qtest <- quantro(object = ex, groupFactor = meta$status)
qtest
registerDoParallel(cores = 7)
qtestPerm <- quantro(ex, groupFactor = meta$status, B = 1000)
qtestPerm
quantroPlot(qtestPerm) #no further normalisation required
ggsave("img/quantro_perm_plot.pdf", width = 8, height = 5)

# limma ####
# Differential expression analysis of time point subgroups of patients vs controls.
# >simple model ####
design <- model.matrix(~ 0 + status, data = meta)
colnames(design) <- gsub("status", "", colnames(design), fixed = T)

fit <- lmFit(geo, design)
contrast.matrix <- makeContrasts(conv-ctl, pre-ctl, post-ctl, post-pre, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
plotMA(fit2, coef = 1)
plotMA(fit2, coef = 2)
plotMA(fit2, coef = 3)
plotMA(fit2, coef = 4)

conv.vs.ctl <- topTable(fit2, coef = 1, adjust.method = "BH", number = Inf)
pre.vs.ctl <- topTable(fit2, coef = 2, adjust.method = "BH", number = Inf)
post.vs.ctl <- topTable(fit2, coef = 3, adjust.method = "BH", number = Inf)
post.vs.pre <- topTable(fit2, coef = 4, adjust.method = "BH", number = Inf)

# >gene selection for plotting ####
# pre-nadir samples
prevc.ils <- pre.vs.ctl[grepl("interleukin", pre.vs.ctl$GENE_NAME) & pre.vs.ctl$adj.P.Val<.05 | grepl("STAT", pre.vs.ctl$GENE_SYMBOL) & pre.vs.ctl$adj.P.Val<.05,
                        c("GENE_SYMBOL", "logFC", "adj.P.Val", "GENE_NAME")]

# post-nadir samples
postvc.ils <- post.vs.ctl[grepl("interleukin", post.vs.ctl$GENE_NAME) & post.vs.ctl$adj.P.Val<.05 | grepl("STAT", post.vs.ctl$GENE_SYMBOL) & post.vs.ctl$adj.P.Val<.05,
                          c("GENE_SYMBOL", "logFC", "adj.P.Val", "GENE_NAME")]

# labels
prevc.ils$group <- "pre"
postvc.ils$group <- "post"
# combine
ils <- rbind(prevc.ils, postvc.ils)

# annotate IL6ST (missing annotation)
features$GENE_SYMBOL[grep("BQ023177", features$CLONE_ACC)] <- "IL6ST"
# collect relevant genes
idx <- grep("interleukin", features$GENE_NAME)
idx <- c(idx, grep("STAT[1-6]", features$GENE_SYMBOL))
idx <- c(idx, grep("BQ023177", features$CLONE_ACC))
features[idx, ]
ids <- features$ID[features$ID %in% idx & features$ID %in% rownames(ils)]

exil6 <- ex[ids, ]
rownames(exil6) <- paste(features$GENE_SYMBOL[features$ID %in% idx & features$ID %in% rownames(ils)], rownames(exil6), sep = "_")

ex_melt <- melt(exil6)
ex_melt$status <- meta$status[match(ex_melt$Var2, rownames(meta))]
ex_melt$status <- factor(ex_melt$status, levels = c("ctl", "pre", "post", "conv"))
# preview plot
ggplot(ex_melt, aes(status, value, fill = status)) + geom_violin() +
  facet_wrap("Var1")

# remove uninformative or erroneous
rem <- c(807, 1814, 2999, 9696)
exil6 <- data.frame(ex[ids, ])
# rename
exil6$name <- paste(features$GENE_SYMBOL[features$ID %in% rownames(exil6) & features$ID %in% rownames(ils)], " (", rownames(exil6), ")", sep = "")
ex_melt <- melt(exil6[!rownames(exil6) %in% rem, ], id.vars = c("name"))
ex_melt$status <- meta$status[match(ex_melt$variable, rownames(meta))]
ex_melt$status <- factor(ex_melt$status, levels = c("ctl", "pre", "post", "conv"))
ex_melt$name <- factor(ex_melt$name, levels = sort(unique(ex_melt$name)))
ggplot(ex_melt, aes(status, value, fill = status)) + geom_violin() +
  facet_wrap("name") + theme_bw()

# >full model for publication####
design <- model.matrix(~ 0 + status + age + sex, data = meta)
colnames(design) <- gsub("status", "", colnames(design), fixed = T)

fit <- lmFit(geo, design)
contrast.matrix <- makeContrasts(conv-ctl, pre-ctl, post-ctl, post-pre, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
plotMA(fit2, coef = 1)
plotMA(fit2, coef = 2)
plotMA(fit2, coef = 3)
plotMA(fit2, coef = 4)

conv.vs.ctl <- topTable(fit2, coef = 1, adjust.method = "BH", number = Inf)
pre.vs.ctl <- topTable(fit2, coef = 2, adjust.method = "BH", number = Inf)
post.vs.ctl <- topTable(fit2, coef = 3, adjust.method = "BH", number = Inf)
post.vs.pre <- topTable(fit2, coef = 4, adjust.method = "BH", number = Inf)

# write out
write.csv(conv.vs.ctl, file = "out/limma_de_conv_vs_ctl.csv", row.names = F, quote = T)
write.csv(pre.vs.ctl, file = "out/limma_de_pre_vs_ctl.csv", row.names = F, quote = T)
write.csv(post.vs.ctl, file = "out/limma_de_post_vs_ctl.csv", row.names = F, quote = T)
write.csv(post.vs.pre, file = "out/limma_de_post_vs_pre.csv", row.names = F, quote = T)

exil6 <- data.frame(ex[ids, ])
# rename
exil6$name <- paste(features$GENE_SYMBOL[features$ID %in% rownames(exil6) & features$ID %in% rownames(ils)], " (", rownames(exil6), ")", sep = "")
ex_melt <- melt(exil6[!rownames(exil6) %in% rem, ], id.vars = c("name"))
ex_melt$status <- meta$status[match(ex_melt$variable, rownames(meta))]
ex_melt$status <- factor(ex_melt$status, levels = c("ctl", "pre", "post", "conv"))
ex_melt$name <- factor(ex_melt$name, levels = sort(unique(ex_melt$name)))
ggplot(ex_melt, aes(status, value, fill = status)) + geom_violin() +
  facet_wrap("name") + theme_bw() #essentially same as simple model
ggsave("img/DE_ILs_violin.pdf", width = 12, height = 10)
