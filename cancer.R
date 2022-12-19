#Import maEndToEnd and do not display the package start up message
suppressPackageStartupMessages({library("maEndToEnd")})
#General Bioconductor packages
library(devtools)
library(remotes)
library(Biobase)
library(oligoClasses)
#Annotation and data import packages
library(ArrayExpress)
library(pd.hugene.1.0.st.v1)
library(hugene10sttranscriptcluster.db)
library(pd.hg.u133.plus.2)
library(hgu133plus2.db)
library(oligo)
library(arrayQualityMetrics)
#Analysis and statistics packages
library(limma)
library(topGO)
library(ReactomePA)
library(clusterProfiler)
#Plotting and color options packages
library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(tidyr)
#Helpers
library(stringr)
library(matrixStats)
library(genefilter)
library(openxlsx)
library("evaluate")
library("hexbin")
library("ggnewscale")
library(gtExtras)
library(gt)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~Import Gene Data~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#read sample and relationship data format
sdrf_location <- file.path('/Users/mahdi/repos/bme460_finalproject/cancer/study2data/E-MEXP-993', "E-MEXP-993.sdrf.txt")
SDRF <- read.delim(sdrf_location)


rownames(SDRF) <- SDRF$Array.Data.File
SDRF <- AnnotatedDataFrame(SDRF)

SDRF

raw_data <- oligo::read.celfiles(filenames = file.path('/Users/mahdi/repos/bme460_finalproject/cancer/study2data/E-MEXP-993/',
                                                       SDRF$Array.Data.File),
                                 verbose = FALSE, phenoData = SDRF)
stopifnot(validObject(raw_data))

raw_data
ncol(Biobase::pData(raw_data))

head(Biobase::pData(raw_data)) %>%
  gt() %>%
  gt_theme_538()

head(Biobase::pData(raw_data))

Biobase::pData(raw_data) <- Biobase::pData(raw_data)[,c('Source.Name',
                                                        'Characteristics..DiseaseState.',
                                                        'Characteristics..Individual.',
                                                        'Factor.Value..DifferentiationState.')]

head(Biobase::pData(raw_data)) %>%
  gt() %>%
  gt_theme_538() %>%
  tab_header(title = "Sample of data")


#spelling fixes 
Biobase::pData(raw_data) <- Biobase::pData(raw_data) %>% 
  mutate(Factor.Value..DifferentiationState. = case_when(
    Factor.Value..DifferentiationState. %in% c("Prostate Epithalial Stem Cells") ~ "Prostate Epithelial Stem Cells",
    Factor.Value..DifferentiationState. %in% c("Prostate Epithalial transit amplifying cells", "Prostate Epithelial Transit amplifying cells", "Prostate Epithelial transit amplifying cells") ~ "Prostate Epithelial Transit Amplifying Cells"
    ,TRUE ~ Factor.Value..DifferentiationState.
  )
  )


head(Biobase::pData(raw_data))


head(Biobase::pData(raw_data)) %>%
  gt() %>%
  gt_theme_538() %>%
  tab_header(title = "Sample of data")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~Data Analysis~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(Biobase::pData(raw_data))

#checking for outliers

Biobase::exprs(raw_data)[1:4, 1:4]



#perform PCA on log2 intesisty scale of expressions
exp_raw <- log2(Biobase::exprs(raw_data))
exp_raw[1:4, 1:4]
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     Disease = pData(raw_data)$Characteristics..DiseaseState.,
                     DiseaseState = pData(raw_data)$Factor.Value..DifferentiationState.,
                     Individual = pData(raw_data)$Characteristics..Individual.)
ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Disease, colour = DiseaseState)) +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) +
  scale_color_manual(values = c("darkorange2", "dodgerblue4"))


#box plot of intensities, one box per individual microarray
oligo::boxplot(raw_data, target = "core",
               main = "Boxplot of log2-intensitites for the raw data")

#perform quality check
# arrayQualityMetrics(expressionset = raw_data,
#                     outdir = '/Users/mahdi/repos/bme460_finalproject/outputData993/',
#                     force = TRUE, do.logtransform = TRUE,
#                     intgroup = c("Characteristics..DiseaseState.", "Factor.Value..DifferentiationState."))

#relative log expression check
palmieri_eset <- oligo::rma(raw_data, normalize = FALSE)

#reshape data
row_medians_assayData <- 
  Biobase::rowMedians(as.matrix(Biobase::exprs(palmieri_eset)))

RLE_data <- sweep(Biobase::exprs(palmieri_eset), 1, row_medians_assayData)

RLE_data <- as.data.frame(RLE_data)

RLE_data_gathered <-
  tidyr::gather(RLE_data, Genes, log2_expression_deviation)

ggplot2::ggplot(RLE_data_gathered, aes(Genes,
                                       log2_expression_deviation)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylim(c(-2, 2)) + 
  theme(axis.text.x = element_text(colour = "aquamarine4",
                                   angle = 60, size = 6.5, hjust = 1 ,
                                   face = "bold"))

#normalize background
palmieri_eset_norm <- oligo::rma(raw_data)


print(palmieri_eset_norm)
exp_palmieri <- Biobase::exprs(palmieri_eset_norm)
PCA <- prcomp(t(exp_palmieri), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

# Disease = pData(raw_data)$Characteristics..DiseaseState.,
# Dstate = pData(raw_data)$Factor.Value..DifferentiationState.,
# Individual = pData(raw_data)$Characteristics..Individual.

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     DISEASESTATE =
                       Biobase::pData(palmieri_eset_norm)$Characteristics..DiseaseState.,
                     DIFFERENTIATIONSTATE =
                       Biobase::pData(palmieri_eset_norm)$Factor.Value..DifferentiationState.)


ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = DISEASESTATE, colour = DIFFERENTIATIONSTATE)) +
  ggtitle("PCA plot of the calibrated, summarized data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) +
  scale_color_manual(values = c("darkorange2", "dodgerblue4"))


#----------

disease_state <- ifelse(str_detect(pData
                                    (palmieri_eset_norm)$Characteristics..DiseaseState.,
                                   "Prostate Adenocarcinoma"),"Prostate Adenocarcinoma","Benign Prostatic Hyperplasia")

print(pData(palmieri_eset_norm)$Characteristics..DiseaseState.)
print('~~~~~~~~~~~~~~~')
print(pData(palmieri_eset_norm)$Factor.Value..DifferentiationState.)


diff_state <- ifelse(str_detect(pData
                                   (palmieri_eset_norm)$Factor.Value..DifferentiationState.,
                                "Prostate Epithelial Stem Cells"), "Prostate Epithelial Stem Cells","Prostate Epithelial Transit Amplifying Cells")


annotation_for_heatmap <-
  data.frame(DiffState = diff_state, Disease = disease_state)

row.names(annotation_for_heatmap) <- row.names(pData(palmieri_eset_norm))


dists <- as.matrix(dist(t(exp_palmieri), method = "manhattan"), rownames.force=NA)

rownames(dists) <- row.names(pData(palmieri_eset_norm))
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))
colnames(dists) <- NULL
diag(dists) <- NA

print(annotation_for_heatmap)


ann_colors <- list(
  DiffState = c("Prostate Epithelial Stem Cells" = "aquamarine", "Prostate Epithelial Transit Amplifying Cells" ="darkgreen"),
  Disease = c( "Prostate Adenocarcinoma" = "blue4",  "Benign Prostatic Hyperplasia"	 = "darkorange2")
)

pheatmap(dists, col = (hmcol),
         annotation_row = annotation_for_heatmap,
         annotation_colors = ann_colors,
         legend = TRUE,
         treeheight_row = 0,
         legend_breaks = c(min(dists, na.rm = TRUE),
                           max(dists, na.rm = TRUE)),
         legend_labels = (c("small distance", "large distance")),
         main = "Clustering heatmap for the calibrated samples")

#soft intensity based filtering
palmieri_medians <- rowMedians(Biobase::exprs(palmieri_eset_norm))

hist_res <- hist(palmieri_medians, 100, col = "cornsilk1", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

#threshhold filtering
man_threshold <- 4

hist_res <- hist(palmieri_medians, 100, col = "cornsilk", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

abline(v = man_threshold, col = "coral4", lwd = 2)

#!!!!!!!!!!
no_of_samples <-
  table(paste0(pData(palmieri_eset_norm)$Characteristics..DiseaseState., "_",
               pData(palmieri_eset_norm)$Factor.Value..DifferentiationState.))
no_of_samples

#check how many are filtered out
samples_cutoff <- min(no_of_samples)
idx_man_threshold <- apply(Biobase::exprs(palmieri_eset_norm), 1,
                           function(x){
                             sum(x > man_threshold) >= samples_cutoff})

idx_man_threshold

#~~~~~~~~~~~~~
palmieri_manfiltered <- subset(palmieri_eset_norm, idx_man_threshold)

anno_palmieri <- AnnotationDbi::select(hgu133plus2.db,
                                       keys = (featureNames(palmieri_manfiltered)),
                                       columns = c("SYMBOL", "GENENAME"),
                                       keytype = "PROBEID")
# columns(hgu133plus2.db)
# keytypes(hgu133plus2.db)
# head(anno_palmieri)
#"16650045" %in% keys(hugene10sttranscriptcluster.db, "GENEID") 

anno_palmieri <- subset(anno_palmieri, !is.na(SYMBOL))

#remove multiple mapping
#palmieri_eset_norm
anno_grouped <- group_by(anno_palmieri, PROBEID)
anno_summarized <-
  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))
head(anno_summarized)

anno_filtered <- filter(anno_summarized, no_of_matches > 1)
head(anno_filtered)

probe_stats <- anno_filtered
nrow(probe_stats)

ids_to_exlude <- (featureNames(palmieri_manfiltered) %in% probe_stats$PROBEID)
table(ids_to_exlude)


palmieri_final <- subset(palmieri_manfiltered, !ids_to_exlude)
validObject(palmieri_final)

head(anno_palmieri)

fData(palmieri_final)$PROBEID <- rownames(fData(palmieri_final))

fData(palmieri_final) <- left_join(fData(palmieri_final), anno_palmieri)

# restore rownames after left_join
rownames(fData(palmieri_final)) <- fData(palmieri_final)$PROBEID
validObject(palmieri_final)

#fitting linear model to data
individual <-
  as.character(Biobase::pData(palmieri_final)$Characteristics..Individual.)

Dstate <- str_replace_all(Biobase::pData(palmieri_final)$Characteristics..DiseaseState.,
                          " ", "_")


Dstate <- ifelse(Dstate == "Prostate_Adenocarcinoma",
                 "PA", "BH")



DiffState <- 
  str_replace_all(Biobase::pData(palmieri_final)$Factor.Value..DifferentiationState.,
                  " ", "_")


DiffState <- ifelse(DiffState == "Prostate_Epithelial_Stem_Cells",
                 "STEM", "TRAN")

# DiffState <-
#   ifelse(str_detect(Biobase::pData(palmieri_final)$Factor.Value..DifferentiationState.,
#                     "Prostate_Epithelial_Stem_Cells"), "STEM", "TRAN")

#~~~~~


# i_STEM <- individual[DiffState == "STEM"]
# design_palmieri_STEM <- model.matrix(~ 0 + Dstate[DiffState == "STEM"] + i_STEM)
# colnames(design_palmieri_STEM)[1:2] <- c("PA", "BH")
# rownames(design_palmieri_STEM) <- i_STEM
# 
# 
# i_TRAN <- individual[DiffState == "TRAN"]
# design_palmieri_TRAN <- model.matrix(~ 0 + Dstate[DiffState == "TRAN"] + i_TRAN )
# colnames(design_palmieri_TRAN)[1:2] <- c("PA", "BH")
# rownames(design_palmieri_TRAN) <- i_TRAN

i_PA <- individual[Dstate == "PA"]
design_palmieri_PA <- model.matrix(~ 0 + DiffState[Dstate == "PA"] + i_PA)
colnames(design_palmieri_PA)[1:2] <- c("STEM", "TRAN")
rownames(design_palmieri_PA) <- i_PA


i_BH <- individual[Dstate == "BH"]
design_palmieri_BH <- model.matrix(~ 0 + DiffState[Dstate == "BH"] + i_BH )
colnames(design_palmieri_BH)[1:2] <- c("STEM", "TRAN")
rownames(design_palmieri_BH) <- i_BH


disease_PA <- DiffState[Dstate == "PA"]
crat_expr <- Biobase::exprs(palmieri_final)["1553016_at", Dstate == "PA"]
crat_data <- as.data.frame(crat_expr)
colnames(crat_data)[1] <- "org_value"
crat_data <- mutate(crat_data, individual = i_PA, disease_PA)

crat_data$disease_PA <- factor(crat_data$disease_PA, levels = c("STEM", "TRAN"))

ggplot(data = crat_data, aes(x = disease_PA, y = org_value,
                             group = individual, color = individual)) +
  geom_line() +
  ggtitle("Expression changes for ADGRF3 gene")


crat_coef <- lmFit(palmieri_final[,Dstate == "PA"],
                   design = design_palmieri_PA)$coefficients["1553016_at",]

crat_coef


crat_fitted <- design_palmieri_PA %*% crat_coef
rownames(crat_fitted) <- names(crat_expr)
colnames(crat_fitted) <- "fitted_value"

crat_fitted

crat_data$fitted_value <- crat_fitted
ggplot(data = crat_data, aes(x = disease_PA, y = fitted_value,
                             group = individual, color = individual)) +
  geom_line() +
  ggtitle("Fitted expression changes for the ADGRF3 gene")

# Differential expression analysis of the CRAT gene. 
# In order to test whether the gene is differentially expressed or not, 
# a t-test with the null hypothesis that there is no difference 
# in the expression between non-inflamed and inflamed tissue is carried out. 
# Our blocking design is conceptually similar to a paired t-test for which the 
# statistic is given by: t = d/(s/n^(1/2))

crat_PA <- na.exclude(crat_data$org_value[Dstate == "PA"])
crat_BH <- na.exclude(crat_data$org_value[Dstate == "BH"])
res_t <- t.test(crat_PA ,crat_BH , paired = FALSE)
res_t


#The result of the eBayes() step is that 
#the individual variances are shrunken towards the prior value.


contrast_matrix_PA <- makeContrasts(STEM-TRAN, levels = design_palmieri_PA)
palmieri_fit_PA <- eBayes(contrasts.fit(lmFit(palmieri_final[,Dstate == "PA"],
                                              design = design_palmieri_PA),
                                        contrast_matrix_PA))

contrast_matrix_BH <- makeContrasts(STEM-TRAN, levels = design_palmieri_BH)
palmieri_fit_BH <- eBayes(contrasts.fit(lmFit(palmieri_final[,Dstate == "BH"],
                                              design = design_palmieri_BH),
                                        contrast_matrix_BH))




table_PA <- topTable(palmieri_fit_PA, number = Inf)
head(table_PA)

hist(table_PA$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "Prostate Epithelial Stem Cells vs Prostate Epithelial Transit Amplifying Cells - Prostate Adenocarcinoma", xlab = "p-values")

table_BH <- topTable(palmieri_fit_BH, number = Inf)
head(table_BH)

hist(table_BH$P.Value, col = brewer.pal(3, name = "Set2")[2],
     main = "Prostate Epithelial Stem Cells vs Prostate Epithelial Transit Amplifying Cells - Benign Prostatic Hyperplasia", xlab = "p-values")
#a p-value of 0.001 was used as a significance cutoff. Using this 
#we get 947 genes identified as differentially expressed for UC:

nrow(subset(table_PA, P.Value < 0.001))


tail ( subset (table_PA, P.Value < 0.001 ))


fpath <- system.file("extdata", "palmieri_DE_res.xlsx", package = "maEndToEnd")
palmieri_DE_res <- sapply(1:4, function(i) read.xlsx(cols = 1, fpath,
                                                     sheet = i, startRow = 4))

names(palmieri_DE_res) <- c("PA_UP", "PA_DOWN", "BH_UP", "BH_DOWN")
palmieri_DE_res <- lapply(palmieri_DE_res, as.character)
paper_DE_genes_PA <- Reduce("c", palmieri_DE_res[1:2])
paper_DE_genes_BH <- Reduce("c", palmieri_DE_res[3:4])

overlap_PA <- length(intersect(subset(table_PA, P.Value < 0.001)$SYMBOL,
                               paper_DE_genes_PA)) / length(paper_DE_genes_PA)


overlap_BH <- length(intersect(subset(table_BH, P.Value < 0.001)$SYMBOL,
                               paper_DE_genes_BH)) / length(paper_DE_genes_BH)
overlap_PA

overlap_BH

total_genenumber_PA <- length(subset(table_PA, P.Value < 0.001)$SYMBOL)
total_genenumber_BH <- length(subset(table_BH, P.Value < 0.001)$SYMBOL)

total_genenumber_PA

total_genenumber_BH

#Visualization of DE analysis results - volcano plot

#PA
volcano_names <- ifelse(abs(palmieri_fit_PA$coefficients)>=1,
                        palmieri_fit_PA$genes$SYMBOL, NA)

volcanoplot(palmieri_fit_PA, coef = 1L, style = "p-value", highlight = 100,
            names = volcano_names,
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)

#BH
volcano_names_BH <- ifelse(abs(palmieri_fit_BH$coefficients)>=1,
                        palmieri_fit_BH$genes$SYMBOL, NA)

volcanoplot(palmieri_fit_BH, coef = 1L, style = "p-value", highlight = 100,
            names = volcano_names_BH,
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)


DE_genes_PA <- subset(table_PA, adj.P.Val < 0.1)$PROBEID

#Matching the background set of genes
back_genes_idx <- genefilter::genefinder(palmieri_final,
                                         as.character(DE_genes_PA),
                                         method = "manhattan", scale = "none")



back_genes_idx <- sapply(back_genes_idx, function(x)x$indices)

back_genes <- featureNames(palmieri_final)[back_genes_idx]
back_genes <- setdiff(back_genes, DE_genes_PA)


intersect(back_genes, DE_genes_PA)

length(back_genes)

multidensity(list(
  all = table_PA[,"AveExpr"] ,
  fore = table_PA[DE_genes_PA , "AveExpr"],
  back = table_PA[rownames(table_PA) %in% back_genes, "AveExpr"]),
  col = c("#e46981", "#ae7ee2", "#a7ad4a"),
  xlab = "mean expression",
  main = "DE genes for Prostate Adenocarcinoma-background-matching")



#running topGO
gene_IDs <- rownames(table_PA)
in_universe <- gene_IDs %in% c(DE_genes_PA, back_genes)
in_selection <- gene_IDs %in% DE_genes_PA
all_genes <- in_selection[in_universe]
all_genes <- factor(as.integer(in_selection[in_universe]))
names(all_genes) <- gene_IDs[in_universe]

top_GO_data <- new("topGOdata", ontology = "BP", allGenes = all_genes,
                   nodeSize = 10, annot = annFUN.db, affyLib = "hgu133plus2.db")


result_top_GO_elim <-
  runTest(top_GO_data, algorithm = "elim", statistic = "Fisher")
result_top_GO_classic <-
  runTest(top_GO_data, algorithm = "classic", statistic = "Fisher")


res_top_GO <- GenTable(top_GO_data, Fisher.elim = result_top_GO_elim,
                       Fisher.classic = result_top_GO_classic,
                       orderBy = "Fisher.elim" , topNodes = 100)
genes_top_GO <- printGenes(top_GO_data, whichTerms = res_top_GO$GO.ID,
                           chip = "hgu133plus2.db", geneCutOff = 1000)
res_top_GO$sig_genes <- sapply(genes_top_GO, function(x){
  str_c(paste0(x[x$'raw p-value' == 2, "Symbol.id"],";"),
        collapse = "")
})
head(res_top_GO[,1:8], 20)


#Visualization of the GO-analysis results
showSigOfNodes(top_GO_data, score(result_top_GO_elim), firstSigNodes = 3,
               useInfo = 'def')


#A pathway enrichment analysis using reactome
entrez_ids <- mapIds(hgu133plus2.db,
                     keys = rownames(table_PA),
                     keytype = "PROBEID",
                     column = "ENTREZID")

reactome_enrich <- enrichPathway(gene = entrez_ids[DE_genes_PA],
                                 universe = entrez_ids[c(DE_genes_PA,
                                                         back_genes)],
                                 organism = "human",
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 0.9,
                                 readable = TRUE)

reactome_enrich@result$Description <- paste0(str_sub(
  reactome_enrich@result$Description, 1, 20),
  "...")

head(summary(reactome_enrich))[1:6]

#Visualizing the reactome based analysis results

barplot(reactome_enrich)

x2 <- pairwise_termsim(reactome_enrich)

emapplot(x2)

