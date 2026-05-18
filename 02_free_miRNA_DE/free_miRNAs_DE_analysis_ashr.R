library(data.table)
library(ggplot2)
library(DESeq2)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

cond_ref <- args[1] # "g1"
cond_test <- args[2] # "g2"
# cond_ref <- "g1"
# cond_test <- "g2"
LFC <- as.numeric(args[3]) # 1
# LFC <- 1

wd <- paste0('/media/jorge/ATTR_1_mac/BIOGUNE_miRNA_FIISA-01+02_smallRNA/DATA ANALYSIS 07122022/Analysis/', cond_ref, '_Vs_', cond_test)
setwd(wd)
# setwd('/media/jorge/ATTR_1_mac/BIOGUNE_miRNA_FIISA-01+02_smallRNA/DATA ANALYSIS 07122022/Analysis/g4_Vs_g2')
# setwd('/media/jorge/ATTR_1_mac/BIOGUNE_miRNA_FIISA-01+02_smallRNA/DATA ANALYSIS 07122022/Analysis/g1_Vs_g2')
# setwd('/media/jorge/ATTR_1_mac/BIOGUNE_miRNA_FIISA-01+02_smallRNA/DATA ANALYSIS 07122022/Analysis/g1_Vs_g3')
# setwd('/media/jorge/ATTR_1_mac/BIOGUNE_miRNA_FIISA-01+02_smallRNA/DATA ANALYSIS 07122022/Analysis/g1_Vs_g4')
# setwd('/media/jorge/ATTR_1_mac/BIOGUNE_miRNA_FIISA-01+02_smallRNA/DATA ANALYSIS 07122022/Analysis/g2_Vs_g3')
# setwd('/media/jorge/ATTR_1_mac/BIOGUNE_miRNA_FIISA-01+02_smallRNA/DATA ANALYSIS 07122022/Analysis/g4_Vs_g3')
# setwd('/media/jorge/ATTR_1_mac/BIOGUNE_miRNA_FIISA-01+02_smallRNA/DATA ANALYSIS 07122022/Analysis/g14_Vs_g23')
# setwd('/media/jorge/ATTR_1_mac/BIOGUNE_miRNA_FIISA-01+02_smallRNA/DATA ANALYSIS 07122022/Analysis/g124_Vs_g3')

## import data
sampleinfo <- read_tsv("samples.txt")
seqdata <- read_tsv("counts.txt", comment="#")

## reformat data
countdata <- seqdata %>% as.data.frame() %>% column_to_rownames("Geneid") %>% dplyr::select(sampleinfo$sample) %>% as.matrix()

class(countdata)  # It is in matrix format
all(countdata == floor(countdata))  # Counts must be integers (they already are). Otherwise, they have to be converted

# Visualize the first lines of the gene count matrix
head(countdata)

summary(rowSums(countdata))
head(rowSums(countdata))

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "hist_rowSums_log10.pdf")
pdf(pdf_file, width = 7, height = 5)
hist(log10(rowSums(countdata) + 1),
     breaks = 100,
     main = "Distribución log10(rowSums)",
     xlab = "log10(suma de counts + 1)")
dev.off()

##### 1. Filtering raw counts

# Before running DESeq2 (or any differential expression analysis) it is useful to pre-filter data. 
# There are computational benefits to doing this as the memory size of the objects within R will decrease 
# and DESeq2 will have less data to work through and will be faster. By removing "low quality" data, 
# it is also reduced the number of statistical tests that will be performed, which in turn reduces 
# the impact of multiple test correction and can lead to more significant genes. The amount of pre-filtering is
# up to the analyst however, it should be done in an unbiased way. DESeq2 recommends removing any gene with 
# less than 10 reads across all samples. Below, we filter a gene if at least 1 sample does not have at least 10 reads. 
# Either way, mostly what is being removed here are genes with very little evidence for expression in any sample 
# (in many cases gene with 0 counts in all samples).


# run a filtering step
# i.e. require that for every gene: at least 1 of all samples must have counts greater than 10
# get index of rows that meet this criterion and use that to subset the matrix
# note the dimensions of the matrix before and after filtering with dim()

tail(countdata)
tail(countdata >= 10)
tail(rowSums(countdata>=10))
tail(rowSums(countdata >= 10) >= 1)
tail(which(rowSums(countdata >= 10) >= 1))

dim(countdata)
filtered_countdata = countdata[which(rowSums(countdata >= 10) >= 1), ]
dim(filtered_countdata)

# Applying the other filter recommended by DEseq2

# tail(rowSums(countdata) >= 10)
# dim(countdata)
# filtered_countdata_alt = countdata[which(rowSums(countdata) >= 10), ]
# dim(filtered_countdata_alt)

##### 2. Specifying the experimental design

# As mentioned above DESeq2 also needs to know the experimental design, that is which samples belong to which condition to test. 
# The experimental design for the example dataset herein is quite simple as there are 6 samples with two conditions to compare 
# (UHR vs HBR), as such we can just create the experimental design right within R. There is one important thing to note, DESeq2 
# does not check sample names, it expects that the column names in the counts matrix we created correspond exactly to the row 
# names we specify in the experimental design.


# construct a mapping of the meta data for our experiment (comparing UHR cell lines to HBR brain tissues)
# this is defining the biological condition/label for each experimental replicate
# create a simple one column dataframe to start

metaData <- sampleinfo %>%
  column_to_rownames("sample")

metaData$status <- factor(metaData$status, levels=c(cond_ref, cond_test))
all(rownames(metaData) == colnames(filtered_countdata))  # This has to be TRUE

##### 3. Construct the DESeq2 object piecing all the data together

# With all the data properly formatted it is now possible to combine all the information required to run differental expression 
# in one object. This object will hold the input data, and intermediary calculations. It is also where the condition to test is specified.
# make deseq2 data sets
# here we are setting up our experiment by supplying: (1) the gene counts matrix, (2) the sample/replicate
# for each column, and (2) the biological conditions we wish to compare.
# this is a simple example that works for many experiments but these can also get more complex
# for example, including designs with multiple variables such as "~ group + condition",
# and designs with interactions such as "~ genotype + treatment + genotype:treatment".

dds = DESeqDataSetFromMatrix(countData = filtered_countdata, colData = metaData, design = ~status)

# the design formula above is often a point of confussion, it is useful to put in words what is happening, when we specify "design = ~Condition" we are saying
# regress gene expression on condition, or put another way model gene expression on condition
# gene expression is the response variable and condition is the explanatory variable
# you can put words to formulas using this [cheat sheet](https://www.econometrics.blog/post/the-r-formula-cheatsheet/)

##### 4. Running DESeq2

# With all the data now in place, DESeq2 can be run. Calling DESeq2 will perform the following actions:
# 
#   - Estimation of size factors. i.e. accounting for differences in sequencing depth (or library size) across samples.
#   - Estimation of dispersion. i.e. estimate the biological variability (over and above the expected variation from sampling) in gene expression across biological replicates. This is needed to assess the significance of differences across conditions. Additional work is performed to correct for higher dispersion seen for genes with low expression.
#   - Perform "independent filtering" to reduce the number of statistical tests performed (see `?results` and [this paper](https://doi.org/10.1073/pnas.0914005107) for details)
#   - Negative Binomial GLM fitting and performing the Wald statistical test
#   - Correct p values for multiple testing using the Benjamini and Hochberg method


# run the DESeq2 analysis on the "dds" object
dds = DESeq(dds)
# view the first 5 lines of the DE results
res = results(dds)
head(res, 5)

##### 5. Log-fold change shrinkage

# It is good practice to shrink the log-fold change values, this does exactly what it sounds like, reducing the 
# amount of log-fold change for genes where there are few counts which create huge variability that is not truly 
# biological signal. Consider for example a gene for two samples, one sample has 1 read, and and one sample has 6 
# reads, that is a 6 fold change, that is likely not accurate. There are a number of algorithms that can be used 
# to shrink log2 fold change, here we will use the apeglm algorithm, which does require the apeglm package to be installed.

# shrink the log2 fold change estimates
# shrinkage of effect size (log fold change estimates) is useful for visualization and ranking of genes

#   In simplistic terms, the goal of calculating "dispersion estimates" and "shrinkage" is also to account for the problem that
#   genes with low mean counts across replicates tend to have higher variability than those with higher mean counts.
#   Shrinkage attempts to correct for this. For a more detailed discussion of shrinkage refer to the DESeq2 vignette

# first get the name of the coefficient (log fold change) to shrink
resultsNames(dds)
library(ashr)

# now apply the shrinkage approach
coef_name <- paste0("status_", cond_test, "_vs_", cond_ref)
resLFC <- lfcShrink(dds, coef= coef_name, type="ashr")

# make a copy of the shrinkage results to manipulate
deGeneResult = resLFC
pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "plotMA_raw.pdf")
pdf(pdf_file, width = 7, height = 5)
plotMA(res, ylim=c(-5,5), main="Raw log2FC")
dev.off()
pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "plotMA_shrinked.pdf")
pdf(pdf_file, width = 7, height = 5)
plotMA(deGeneResult, ylim=c(-5,5), main="Shrinked log2FC")
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "raw_vs_shrinked")
pdf(pdf_file, width = 7, height = 5)
plot(
  res$log2FoldChange,
  deGeneResult$log2FoldChange,
  pch=16, cex=0.4,
  xlab="log2FC raw",
  ylab="log2FC shrinked"
)
abline(0, 1, col="red", lwd=2)
dev.off()

#### Metrics for shrinkage

raw <- res            
shr <- deGeneResult

common <- complete.cases(raw$log2FoldChange, shr$log2FoldChange)

raw_lfc <- raw$log2FoldChange[common]
shr_lfc <- shr$log2FoldChange[common]

prop_increase <- mean(abs(shr_lfc) > abs(raw_lfc))
prop_increase

df <- data.frame(
  baseMean = raw$baseMean[common],
  rawLFC   = raw_lfc,
  shrLFC   = shr_lfc
)

# Definir bins de expresión
df$bin <- cut(
  log10(df$baseMean + 1e-8), 
  breaks = seq(-2, 6, by = 0.5), 
  include.lowest = TRUE
)

summary_bins <- df %>%
  group_by(bin) %>%
  summarize(
    n = n(),
    prop_increase = mean(abs(shrLFC) > abs(rawLFC)),
    sd_raw = sd(rawLFC, na.rm = TRUE),
    sd_shr = sd(shrLFC, na.rm = TRUE),
    .groups = "drop"
  )

summary_bins

### Annotate gene symbols onto the DE results

library(org.Hs.eg.db)
library(biomaRt)

# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror='useast')

genes <- rownames(deGeneResult)   # O la lista de genes que quieras anotar

annot <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype",  "mirbase_id",
                 "mirbase_accession", "chromosome_name", "start_position", "end_position", "strand"),
  filters = "ensembl_gene_id",
  values = genes,
  mart = mart
)

deGeneResult <- as.data.frame(deGeneResult)
deGeneResult$ensembl_gene_id <- rownames(deGeneResult)

deGeneResult <- deGeneResult %>% left_join(annot, by = "ensembl_gene_id")

head(deGeneResult)

# view the top genes according to adjusted p-value
deGeneResult[order(deGeneResult$padj), ]

# order the DE results by adjusted p-value
deGeneResultSorted = deGeneResult[order(deGeneResult$padj), ]
deGeneResultSorted <- deGeneResultSorted %>%
  filter(tolower(gene_biotype) == "mirna")
dim(deGeneResultSorted)

# create a filtered data frame that limits to only the significant DE genes (adjusted p.value < 0.05)
deGeneResultSignificant = deGeneResultSorted[!is.na(deGeneResultSorted$padj) & deGeneResultSorted$padj < 0.05, ]
dim(deGeneResultSignificant)

miRNAs_sig <- deGeneResultSignificant %>%
  filter(tolower(gene_biotype) == "mirna")
miRNAs_sig  # This order matches the order of miRNAs in the BIOGUNE files

### Save in an xlsx file
# Name the directory and file where we will save the results 

file_name <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", cond_test, "_vs_", cond_ref, "_ashr_shr.xlsx")
library(openxlsx)
write.xlsx(miRNAs_sig, file = file_name, rowNames = TRUE)
file_name <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", cond_test, "_vs_", cond_ref, "_ashr_shr_all_DE_genes.xlsx")
write.xlsx(deGeneResultSorted, file = file_name, rowNames = TRUE)

### Briefly examine the top over-expressed genes

# For both conditions (HBR and UHR) lets take a look at the top n genes but this time according to fold-change instead of p-value.

# create a new copy of the data frame, sorted by log2 fold change
miRNAs_sigSortedFoldchange = miRNAs_sig %>%
  filter(abs(log2FoldChange) >= LFC) %>%
  arrange(desc(abs(log2FoldChange)))
miRNAs_sigSortedFoldchange

# Compare dimensions between the statistically significant genes before and after applying the logfoldchange filter
dim(miRNAs_sig)
dim(miRNAs_sigSortedFoldchange)

file_name <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", cond_test, "_vs_", cond_ref, "_ashr_shr_lfc_filtered.xlsx")
write.xlsx(miRNAs_sigSortedFoldchange, file = file_name, rowNames = TRUE)

# Select top 50 genes with differential expression

top50_miRNAs <- miRNAs_sigSortedFoldchange %>%
  slice_head(n = 50)
top50_miRNAs

file_name <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", cond_test, "_vs_", cond_ref, "_ashr_shr_top50.xlsx")
write.xlsx(top50_miRNAs, file = file_name, rowNames = TRUE)

# Check the intersection of which miRNAs are shared in both analyses (with shrinkage vs no shrinkage) 

directory <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/")
setwd(directory)
library(readxl)

file_name_ashr <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", cond_test, "_vs_", cond_ref, "_ashr_shr_top50.xlsx")
file_name_no_shr <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", cond_test, "_vs_", cond_ref, "_no_shr_top50.xlsx")
top50_miRNAs_ashr <- read_excel(file_name_ashr)
top50_miRNAs_no_shr <- read_excel(file_name_no_shr)
top50_miRNAs_ashr <- top50_miRNAs_ashr[, -1]
top50_miRNAs_no_shr <- top50_miRNAs_no_shr[, -1]
top50_miRNAs_ashr <- as.data.frame(top50_miRNAs_ashr)
top50_miRNAs_no_shr <- as.data.frame(top50_miRNAs_no_shr)

intersection_miRNAs <- intersect(top50_miRNAs_no_shr$ensembl_gene_id, top50_miRNAs_ashr$ensembl_gene_id)
length(intersection_miRNAs)
intersection_miRNAs

intersection_miRNAs_complete <- top50_miRNAs_ashr %>%
  filter(top50_miRNAs_ashr$ensembl_gene_id %in% intersection_miRNAs)
intersection_miRNAs_complete

shared_top_DE_genes <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", cond_test, "_vs_", cond_ref, "_shared_top_DE_genes.xlsx")
write.xlsx(intersection_miRNAs_complete, file = shared_top_DE_genes, rowNames = TRUE)

# create a convenient data structure with just the top n genes from each condition
top_bottom = bind_rows(
  head(miRNAs_sigSortedFoldchange, 10) %>% mutate(Set = "Bottom 10"),
  tail(miRNAs_sigSortedFoldchange, 10) %>% mutate(Set = "Top 10")
)
top_bottom

### Viewing pairwise sample clustering

# It may often be useful to view inter-sample relatedness. In other words, how similar or disimilar samples are to one another overall. 
# While not part of the DESeq2 package, there is a convenient library that can easily construct a hierarchically clustered heatmap from our DESeq2 data. 
# It should be noted that when doing a distance calculation using "raw count" data is not ideal, the count data should be transformed using 
# `vst()` or `rlog()` which can be performed directly on the dds object. The reason for this is described in detail in the DESeq2 manuscript, 
# suffice it to say that transforming gene variance to be more homoskedastic will make inferences of sample relatedness more interpretable.

# note that we use rlog because we don't have a large number of genes, for a typical DE experiment with 1000's of genes use the vst() function
vst <- vst(dds, blind = FALSE)
vst

# compute sample distances (the dist function uses the euclidean distance metric by default)
# in this command we will pull the rlog transformed data ("regularized" log2 transformed, see ?rlog for details) using "assay"
# then we transpose that data using t()
# then we calculate distance values using dist() 
# the distance is calculated for each vector of sample gene values, in a pairwise fashion comparing all samples

# view the first few lines of raw data
head(assay(dds))

# see the vst transformed data
head(assay(vst))

# see the impact of transposing the matrix
t(assay(vst))[1:6, 1:5]

# see the distance values
dist(t(assay(vst)))

# put it all together and store the result
sampleDists <- dist(t(assay(vst)))

# convert the distance result to a matrix
sampleDistMatrix <- as.matrix(sampleDists)

# construct clustered heatmap, important to use the computed sample distances for clustering
# install.packages("pheatmap")
library(pheatmap)

# view the distance numbers directly in the pairwise distance matrix
head(sampleDistMatrix)
pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "distance_sample_heatmap")
pdf(pdf_file, width = 7, height = 5)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists)
dev.off()
# Instead of a distance metric we could also use a similarity metric such as a Pearson correlation
# 
# There are many correlation and distance options:
# 
# Correlation: "pearson", "kendall", "spearman"
# Distance: "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"

sampleCorrs = cor(assay(vst), method = "pearson")
sampleCorrMatrix = as.matrix(sampleCorrs)
head(sampleCorrMatrix)

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "similarity_sample_heatmap")
pdf(pdf_file, width = 7, height = 5)
pheatmap(sampleCorrMatrix)
dev.off()

# Instead of boiling all the gene count data for each sample down to a distance metric you can 
# get a similar sense of the pattern by just visualizing all the genes at once
# because there are so many gene we choose not to display them
# pheatmap(mat = t(assay(vst)), show_colnames = FALSE)  # Generates fatal error

#### Supplementary R DE Visualization

# Occasionally you may wish to reformat and work with expression estimates in R in an ad hoc way. 
# Here, we provide an optional/advanced tutorial on how to visualize your results for R and perform "old school" (non-ballgown, non-DESeq2) visualization of your data.
# 
# In this tutorial you will:
#   
# * Learn basic R usage and commands (common plots, and data manipulation tasks)
# * Examine the expression estimates
# * Create an MDS plot to visualize the differences between/among replicates, library prep methods and UHR versus HBR
# * Examine the differential expression estimates
# * Visualize the expression estimates and highlight those genes that appear to be differentially expressed
# * Ask how reproducible technical replicates are.
# 
# Expression and differential expression files will be read into R. The R analysis will make use of the gene-level expression estimates from 
# differential expression results from HISAT2/htseq-count/DESeq2 (fold-changes and p-values).
# 
# Start RStudio, or launch a posit Cloud session, or if you are using AWS, navigate to the correct directory and then launch R:

#Load your libraries
library(ggplot2)
library(gplots)
library(GenomicRanges)
library(ggrepel)

#Import DE results from the HISAT2/htseq-count/DESeq2 pipeline (http://genomedata.org/cri-workshop/deseq2/DE_all_genes_DESeq2.tsv)
file_name_all_DE_genes <- file_name <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", cond_test, "_vs_", cond_ref, "_ashr_shr_all_DE_genes.xlsx")
results_genes = read_excel(file_name_all_DE_genes)

#Now, exlore the differential expression (DESeq2 results) 
head(results_genes)
dim(results_genes)

# dev.off()
#### Plot #6 - View the distribution of differential expression values as a histogram
#Display only those results that are significant according to DESeq2 (loaded above)

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "DE_distribution")
pdf(pdf_file, width = 7, height = 5)
xlab <- paste0("log2(Fold change) ", cond_test, "_vs_", cond_ref)
fc <- miRNAs_sig[,"log2FoldChange"]
lim <- max(abs(fc), na.rm = TRUE)
hist(miRNAs_sig[,"log2FoldChange"], breaks = 50, col = "seagreen", xlab = xlab, main = "Distribution of differential expression values",  xaxt = "n")
axis(
  side = 1,
  at = seq(-10, 10, by = 2)   # cambia el paso si quieres
)
abline(v = -1, col = "black", lwd = 2, lty = 2)
abline(v = 1, col = "black", lwd = 2, lty = 2)
legend("topleft", "Fold-change >= 1", lwd = 2, lty = 2)
dev.off()

#### Plot #9 - Volcano plot

# Set differential expression status for each gene - default all genes to "no change"
under_exp <- paste0("Underexpressed in ", cond_test, " vs ", cond_ref)
over_exp <- paste0("Overexpressed in ", cond_test, " vs ", cond_ref)
results_genes$diffexpressed = "No"
results_genes$diffexpressed[results_genes$log2FoldChange >= 1 & results_genes$padj <= 0.05] = over_exp
results_genes$diffexpressed[results_genes$log2FoldChange <= -1 & results_genes$padj <= 0.05] = under_exp

# write the gene names of those significantly upregulated/downregulated to a new column
results_genes$gene_label = NA
results_genes$gene_label[results_genes$diffexpressed != "No"] = results_genes$external_gene_name[results_genes$diffexpressed != "No"]

library(multiMiR)
library(miRBaseConverter)

# mir_genes <- miRNAs_sigSortedFoldchange$ensembl_gene_id
# mir_genes

mir_map <- results_genes

mir_multimir_df <- mir_map %>%
  mutate(
    # Paso 1: poner todo en minúscula temporalmente
    mir_lower = tolower(mirbase_id),
    
    # Paso 2: quitar solo loci numéricos al final (-1, -2, etc.)
    mir_step1 = gsub("-[12]$", "", mir_lower),
    
    # Paso 3: quitar a/b solo en miRNAs con 4 o más dígitos en el nombre
    mir_step2 = gsub("(mir-[0-9]{4,})[ab]$", "\\1", mir_step1),
    
    # Paso 4: volver a formato correcto Hsa-miR
    mir_canonical = gsub("^hsa-mir", "hsa-miR", mir_step2)
  ) %>%
  dplyr::select(-mir_lower, -mir_step1, -mir_step2)
mir_multimir_df

####

mir_multimir <- unique(mir_multimir_df$mir_canonical)
length(mir_multimir)

volcano_df <- mir_multimir_df %>%
  group_by(mir_canonical) %>%
  slice_min(padj, n = 1, with_ties = FALSE) %>% # elegir el precursor con menor padj
  ungroup()

# print(volcano_df[volcano_df$external_gene_name == "MIR4449",], n = 6, width=Inf)

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "DEG_volcano_plot")
pdf(pdf_file, width = 7, height = 5)
ggplot(data = volcano_df[volcano_df$diffexpressed != "No",], aes(x = log2FoldChange, y = -log10(padj), label = mir_canonical, color = diffexpressed)) +
  xlab("log2Foldchange") +
  scale_color_manual(name = "Differentially expressed", values=c("red", "blue")) +
  geom_point() +
  theme_minimal() +
  geom_text_repel(force = 12, max.overlaps = 10) +
  geom_vline(xintercept = c(-1, 1), col = "red") +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  geom_point(data = volcano_df[volcano_df$diffexpressed == "No",], aes(x = log2FoldChange, y = -log10(padj)), colour = "black")
dev.off()

sum(is.na(results_genes$log2FoldChange))
sum(is.na(results_genes$padj))
sum(results_genes$padj == 0, na.rm = TRUE)

#### Importing DE results for gage

# Before we perform the pathway analysis we need to read in our differential expression results from the previous DE analysis.
# Now lets go ahead and load [GAGE](https://bioconductor.org/packages/release/bioc/html/gage.html) and some other useful packages. 

library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)
library(gage)

#### Setting up gene set databases

# In order to perform our pathway analysis we need a list of pathways and their respective genes. 
# There are many databases that contain collections of genes (or gene sets) that can be used to understand whether a set of mutated 
# or differentially expressed genes are functionally related.  Some of these resources include: [GO](http://www.geneontology.org/), 
# [KEGG](https://www.kegg.jp), [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb), and [WikiPathways](https://www.wikipathways.org/index.php/WikiPathways). 
# For this exercise we are going to investigate [GO](http://www.geneontology.org/) and [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb).  
# The [GAGE](https://bioconductor.org/packages/release/bioc/html/gage.html) package has a function for querying [GO](http://www.geneontology.org/) 
# in real time: [go.gsets()](https://www.rdocumentation.org/packages/gage/versions/2.22.0/topics/go.gsets). This function takes a species as an argument
# and will return a list of gene sets and some helpful meta information for subsetting these lists. If you are unfamiliar with [GO](http://www.geneontology.org/), 
# it is helpful to know that GO terms are categorized into three gene ontologies: "Biological Process", "Molecular Function", and "Cellular Component". 
# This information will come in handy later in our exercise. GAGE does not provide a similar tool to investigate the gene sets available in MSigDB.
# Fortunately, MSigDB provides a  downloadable `.gmt` file for all gene sets. This format is easily read into GAGE using a function called 
# [readList()](https://www.rdocumentation.org/packages/gage/versions/2.22.0/topics/readList). If you check out [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb) 
# you will see that there are 8 unique gene set collections, each with slightly different features. For this exercise we will use the
# [C8 - cell type signature gene sets collection](https://www.gsea-msigdb.org/gsea/msigdb/collection_details.jsp#C8), which is a collection of gene sets 
# that contain cluster markers for cell types identified from single-cell sequencing studies of human tissue.

# Get target genes

precursor_miRNAs <- mir_map %>%
  filter(abs(log2FoldChange) >= LFC & padj <= 0.05) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  dplyr::select(mirbase_id)
precursor_miRNAs

# DE_miRNAs_canonical <- volcano_df %>%
#   filter(abs(log2FoldChange) >= LFC & padj <= 0.05) %>%
#   arrange(desc(abs(log2FoldChange))) %>%
#   dplyr::select(mir_canonical)
# DE_miRNAs_canonical

miRNANames <- c(precursor_miRNAs$mirbase_id)
mature <- miRNA_PrecursorToMature(miRNANames)
mature_mirnas <- unique(
  na.omit(
    c(mature$Mature1, mature$Mature2)
  )
)
mature_mirnas
length(mature_mirnas)

targets <- get_multimir(
  mirna = mature_mirnas, # DE_miRNAs_canonical,
  summary = TRUE,
  table = "validated"   #option: all --> validated + predicted
)

### check the miRNAs that do not have validated targets

# no_targets <- get_multimir(
#   mirna = "hsa-miR-3135a",
#   summary = TRUE,
#   table = "validated"   # validated + predicted
# )

### check the miRNAs that do not have validated targets

targets_df <- targets@data
head(targets_df)

genes_target <- unique(targets_df$target_symbol)
length(genes_target)
genes_target[1:10]
unique(targets_df$mature_mirna_id)

### Define background genes

background_genes <- rownames(dds)
library(AnnotationDbi)
library(org.Hs.eg.db)

bg_entrez <- mapIds(
  org.Hs.eg.db,
  keys = background_genes,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

bg_entrez <- bg_entrez[!is.na(bg_entrez)]
bg_entrez

gene_entrez <- unique(targets_df$target_entrez)
gene_entrez <- gene_entrez[!is.na(gene_entrez)]
gene_entrez

library(clusterProfiler)

ego_bp <- enrichGO(
  gene          = gene_entrez,
  universe      = bg_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_mf <- enrichGO(
  gene = gene_entrez,
  universe = bg_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_cc <- enrichGO(
  gene = gene_entrez,
  universe = bg_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

head(ego_bp@result)

# ego_bp <- ego_bp[ego_bp@result$p.adjust < 0.05]
# ego_cc <- ego_cc[ego_cc@result$p.adjust < 0.05]
# ego_mf <- ego_mf[ego_mf@result$p.adjust < 0.05]

ego_bp_simplified <- clusterProfiler::simplify(ego_bp) # Filter redundant GO terms
ego_mf_simplified <- clusterProfiler::simplify(ego_mf) # Filter redundant GO terms
ego_cc_simplified <- clusterProfiler::simplify(ego_cc) # Filter redundant GO terms

bp_all <- ego_bp@result$Description
bp_simpl <- ego_bp_simplified@result$Description
length(bp_all)
length(bp_simpl)

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "BP_miRNA_targets")
pdf(pdf_file, width = 7, height = 9)
dotplot(ego_bp_simplified, showCategory = 15) +
  ggtitle("GO Biological Process – miRNA targets")
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "MF_miRNA_targets")
pdf(pdf_file, width = 7, height = 9)
dotplot(ego_mf_simplified, showCategory = 15) +
  ggtitle("GO Molecular Function – miRNA targets")
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "CC_miRNA_targets")
pdf(pdf_file, width = 7, height = 9)
dotplot(ego_cc_simplified, showCategory = 15) +
  ggtitle("GO Cellular Component – miRNA targets")
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "BP_miRNA_targets_barplot")
pdf(pdf_file, width = 7, height = 9)
barplot(ego_bp_simplified, showCategory = 15)
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "MF_miRNA_targets_barplot")
pdf(pdf_file, width = 7, height = 9)
barplot(ego_mf_simplified, showCategory = 15)
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "CC_miRNA_targets_barplot")
pdf(pdf_file, width = 7, height = 9)
barplot(ego_cc_simplified, showCategory = 15)
dev.off()

# Separated GO enrichment
miRNA_up <- mir_map[mir_map$log2FoldChange >= 1,] 
miRNA_down <- mir_map[mir_map$log2FoldChange <= -1,] 

# Upregulated

miRNA_up <- mir_map %>%
  filter(log2FoldChange >= LFC & padj <= 0.05) %>%
  arrange(desc(log2FoldChange))

# mir_multimir_up <- miRNA_up$mir_canonical
# mir_multimir_up
# length(mir_multimir_up)

miRNANames_up <- c(miRNA_up$mirbase_id)
mature_up <- miRNA_PrecursorToMature(miRNANames_up)
mature_mirnas_up <- unique(
  na.omit(
    c(mature_up$Mature1, mature_up$Mature2)
  )
)
mature_mirnas_up
length(mature_mirnas_up)

targets_up <- get_multimir(
  mirna = mature_mirnas_up,
  summary = TRUE,
  table = "validated"   # validated + predicted
)

targets_df_up <- targets_up@data
head(targets_df_up)

genes_target_up <- unique(targets_df_up$target_symbol)
length(genes_target_up)
genes_target_up[1:10]
unique(targets_df_up$mature_mirna_id)

### Define background genes

gene_entrez_up <- unique(targets_df_up$target_entrez)
gene_entrez_up <- gene_entrez_up[!is.na(gene_entrez_up)]
gene_entrez_up

ego_bp_up <- enrichGO(
  gene          = gene_entrez_up,
  universe      = bg_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_mf_up <- enrichGO(
  gene = gene_entrez_up,
  universe = bg_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_cc_up <- enrichGO(
  gene = gene_entrez_up,
  universe = bg_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

head(ego_bp@result)

ego_bp_up_simplified <- clusterProfiler::simplify(ego_bp_up) # Filter redundant GO terms
ego_mf_up_simplified <- clusterProfiler::simplify(ego_mf_up) # Filter redundant GO terms
ego_cc_up_simplified <- clusterProfiler::simplify(ego_cc_up) # Filter redundant GO terms

bp_all <- ego_bp_up@result$Description
bp_simpl <- ego_bp_up_simplified@result$Description
length(bp_all)
length(bp_simpl)

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "BP_miRNA_targets_up")
pdf(pdf_file, width = 7, height = 9)
dotplot(ego_bp_up_simplified, showCategory = 15) +
  ggtitle("GO Biological Process – miRNA targets")
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "MF_miRNA_targets_up")
pdf(pdf_file, width = 7, height = 9)
dotplot(ego_mf_up_simplified, showCategory = 15) +
  ggtitle("GO Molecular Function – miRNA targets")
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "CC_miRNA_targets_up")
pdf(pdf_file, width = 7, height = 9)
dotplot(ego_cc_up_simplified, showCategory = 15) +
  ggtitle("GO Cellular Component – miRNA targets")
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "BP_miRNA_targets_barplot_up")
pdf(pdf_file, width = 7, height = 9)
barplot(ego_bp_up_simplified, showCategory = 15)
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "MF_miRNA_targets_barplot_up")
pdf(pdf_file, width = 7, height = 9)
barplot(ego_mf_up_simplified, showCategory = 15)
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "CC_miRNA_targets_barplot_up")
pdf(pdf_file, width = 7, height = 9)
barplot(ego_cc_up_simplified, showCategory = 15)
dev.off()

# Downregulated  --> in g1 vs g2 there is no significant GO terms because padjust > 0.05 for the only miRNA DE

miRNA_down <- mir_map %>%
  filter(log2FoldChange <= -LFC & padj <= 0.05) %>%
  arrange(desc(log2FoldChange))

# mir_multimir_down <- miRNA_down$mir_canonical
# mir_multimir_down
# length(mir_multimir_down)

miRNANames_down <- c(miRNA_down$mirbase_id)
mature_down <- miRNA_PrecursorToMature(miRNANames_down)
mature_mirnas_down <- unique(
  na.omit(
    c(mature_down$Mature1, mature_down$Mature2)
  )
)
mature_mirnas_down
length(mature_mirnas_down)

targets_down <- get_multimir(
  mirna = mature_mirnas_down,
  summary = TRUE,
  table = "validated"   # validated + predicted
)

targets_df_down <- targets_down@data
head(targets_df_down)

genes_target_down <- unique(targets_df_down$target_symbol)
length(genes_target_down)
genes_target_down[1:10]
unique(targets_df_down$mature_mirna_id)

gene_entrez_down <- unique(targets_df_down$target_entrez)
gene_entrez_down <- gene_entrez_down[!is.na(gene_entrez_down)]
gene_entrez_down

ego_bp_down <- enrichGO(
  gene          = gene_entrez_down,
  universe      = bg_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_mf_down <- enrichGO(
  gene = gene_entrez_down,
  universe = bg_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_cc_down <- enrichGO(
  gene = gene_entrez_down,
  universe = bg_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

head(ego_bp@result)

ego_bp_down_simplified <- clusterProfiler::simplify(ego_bp_down) # Filter redundant GO terms
ego_mf_down_simplified <- clusterProfiler::simplify(ego_mf_down) # Filter redundant GO terms
ego_cc_down_simplified <- clusterProfiler::simplify(ego_cc_down) # Filter redundant GO terms

bp_all <- ego_bp@result$Description
bp_simpl <- ego_bp_simplified@result$Description
length(bp_all)
length(bp_simpl)

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "BP_miRNA_targets_down")
pdf(pdf_file, width = 7, height = 9)
dotplot(ego_bp_down, showCategory = 15) +
  ggtitle("GO Biological Process – miRNA targets")
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "MF_miRNA_targets_down")
pdf(pdf_file, width = 7, height = 9)
dotplot(ego_mf_down, showCategory = 15) +
  ggtitle("GO Molecular Function – miRNA targets")
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "CC_miRNA_targets_down")
pdf(pdf_file, width = 7, height = 9)
dotplot(ego_cc_down, showCategory = 15) +
  ggtitle("GO Cellular Component – miRNA targets")
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "BP_miRNA_targets_barplot_down")
pdf(pdf_file, width = 7, height = 9)
barplot(ego_bp_down, showCategory = 15)
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "MF_miRNA_targets_barplot_down")
pdf(pdf_file, width = 7, height = 9)
barplot(ego_mf_down, showCategory = 15)
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "CC_miRNA_targets_barplot_down")
pdf(pdf_file, width = 7, height = 9)
barplot(ego_cc_down, showCategory = 15)
dev.off()

### KEGG 

# Upregulated

ekegg_up <- enrichKEGG(
  gene = gene_entrez_up,
  organism = "hsa",
  pvalueCutoff = 0.05
)

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "kegg_up")
pdf(pdf_file, width = 7, height = 5)
dotplot(ekegg_up)
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "kegg_up_barplot")
pdf(pdf_file, width = 7, height = 5)
barplot(ekegg_up)
dev.off()

genes_kegg_up <- strsplit(ekegg_up@result$geneID, "/")  # Counts the number of enriched pathways that contain that gene (e.g., 10000 --> 102, gene 10000 appears in 102 enriched pathways by DE miRNAs)
hub_kegg_genes_up <- table(unlist(genes_kegg_up)) %>%
  sort(decreasing = TRUE)

# Downregulated

ekegg_down <- enrichKEGG(
  gene = gene_entrez_down,
  organism = "hsa",
  pvalueCutoff = 0.05
)

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "kegg_down")
pdf(pdf_file, width = 7, height = 5)
dotplot(ekegg_down)
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", "kegg_down_barplot")
pdf(pdf_file, width = 7, height = 5)
barplot(ekegg_down)
dev.off()

genes_kegg_down <- strsplit(ekegg_down@result$geneID, "/")  # Counts the number of enriched pathways that contain that gene (e.g., 10000 --> 102, gene 10000 appears in 102 enriched pathways by DE miRNAs)
hub_kegg_genes_down <- table(unlist(genes_kegg_down)) %>%
  sort(decreasing = TRUE)

# More plots 

### Upregulated

cnetplot(
  ego_bp_up,
  showCategory = 5
)

ego_bp_up_top <- ego_bp_up_simplified
ego_bp_up_top@result <- ego_bp_up_top@result[
  order(ego_bp_up_top@result$p.adjust),
][1:20, ]

ego_bp_up_top <- enrichplot::pairwise_termsim(ego_bp_up_top)
emapplot(
  ego_bp_up_top,
  showCategory = 20,
  layout = "kk",
  cluster.params = list(label_format = 0))

### Downregulated  --> In this case the number of genes is very low, that's why some plots do not work.

cnetplot(
  ego_bp_down,
  showCategory = 5
)

ego_bp_down_top <- ego_bp_down_simplified
ego_bp_down_top@result <- ego_bp_down_top@result[
  order(ego_bp_down_top@result$p.adjust),
][1:20, ]

ego_bp_down_top <- enrichplot::pairwise_termsim(ego_bp_down_top)
emapplot(
  ego_bp_down_top,
  showCategory = 20,
  layout = "kk",
  cluster.params = list(label_format = 0))

# Identify HUB genes

# Upregulated

hub_genes_up <- targets_df_up %>% 
  count(target_entrez)%>%
  arrange(desc(n))%>%    # n gives the number of evidences of a gene being regulated by an miRNA
  filter(n >= 3)
hub_genes_up

biologically_relevant_genes_up <- intersect(hub_genes_up$target_entrez, names(hub_kegg_genes_up)) # Genes that are regulated by several miRNAs and also implied in many pathways
length(biologically_relevant_genes_up)
biologically_relevant_genes_up_name <- targets_df_up[targets_df_up$target_entrez %in% biologically_relevant_genes_up, "target_symbol"]
length(biologically_relevant_genes_up_name)
biologically_relevant_genes_up_name <- unique(biologically_relevant_genes_up_name)
# biologically_relevant_genes_up_name
length(biologically_relevant_genes_up_name)

# Downregulated

hub_genes_down <- targets_df_down %>%
  count(target_entrez)%>%
  arrange(desc(n))%>%
  filter(n >= 3)
hub_genes_down

biologically_relevant_genes_down <- intersect(hub_genes_down$target_entrez, names(hub_kegg_genes_down)) # Genes that are regulated by several miRNAs and also implied in many pathways
length(biologically_relevant_genes_down)
biologically_relevant_genes_down_name <- targets_df_down[targets_df_down$target_entrez %in% biologically_relevant_genes_down, "target_symbol"]
length(biologically_relevant_genes_down_name)
biologically_relevant_genes_down_name <- unique(biologically_relevant_genes_down_name)
# biologically_relevant_genes_up_name
length(biologically_relevant_genes_down_name)

### Integration with gene expression
# Upregulated

targets_df_up_DE <- targets_df_up %>%
  left_join(deGeneResult[, c("ensembl_gene_id", "log2FoldChange")], 
            by = c("target_ensembl" = "ensembl_gene_id"))
head(targets_df_up_DE)

targets_df_up_DE <-unique(targets_df_up_DE[, c("target_entrez","target_symbol","log2FoldChange")])

targets_df_up_DE <- targets_df_up_DE %>%
  mutate(direction = ifelse(log2FoldChange < 0, "expected", "unexpected"))

table(targets_df_up_DE$direction)
wilcox.test(targets_df_up_DE$log2FoldChange, mu = 0)

wilcox.test(
  targets_df_up_DE$log2FoldChange,
  deGeneResult$log2FoldChange
)


