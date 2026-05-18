library(openxlsx)
library(biomaRt)
library(data.table)
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(tibble)
library(multiMiR)
library(miRBaseConverter)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DESeq2)

cond_test <- "g4"
cond_ref <- c("g2", "g3")
LFC <- 1

# Read differentially expressed genes from exosomes (ATTRv vs CC and ATTRv vs HC) and circulating (ATTRv vs CC and ATTRv vs HC) miRNAs
files_exo <- paste0("/home/jorge/Escritorio/RNA_seq_luca/exosomes/results/LFC", LFC, "/", cond_test, cond_ref, "/", cond_test, "_vs_", cond_ref, "_ashr_shr_lfc_filtered.xlsx")
files_circ <- paste0("/home/jorge/Escritorio/RNA_seq_luca/results/LFC", LFC, "/", cond_test, cond_ref, "/", cond_test, "_vs_", cond_ref, "_ashr_shr_lfc_filtered.xlsx")
exo_V_CC <- read.xlsx(files_exo[1])
exo_V_HC <- read.xlsx(files_exo[2])
circ_V_CC <- read.xlsx(files_circ[1])
circ_V_HC <- read.xlsx(files_circ[2])

### Free
# Find the common DE miRNAs
common_DE_genes_free <- circ_V_CC %>%
  inner_join(circ_V_HC, by = "mirbase_id") %>%
  pull(mirbase_id) %>%
  unique()
file_name <- paste0("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/free/DE_miRNAs_free_comparisons_precursors.xlsx")
write.xlsx(as.data.frame(common_DE_genes_free), file = file_name)

# Mature miRNAs with 5p and 3p arms 
mature_free <- miRNA_PrecursorToMature(common_DE_genes_free)
mature_free <- mature_free %>%
  filter(!is.na(Mature1) | !is.na(Mature2))

mature_mirnas_free <- c(mature_free$Mature1, mature_free$Mature2)
mature_mirnas_free <- unique(mature_mirnas_free[!is.na(mature_mirnas_free)])
file_name <- paste0("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/free/DE_miRNAs_free_comparisons.xlsx")
write.xlsx(as.data.frame(mature_mirnas_free), file = file_name)

# Find DE miRNA targets 
targets_free <- get_multimir(
  mirna = mature_mirnas_free, # DE_miRNAs_canonical,
  summary = TRUE,
  table = "validated"   #option: all --> validated + predicted
)

file_name <- paste0("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/free/targets_df_free.xlsx")
write.xlsx(as.data.frame(targets_free@data), file = file_name)

### Exosomes

# Find the common DE miRNAs
common_DE_genes_exo <- exo_V_CC %>%
  inner_join(exo_V_HC, by = "mirbase_id") %>%
  pull(mirbase_id) %>%
  unique()
file_name <- paste0("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/exosomes/DE_miRNAs_exo_comparisons_precursors.xlsx")
write.xlsx(as.data.frame(common_DE_genes_exo), file = file_name)

# Mature miRNAs with 5p and 3p arms 
mature_exo <- miRNA_PrecursorToMature(common_DE_genes_exo)
mature_exo <- mature_exo %>%
  filter(!is.na(Mature1) | !is.na(Mature2))

mature_mirnas_exo <- c(mature_exo$Mature1, mature_exo$Mature2)
mature_mirnas_exo <- unique(mature_mirnas_exo[!is.na(mature_mirnas_exo)])
file_name <- paste0("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/exosomes/DE_miRNAs_exo_comparisons.xlsx")
write.xlsx(as.data.frame(mature_mirnas_exo), file = file_name)

# Find DE miRNA targets 
targets_exo <- get_multimir(
  mirna = mature_mirnas_exo, # DE_miRNAs_canonical,
  summary = TRUE,
  table = "validated"   #option: all --> validated + predicted
)

file_name <- paste0("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/exosomes/targets_df_exo.xlsx")
write.xlsx(as.data.frame(targets_exo@data), file = file_name)

targets_read_free <- read.xlsx('/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/free/targets_df_free.xlsx')
targets_read_exo <- read.xlsx('/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/exosomes/targets_df_exo.xlsx')
targets_read_free
targets_read_exo

read_mature_miRNAs_free <- read.xlsx('/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/free/DE_miRNAs_free_comparisons.xlsx')
read_mature_miRNAs_exo <- read.xlsx('/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/exosomes/DE_miRNAs_exo_comparisons.xlsx')

read_precursor_miRNAs_free <- read.xlsx('/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/free/DE_miRNAs_free_comparisons_precursors.xlsx')
read_precursor_miRNAs_exo <- read.xlsx('/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/exosomes/DE_miRNAs_exo_comparisons_precursors.xlsx')

targets_combined <- rbind(targets_read_free, targets_read_exo)

# Rename columns to bind dataframes
read_precursor_miRNAs_free <- read_precursor_miRNAs_free %>%
  rename(common_DE_genes = common_DE_genes_free)
read_precursor_miRNAs_exo <- read_precursor_miRNAs_exo %>%
  rename(common_DE_genes = common_DE_genes_exo)
read_mature_miRNAs_free <- read_mature_miRNAs_free %>%
  rename(common_DE_genes = mature_mirnas_free)
read_mature_miRNAs_exo <- read_mature_miRNAs_exo %>%
  rename(common_DE_genes = mature_mirnas_exo)

mature_miRNAs_combined <- rbind(read_mature_miRNAs_free, read_mature_miRNAs_exo)
precursor_miRNAs_combined <- bind_rows(read_precursor_miRNAs_free, read_precursor_miRNAs_exo)

targets_read <- targets_combined[, c(
  "mature_mirna_id",
  "target_symbol",
  "target_ensembl"
)]
targets_read <- unique(targets_read)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

mir_annot <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "gene_biotype",
    "mirbase_id",
    "mirbase_accession"
  ),
  filters = "mirbase_id",
  values = tolower(precursor_miRNAs_combined$common_DE_genes),
  mart = mart
)

mir_map <- mir_annot[, c("ensembl_gene_id", "mirbase_id")]
mir_map <- mir_map %>%
  distinct(mirbase_id, .keep_all = TRUE)

mirnas_ensembl <- mir_map$ensembl_gene_id
targets_ensembl <- unique(targets_read$target_ensembl)
# mat_mirnas <- mat[intersect(rownames(mat), mirnas_ensembl), ]
# mat_mirnas
# mat_targets <- mat[intersect(rownames(mat), targets_ensembl), ]
# mat_targets

targets_read <- targets_read %>%
  mutate(precursor = miRNA_MatureToPrecursor(mature_mirna_id)$Precursor)
targets_read

precursors <- unique(targets_read$precursor)
length(precursors)

mir_ensembl <- getBM(
  attributes = c("mirbase_id", "ensembl_gene_id",
                 "external_gene_name",
                 "gene_biotype"),
  filters = "mirbase_id",
  values = precursors,
  mart = mart
)

mir_ensembl
mir_ensembl <- mir_ensembl %>%
  distinct(mirbase_id, .keep_all = TRUE)

targets_read <- targets_read %>%
  left_join(mir_ensembl, by = c("precursor" = "mirbase_id"))
targets_read  

####

samples_exo_V_CC$fraction  <- "exo"
samples_exo_V_HC$fraction  <- "exo"
samples_circ_V_CC$fraction <- "free"
samples_circ_V_HC$fraction <- "free"

# Collect all the samples in a single data.frame
coldata_all <- do.call(rbind, list(
  samples_exo_V_CC,
  samples_exo_V_HC,
  samples_circ_V_CC,
  samples_circ_V_HC
))

# Remove duplicates by sample name
coldata_all <- coldata_all[!duplicated(coldata_all$sample), ]

# Assure it's a dataframe not a tibble
coldata_all <- as.data.frame(coldata_all)

# Fix rownames as sample names
rownames(coldata_all) <- coldata_all$sample

# Status as a factor
coldata_all$status <- factor(coldata_all$status, levels = c(cond_ref, cond_test))

# Fraction as a factor
coldata_all$fraction <- factor(coldata_all$fraction, levels = c("free", "exo"))

# Verify columns coincide with counts
all(rownames(coldata_all) == colnames(filtered_countdata))

dds_all <- DESeqDataSetFromMatrix(
  countData = filtered_countdata,
  colData = coldata_all,
  design = ~ fraction + status
)

dds_all <- estimateSizeFactors(dds_all)

# Transformation of variance-stabilizing
vsd <- vst(dds_all, blind = FALSE)

# Transformed expression matrix 
mat <- assay(vsd)

# Free samples
samples_free <- rownames(coldata_all)[coldata_all$fraction == "free"]
mat_free <- mat[, samples_free]

# Exosomes samples
samples_exo <- rownames(coldata_all)[coldata_all$fraction == "exo"]
mat_exo <- mat[, samples_exo]

dim(mat_free)
dim(mat_exo)

# Subset in mRNA targets' and miRNAs' matrices
mir_exo   <- mat_exo[intersect(rownames(mat_exo), mir_ensembl$ensembl_gene_id), ]  # salen 95 filas porque de 97 mir_ensembl$ensembl_gene_id hay dos que están duplicados
mir_free  <- mat_free[intersect(rownames(mat_free), mir_ensembl$ensembl_gene_id), ]

mrna_exo  <- mat_exo[intersect(rownames(mat_exo), targets_ensembl), ]
mrna_free <- mat_free[intersect(rownames(mat_free), targets_ensembl), ]

dim(mir_exo)
dim(mir_free)
dim(mrna_exo)
dim(mrna_free)

calc_cor <- function(mir_mat, mrna_mat, pairs) {
  
  out <- lapply(seq_len(nrow(pairs)), function(i) {
    
    mir_id  <- pairs$ensembl_gene_id[i]
    gene_id <- pairs$target_ensembl[i]
    
    if (!mir_id %in% rownames(mir_mat) ||
        !gene_id %in% rownames(mrna_mat)) {
      return(c(NA, NA))
    }
    
    mir_vals  <- as.numeric(mir_mat[mir_id, ])
    gene_vals <- as.numeric(mrna_mat[gene_id, ])
    
    test <- cor.test(mir_vals, gene_vals, method="spearman")
    
    c(rho = test$estimate, p = test$p.value)
  })
  
  as.data.frame(do.call(rbind, out))
}

res_exo  <- calc_cor(mir_exo, mrna_exo, targets_read)
res_free <- calc_cor(mir_free, mrna_free, targets_read)

targets_read$rho_exo  <- res_exo$rho
targets_read$p_exo    <- res_exo$p

targets_read$rho_free <- res_free$rho
targets_read$p_free   <- res_free$p

targets_read$FDR_exo  <- p.adjust(targets_read$p_exo, method="BH")
targets_read$FDR_free <- p.adjust(targets_read$p_free, method="BH")

validated_exo <- subset(targets_read,
                        rho_exo < -0.3 & FDR_exo < 0.05)

validated_free <- subset(targets_read,
                         rho_free < -0.3 & FDR_free < 0.05)

targets_read$type <- "none"

targets_read$type[
  targets_read$rho_exo < -0.3 & targets_read$FDR_exo < 0.05 &
    !(targets_read$rho_free < -0.3 & targets_read$FDR_free < 0.05)
] <- "exo_specific"

targets_read$type[
  targets_read$rho_free < -0.3 & targets_read$FDR_free < 0.05 &
    !(targets_read$rho_exo < -0.3 & targets_read$FDR_exo < 0.05)
] <- "free_specific"

targets_read$type[
  targets_read$rho_exo < -0.3 & targets_read$FDR_exo < 0.05 &
    targets_read$rho_free < -0.3 & targets_read$FDR_free < 0.05
] <- "shared"

interest <- targets_read %>%
  filter(type!="NA")

interest <- targets_read %>%
  filter(type=="shared")
interest_exo <- targets_read %>%
  filter(type=="exo_specific")
interest_free <- targets_read %>%
  filter(type=="free_specific")

file_name <- paste0("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/shared/anticorrelated_miRNA_mRNA_target_pairs_shared.xlsx")
write.xlsx(as.data.frame(interest), file = file_name)
file_name <- paste0("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/exosomes/anticorrelated_miRNA_mRNA_target_pairs_exosomes.xlsx")
write.xlsx(as.data.frame(interest_exo), file = file_name)
file_name <- paste0("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/free/anticorrelated_miRNA_mRNA_target_pairs_free.xlsx")
write.xlsx(as.data.frame(interest_free), file = file_name)

interest_entrez <- mapIds(
  org.Hs.eg.db,
  keys = interest_exo$target_ensembl, # chnage between interest_free and interest_exo
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

background_genes <- rownames(dds_all)

bg_entrez <- mapIds(
  org.Hs.eg.db,
  keys = background_genes,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

bg_entrez <- bg_entrez[!is.na(bg_entrez)]
bg_entrez

# Exosomes
ekegg <- enrichKEGG(
  gene = interest_entrez,
  universe = bg_entrez,
  organism = "hsa",
  pvalueCutoff = 0.05
)

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/exosomes/kegg")
pdf(pdf_file, width = 7, height = 9)
dotplot(ekegg)
dev.off()

ego_bp <- enrichGO(
  gene          = interest_entrez,
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
  gene          = interest_entrez,
  universe      = bg_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_cc <- enrichGO(
  gene          = interest_entrez,
  universe      = bg_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/exosomes/BP")
pdf(pdf_file, width = 7, height = 9)
dotplot(ego_bp, showCategory = 15) +
  ggtitle("GO Biological Process")
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/exosomes/MF")
pdf(pdf_file, width = 7, height = 9)
dotplot(ego_mf, showCategory = 15) +
  ggtitle("MF Biological Process")
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/exosomes/CC")
pdf(pdf_file, width = 7, height = 9)
dotplot(ego_cc, showCategory = 15) +
  ggtitle("CC Biological Process")
dev.off()

# Free

interest_entrez <- mapIds(
  org.Hs.eg.db,
  keys = interest_free$target_ensembl, # chnage between interest_free and interest_exo
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

ekegg <- enrichKEGG(
  gene = interest_entrez,
  universe = bg_entrez,
  organism = "hsa",
  pvalueCutoff = 0.05
)

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/free/kegg")
pdf(pdf_file, width = 7, height = 9)
dotplot(ekegg)
dev.off()

ego_bp <- enrichGO(
  gene          = interest_entrez,
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
  gene          = interest_entrez,
  universe      = bg_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_cc <- enrichGO(
  gene          = interest_entrez,
  universe      = bg_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/free/BP")
pdf(pdf_file, width = 7, height = 9)
dotplot(ego_bp, showCategory = 15) +
  ggtitle("GO Biological Process")
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/free/MF")
pdf(pdf_file, width = 7, height = 9)
dotplot(ego_mf, showCategory = 15) +
  ggtitle("MF Biological Process")
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/free/CC")
pdf(pdf_file, width = 7, height = 9)
dotplot(ego_cc, showCategory = 15) +
  ggtitle("CC Biological Process")
dev.off()

# Shared

interest_entrez <- mapIds(
  org.Hs.eg.db,
  keys = interest$target_ensembl, # chnage between interest_free and interest_exo
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

ekegg <- enrichKEGG(
  gene = interest_entrez,
  universe = bg_entrez,
  organism = "hsa",
  pvalueCutoff = 0.05
)

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/shared/kegg")
pdf(pdf_file, width = 7, height = 9)
dotplot(ekegg)
dev.off()

ego_bp <- enrichGO(
  gene          = interest_entrez,
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
  gene          = interest_entrez,
  universe      = bg_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_cc <- enrichGO(
  gene          = interest_entrez,
  universe      = bg_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/shared/BP")
pdf(pdf_file, width = 7, height = 9)
dotplot(ego_bp, showCategory = 15) +
  ggtitle("GO Biological Process")
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/shared/MF")
pdf(pdf_file, width = 7, height = 9)
dotplot(ego_mf, showCategory = 15) +
  ggtitle("MF Biological Process")
dev.off()

pdf_file <- paste0("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/shared/CC")
pdf(pdf_file, width = 7, height = 9)
dotplot(ego_cc, showCategory = 15) +
  ggtitle("CC Biological Process")
dev.off()

igsf <- strsplit(ekegg$geneID[1], "/")[[1]]
igsf_symbols <- bitr(
  igsf,
  fromType = "ENTREZID",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
)

adherens <- strsplit(ekegg$geneID[2], "/")[[1]]
adherens_symbols <- bitr(
  adherens,
  fromType = "ENTREZID",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
)

# All comparisons

DE_mRNAs_all_comparisons <- read.xlsx("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_analysis/all_comparisons_DE_mRNA_targets.xlsx")
common <- interest %>%
  filter(target_symbol %in% DE_mRNAs_all_comparisons$target_symbol)

# Free

DE_mRNAs_free_g4g2 <- read.xlsx("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_analysis/free/g4g2_protein_coding_DE.xlsx")
DE_mRNAs_free_g4g3 <- read.xlsx("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_analysis/free/g4g3_protein_coding_DE.xlsx")
DE_mRNAs_free <- Reduce(intersect, list(DE_mRNAs_free_g4g2$external_gene_name, DE_mRNAs_free_g4g3$external_gene_name))
common_free <- interest_free %>%
  filter(target_symbol %in% DE_mRNAs_free)

# Exosomes

DE_mRNAs_exo_g4g2 <- read.xlsx("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_analysis/exosomes/g4g2_protein_coding_DE.xlsx")
DE_mRNAs_exo_g4g3 <- read.xlsx("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_analysis/exosomes/g4g3_protein_coding_DE.xlsx")
DE_mRNAs_exo <- Reduce(intersect, list(DE_mRNAs_exo_g4g2$external_gene_name, DE_mRNAs_exo_g4g3$external_gene_name))
common_exo <- interest_exo %>%
  filter(target_symbol %in% DE_mRNAs_exo)

file_name <- paste0("/home/jorge/Escritorio/RNA_seq_luca/protein_coding_new_analysis/exosomes/anticorrelated_miRNA_and_mRNA_DE_validated_pairs.xlsx")
write.xlsx(as.data.frame(common_exo), file = file_name)



