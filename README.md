# miRNA Differential Expression & miRNA–mRNA Correlation Analysis

Differential expression analysis of **exosomal and circulating free miRNAs** in ATTRv patients versus healthy controls (HC) and cardiac controls, or patients with other cardiac affections (CC), followed by target identification and miRNA–mRNA correlation and functional enrichment analysis.

---

## Project structure

```
miRNA-mRNA-DE-analysis/
│
├── 01_exosomal_miRNA_DE/
│   └── exo_DE_analysis_ashr.R
│
├── 02_free_miRNA_DE/
│   └── free_miRNAs_DE_analysis_ashr.R
│
└── 03_correlation_analysis/
    └── correlation_miRNA_mRNA_new_analysis.R
```

---

## Pipeline overview

Steps 01 and 02 are **independent and can be run in parallel**. Step 03 depends on the outputs of both.

```
01_exosomal_miRNA_DE    ──┐
                           ├──► 03_correlation_analysis
02_free_miRNA_DE        ──┘
```

---

## Script descriptions

### `01_exosomal_miRNA_DE/exo_DE_analysis_ashr.R`
Differential expression analysis of **exosomal miRNAs**, with functional enrichment. Run via command-line arguments:

```bash
Rscript exo_DE_analysis_ashr.R <cond_ref> <cond_test> <LFC>
# e.g.:
Rscript exo_DE_analysis_ashr.R g2 g4 1
```

- Loads raw counts and sample metadata
- Filters low-count genes
- Runs DESeq2 with `~status` design
- Applies `ashr` log2 fold-change shrinkage
- Annotates results with Ensembl gene symbols and miRBase IDs via `biomaRt`
- Exports significant DE miRNAs to `.xlsx` (full, LFC-filtered, top 50, and shared with no-shrinkage results)
- Converts precursor to mature miRNA IDs via `miRBaseConverter`
- Retrieves validated miRNA–mRNA targets with `multiMiR`
- GO enrichment (BP, MF, CC) and KEGG via `clusterProfiler`, separately for up- and down-regulated miRNAs
- Identifies hub genes (targets regulated by multiple miRNAs and implicated in multiple pathways)
- Integrates target expression with mRNA log2FC to assess directionality
- Visualisations: MA plots, volcano plot, sample distance/correlation heatmaps, dotplots, barplots, cnetplot, emapplot

**Input:** `counts.txt`, `samples.txt`  
**Output:** `exosomes/results/LFC{n}/{cond_test}{cond_ref}/` → `.xlsx` and `.pdf` files

---

### `02_free_miRNA_DE/free_miRNAs_DE_analysis_ashr.R`
Differential expression analysis of **circulating free miRNAs**. Same pipeline as script 01 but parameterised via command-line arguments:

```bash
Rscript free_miRNAs_DE_analysis_ashr.R <cond_ref> <cond_test> <LFC>
# e.g.:
Rscript free_miRNAs_DE_analysis_ashr.R g1 g2 1
```

**Input:** `counts.txt`, `samples.txt`  
**Output:** `results/LFC{n}/{cond_test}{cond_ref}/` → several `.xlsx` files

---

### `03_correlation_analysis/correlation_miRNA_mRNA_new_analysis.R`
Integrates exosomal and free miRNA DE results to identify miRNA–mRNA regulatory relationships.

- Finds DE miRNAs common across comparisons (ATTRv vs CC and ATTRv vs HC)
- Converts precursor miRNA IDs to mature forms via `miRBaseConverter`
- Retrieves validated miRNA–mRNA targets with `multiMiR`
- Computes Spearman correlations between miRNA and mRNA expression
- Performs functional enrichment: GO (BP, MF, CC) and KEGG via `clusterProfiler`
- Exports results and plots to `.xlsx` and `.pdf`

**Input:** `.xlsx` files from steps 01 and 02  
**Output:** `protein_coding_new_analysis/{free,exosomes}/` → `.xlsx` and `.pdf` files

---

## Dependencies

| Package | Purpose |
|---|---|
| `DESeq2` | Differential expression |
| `ashr` | Log2FC shrinkage |
| `biomaRt` | Gene annotation |
| `multiMiR` | miRNA target retrieval |
| `miRBaseConverter` | Precursor → mature miRNA conversion |
| `clusterProfiler` | GO and KEGG enrichment |
| `org.Hs.eg.db` | Human gene annotation |
| `ggplot2`, `tidyverse` | Data wrangling and visualisation |
| `openxlsx`, `readxl` | Excel I/O |

Install all dependencies from Bioconductor and CRAN:

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "ashr", "biomaRt", "multiMiR",
                        "miRBaseConverter", "clusterProfiler", "org.Hs.eg.db"))

install.packages(c("tidyverse", "ggplot2", "openxlsx", "readxl", "data.table"))
```

---

## Notes

- Before running, update the working directory paths in each script to match your local setup.
- Scripts 01 and 02 expect `counts.txt` and `samples.txt` in the working directory.
- Group labels (`g1`–`g4`) correspond to ATTRv, CC, and HC sample groups as defined in the original study design.
