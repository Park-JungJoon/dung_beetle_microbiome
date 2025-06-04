# Dung Beetle Gut Microbiome Analysis

This repository contains the code and data used to analyze gut microbiome profiles of various dung beetle species collected from contrasting ecological environments, specifically **grasslands** and **deserts**. The study aims to understand how environmental factors influence microbial community structure and identify key microbial taxa that differentiate habitats.

---

## Background

Dung beetles (family: Scarabaeidae) play a crucial role in nutrient cycling and ecosystem functioning. Their gut microbiome is shaped by host species, diet, and environmental conditions. However, the influence of ecological context (e.g., desert vs. grassland) on gut microbiome composition remains underexplored.

This study addresses that gap by comparing gut microbiomes across different beetle species and habitats using both **compositional-aware** and **count-based** differential abundance analyses.

---

## Objectives

- Assess **beta diversity** across habitat types using Bray–Curtis distances and **PERMANOVA**
- Identify **differentially abundant sequence variants (SVs)** between desert and grassland environments
- Evaluate both **statistical significance** and **biological relevance** of microbial changes
- Compare two complementary approaches: **DESeq2** (count-based) and **ALDEx2** (composition-aware)

---

## Analyzed Beetle Species

- `Onthophagus nuchicornis`
- `Bodilus sordescens`
- `Onthophagus gibbulus`
- `Gymnopleurus mopsus`

Each species has associated data files that include SV count matrices and sample metadata.
- `Onthophagus laticornis` species were excluded from the analysis after preprocessing because the data was too small, with one Desert sample in six Grassland samples.
---

## Repository Structure
```
dung_beetle_microbiome/
│
├── calculate_1/
│ ├── Python scripts for PERMANOVA analysis
│ └── Outputs: JSON files with test statistics
│
├── calculate_3/
│ ├── R scripts for DESeq2 and ALDEx2 analysis
│ └── Outputs: TSV files of differential abundance results
│
├── output/
│ └── Intermediate processed data and SV matrices
```
---
## Methods

### 1. Data Preprocessing

- Raw `.tsv` files include SV counts and sample class labels in the first row (`Class`)
- Samples with zero total SV counts are excluded
- SVs are **filtered by abundance**: only the **top 10% SVs by mean relative abundance** across samples are retained for PERMANOVA
- All SV counts are **normalized to relative abundance** for Bray–Curtis calculations

---

### 2. Beta Diversity Analysis (PERMANOVA)

- Implemented using **Python (scikit-bio)**
- Bray–Curtis distance matrices are calculated from relative abundance tables
- PERMANOVA tests for significant compositional differences by `SiteType` (desert vs. grassland)
- **999 permutations** used for robust p-value estimation


---

### 3. Differential Abundance Analysis
#### ➤ DESeq2 (Count-based method)
- Conducted using R (DESeq2 via phyloseq)
- Size factors estimated using the poscounts method (robust for sparse microbiome data)
- Log2 fold changes and adjusted p-values computed using Wald tests
- SVs with adjusted `padj < 0.05` and `abs(Log2FC(desert/grassland)) > 1` are considered statistically significant

#### ➤ ALDEx2 (Composition-aware method)

- Conducted using **R (ALDEx2)**
- SVs are CLR-transformed via Monte Carlo Dirichlet sampling (128 iterations)
- Both Welch’s t-test and Wilcoxon test are applied
- SVs with **`we.eBH` or `wi.eBH` < 0.05** and **`|effect| ≥ 1`** are interpreted as significantly different


## Output
Each species yields:

*_permanova.json: PERMANOVA test results

*_DESeq2_result.tsv: Differentially abundant SVs (DESeq2)

*_ALDEx2_all.tsv: All ALDEx2 results

*_ALDEx2_significant.tsv: Filtered significant ALDEx2 results


## Requirements
Python

-  python >= 3.8

-  pandas v.2.2.3

-  scikit-bio v.1.6.1


R (for DESeq2 / ALDEx2)

-  R >= 4.0

-  DESeq2 v.1.48.1

-  phyloseq v.1.52.0

-  ALDEx2 v.1.40.0
-  ggplot2 v.3.5.2
-  vegan v.2.6-10
