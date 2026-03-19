# lncRNA Signatures of Temozolomide Response in Glioblastoma

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the analysis code for the manuscript:

> **Long Non-Coding RNA Signatures of Temozolomide Response in IDH-Wildtype Glioblastoma: Identification of H19 and LINC01936 as Independent Prognostic Factors**
>
> Jaeil Han et al. *(Neuro-Oncology Advances, submitted)*

We identified 16 differentially expressed lncRNAs in TMZ responders vs. non-responders using TCGA-GBM RNA-seq data, and characterized H19 and LINC01936 as independent prognostic factors, with DNM3OS as a ceRNA network hub.

---

## Data Availability

- **TCGA-GBM RNA-seq**: Available via [GDC Data Portal](https://portal.gdc.cancer.gov/) (project: TCGA-GBM)
- **CGGA validation cohort**: Available at [CGGA](http://www.cgga.org.cn/) (mRNAseq_693)
- **Reference genome**: GRCh38 + GENCODE v44

> Raw data cannot be redistributed per TCGA data use agreement. Processed count matrices are available upon reasonable request.

---

## Repository Structure

```
.
├── config/
│   └── analysis_config.yaml       # All analysis parameters
├── envs/
│   └── environment.yaml           # Conda environment specification
└── scripts/
    ├── 00_setup/
    │   ├── setup_environment.sh   # Environment setup
    │   └── download_reference.sh  # GRCh38 + GENCODE v44 download
    ├── 01_data_acquisition/
    │   ├── 01_query_gdc_metadata.py     # GDC API query
    │   ├── 02_download_molecular_annotations.py
    │   └── 03_download_fastq.sh         # TCGA FASTQ download
    ├── 02_preprocessing/
    │   ├── 01_fastqc_raw.sh             # Raw QC (FastQC + MultiQC)
    │   ├── 02_trim_galore.sh            # Adapter trimming
    │   ├── 03_star_alignment.sh         # STAR 2-pass alignment
    │   └── 04_alignment_qc.sh           # Alignment QC (RSeQC)
    ├── 03_quantification/
    │   ├── 01_featurecounts.sh          # featureCounts
    │   └── 02_separate_lncrna_mrna.R   # lncRNA/mRNA separation
    ├── 04_differential_expression/
    │   └── 01_deseq2_analysis.R         # DESeq2 + SVA
    ├── 05_functional_analysis/
    │   ├── 01_coexpression_analysis.R   # Spearman co-expression
    │   ├── 02_pathway_enrichment.R      # GO/KEGG/GSEA
    │   ├── 03_wgcna_network.R           # WGCNA
    │   ├── 04_cerna_network.R           # ceRNA network
    │   └── 04b_cerna_with_encori.R      # ceRNA with ENCORI data
    ├── 06_survival_analysis/
    │   └── 01_survival_analysis.R       # KM, Cox, LASSO, Bootstrap
    ├── 07_validation/
    │   ├── 01_external_validation.R     # External cohort validation
    │   └── 02_cgga_validation.R         # CGGA-specific validation
    ├── 08_figures/
    │   └── 01_generate_main_figures.R   # Publication figures (Fig 2–4)
    └── 09_final_integration/
        ├── 01_final_summary.R           # Final summary + Fig 5–7
        ├── 02_fig1_study_design.R       # Study design figure
        └── 03_prepare_manuscript_figures.R
```

---

## Requirements

### Conda Environment

```bash
conda env create -f envs/environment.yaml
conda activate lncrna-gbm-tmz
```

### Key Software Versions

| Tool | Version |
|------|---------|
| STAR | 2.7.10a |
| featureCounts | 2.0.3 |
| FastQC | 0.11.9 |
| Trim Galore | 0.6.7 |
| R | 4.3.x |
| DESeq2 | 1.42.x |
| clusterProfiler | 4.10.x |
| survival | 3.5.x |

---

## How to Reproduce

### 1. Setup

```bash
# Clone this repository
git clone https://github.com/jaeilhan-sketch/lncRNA-GBM-TMZ-code.git
cd lncRNA-GBM-TMZ-code

# Create conda environment
conda env create -f envs/environment.yaml
conda activate lncrna-gbm-tmz

# Download reference genome and annotations
bash scripts/00_setup/download_reference.sh
```

### 2. Data Acquisition

```bash
# Query GDC metadata (94 IDH-wt primary GBM, TMZ-treated)
python scripts/01_data_acquisition/01_query_gdc_metadata.py

# Download FASTQ files (requires GDC token)
bash scripts/01_data_acquisition/03_download_fastq.sh
```

### 3. Preprocessing & Alignment

```bash
bash scripts/02_preprocessing/01_fastqc_raw.sh
bash scripts/02_preprocessing/02_trim_galore.sh
bash scripts/02_preprocessing/03_star_alignment.sh   # STAR 2-pass
bash scripts/02_preprocessing/04_alignment_qc.sh
```

### 4. Quantification & Differential Expression

```bash
bash scripts/03_quantification/01_featurecounts.sh
Rscript scripts/03_quantification/02_separate_lncrna_mrna.R
Rscript scripts/04_differential_expression/01_deseq2_analysis.R
```

### 5. Functional & Network Analysis

```bash
Rscript scripts/05_functional_analysis/01_coexpression_analysis.R
Rscript scripts/05_functional_analysis/02_pathway_enrichment.R
Rscript scripts/05_functional_analysis/04b_cerna_with_encori.R
```

### 6. Survival Analysis

```bash
Rscript scripts/06_survival_analysis/01_survival_analysis.R
```

### 7. Validation & Figures

```bash
Rscript scripts/07_validation/02_cgga_validation.R
Rscript scripts/08_figures/01_generate_main_figures.R
Rscript scripts/09_final_integration/01_final_summary.R
```

---

## Key Analysis Parameters

See `config/analysis_config.yaml` for all parameters. Key settings:

- **TMZ response definition**: PFS ≥ 6 months (Responder) vs. < 6 months (Non-Responder)
- **Inclusion criteria**: IDH-wildtype, primary GBM, TMZ-treated (n = 94)
- **lncRNA filtering**: CPM > 0.5 in ≥ 30% of samples
- **DE threshold**: |log2FC| > 1, adjusted p-value < 0.05 (Benjamini-Hochberg)
- **Reference**: GRCh38 + GENCODE v44
- **Alignment**: STAR 2-pass mode

---

## Citation

If you use this code, please cite:

> Han J, et al. Long Non-Coding RNA Signatures of Temozolomide Response in IDH-Wildtype Glioblastoma. *Neuro-Oncology Advances*, 2025. (in press)

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

## Contact

Jaeil Han — Han Lab
For questions about the analysis pipeline, please open a [GitHub Issue](https://github.com/jaeilhan-sketch/lncRNA-GBM-TMZ-code/issues).
