# Dual RNA-seq Analysis of *Glycine max* and *Alternaria alternata*
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.32.4-brightgreen.svg)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains a reproducible Snakemake pipeline for analyzing dual RNA-seq data from *Glycine max* (soybean) following infection by the fungal pathogen, *Alternaria alternata*. The project aims to identify defense responses in soybean and virulence mechanisms in the pathogen across three post-inoculation time points (1dpi, 7dpi, and 14dpi). This work is part of my MSc thesis research at Brandon University.

## 📊 Project Overview
*Alternaria alternata* is a fungal pathogen that causes leaf spot disease in soybean, resulting in considerable yield loss. Yet, the gene expression changes it induces in soybean are not completely understood. This project uses a multi-omics approach to unravel these plant-pathogen interactions.

**Key Goals:**
- Identify differentially expressed genes (DEGs) in soybean during infection.
- Characterize the expression of virulence factors and effector genes in *A. alternata*.
- Understand the metabolic reprogramming of both organisms during infection.

## 📁 Project Structure

```
soybean-alt-rnaseq/
├── config/
│ └── config.yml     # Project configuration files (e.g., `config.yml` with paths to reference data)
├── data/
│ └── raw/           # Place raw FASTQ files here
├── resources/       # Reference genomes and annotations
│ ├── genomes/
│ └── annotations/
├── workflow/        # Snakemake workflow definition
│ └── Snakefile
├── scripts/         # Helper scripts for analysis
├── notebooks/       # Jupyter notebooks for visualization
├── results/         # Pipeline output (DE lists, plots, etc.)
├── logs/            # Log files from pipeline execution
├── .gitignore
├── environment.yml
└── README.md
```

## ⚙️ Installation & Setup
1. **Snakemake** are best installed via the [Conda](https://docs.conda.io/en/latest/). It is recommended to install **conda** via **Miniforge**. To install run:
 ```bash
 conda create -c conda-forge -c bioconda -c nodefaults --name snakemake snakemake
 conda activate snakemake
 ```

2. **Install Git**: If you don't have git, install it first. For Debian-based distributions like Ubuntu, Mint use:
 ```bash
sudo apt-get update
sudo apt-get install git
```

3.  **Clone the repository:**
 ```bash
 git clone https://github.com/Achukwunta/soybean-alt-rnaseq.git
 cd soybean-alt-rnaseq
 ```

4. **Download Reference Data:**  
Reference genomes and annotations from Ensembl must be downloaded and placed in the `resources/` directory. The `config.yml` file expects them at specific paths. See the **Data Sources** section below.

5.  **Create and activate the Conda environment:**
Install any required dependencies from `environment.yml` if you're not using **conda** to install everything initially:
 ```bash
 conda env create -f environment.yml
 conda activate soybean-alt-rnaseq
 ```

## 🚀 Running the Pipeline
The pipeline is managed with Snakemake. You need to first install:
```bash
#run the full analysis with:
snakemake --cores all --use-conda

#run a specific rule, like FastQC:
snakemake fastqc_raw --cores 2 --use-conda
```
Use `--use-conda` to ensure all dependencies are managed in isolated environments


## 🔧 Configuration

Project-specific variables are set in `config/config.yml`. You can modify these values according to your needs.

**Example**:
```yaml
samples:
  - Control_T14_1
  - Control_T14_2
  - Control_T14_3
  - Alternaria_T14_1
  - Alternaria_T14_2

params:
  trimmomatic: "ILLUMINACLIP:resources/adapters/NexteraPE-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:5:18 MINLEN:30"

references:
  glygenome: "resources/genomes/Gmax_880_v6.0.fa"
  glyannotation: "resources/annotations/Gmax_880_Wm82.a6.v1.gene.gff3"


**params**:
trimmomatic: "ILLUMINACLIP:resources/adapters/NexteraPE-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:5:18 MINLEN:30"

**references**:
  glygenome: "resources/genomes/Gmax_880_v6.0.fa"
  glyannotation: "resources/annotations/Gmax_880_Wm82.a6.v1.gene.gff3"
```
- **Samples**: You can add new samples in this list, making sure they match the expected naming convention.
- **Params**: Adjust trimming parameters here if needed (for example, changing adapter files or quality filtering criteria).
- **References**: Make sure the genome and annotation files are correctly placed in the `resources/` folder. See the **Data Sources** section for more info.

## 📊 Data Sources
- **Raw Reads**: Generated in Cassone's Lab, Department of Biology, Brandon University
- **Reference Genomes & Annotations**: Downloaded from Phytozome. See the download link below.

**Glycine max (Soybean):**
- **Assembly**: Gmax_880_v6.0.fa  
- **Source**: [Phytozome v14](https://phytozome-next.jgi.doe.gov/info/Gmax_Wm82_a6_v1) 
- **Note**: The reference genome and annotation files were downloaded from Phytozome. Due to access restrictions, direct download links are not available. Users can freely download the files from the Phytozome portal.

## 👨‍💻 Author
**Augustine Chukwunta**  
MSc Biology, Brandon University, Manitoba, Canada
Specialization in Bioinformatics  
Thesis: *"Multi-Omics Approaches to Unravel Plant-Pathogen Interactions: Transcriptomics and Microbiome Analysis."*  
Advisor: Dr. Bryan Cassone  
[![GitHub: Achukwunta](https://img.shields.io/badge/GitHub-Achukwunta-blue?logo=github)](https://github.com/Achukwunta)

## 📜 Citation
This work is part of an ongoing study. If you use this pipeline or results, please cite:  

Chukwunta A., & Cassone B. (2025). *Dual RNA-seq analysis reveals metabolic reprogramming of Alternaria alternata and defense activation in Glycine max during infection.* (Target Journal: Molecular Plant-Microbe Interactions).

## 🤝 Contributing & Contact
For questions, suggestions, or collaboration, please contact Augustine Chukwunta or open an issue on this repository.

Thank you and Happy research!


