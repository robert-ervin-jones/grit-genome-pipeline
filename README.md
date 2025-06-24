# Genome Assembly and Annotation Pipeline

This repository contains a Snakemake-based workflow for long-read genome assembly, annotation, and quality assessment. It is designed for Oxford Nanopore data and supports both eukaryotic and prokaryotic genome analysis.

## Features
- Read cleaning and filtering
- Assembly with Flye
- Eukaryote/prokaryote separation with Whokaryote
- Gene prediction with BRAKER
- Functional annotation with eggNOG-mapper
- Quality assessment with BUSCO
- Automated environment management with Conda
- Wrapper script for easy pipeline execution

## Directory Structure
- `data/` — Raw input FASTQ files (not tracked by git)
- `results/` — All pipeline outputs (not tracked by git)
- `workflow/` — Snakemake workflow, environments, and scripts
- `resources/` — Reference data and configuration files
 

## Quick Start
1. **Install dependencies:**
   - [Mamba](https://github.com/conda-forge/miniforge)
   - [Apptainer](https://apptainer.org/docs/user/main/quick_start.html)
2. **Clone this repository:**
   ```sh
   git clone <repo-url>
   cd genome-pipeline
   ```
3. **If this is a fresh mamba installation, add bioconda to `.condarc`:**
   ```sh
   conda config --add channels bioconda
   conda config --add channels conda-forge
   ```
4. **Setup environment for resource downloads:**
   ```sh
   mamba create -n bakta bakta=1.11.0
   mamba create -n eggnog eggnog-mapper=2.1.13
   ```
5. **Download resources:**
   ```sh
   conda activate bakta
   bakta_db download --output resources/bakta-light --type light
   conda deactivate
   git clone https://github.com/nextgenusfs/augustus.git
   cp -r augustus/config resources/augustus_config
   rm -rf augustus
   conda activate eggnog
   download_eggnog_data.py --data_dir resources/eggnog_dbs
   conda deactivate
   ```
6. **Add your FASTQ files to `data/`**
    - Place your raw FASTQ files in the `data/` directory. The files should be named as `{sample}.fastq.gz`.
7. **Make conda environment:**
   ```sh
   mamba create -n snakemake snakemake=9.6.0 biopython=1.85 r=4.4 r-ggplot2=3.5.2 r-ggpubr=0.6.0
   ```
8. **Run the workflow using the wrapper script:**
   ```sh
   conda activate snakemake
   python run_snakemake.py --threads 8 --q_filter <N> --lengths <N>
   ```
   - Use `-n` or `--dry_run` for a dry run.
   - See `python run_snakemake.py -h` for all options.

## Output
Final results for each sample are collected in `results/final/{sample}/` and include:
- Assembly files
- Gene predictions
- BUSCO summary
- eggNOG annotations
- Plots



## Citation
If you use this pipeline, please cite the relevant tools and this repository.

TODO: Add citations for all the tools used in the pipeline.