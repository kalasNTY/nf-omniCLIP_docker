# omniCLIP Nextflow Pipeline (Docker Version)

## Overview
This repository contains a Nextflow pipeline for running **omniCLIP**, a computational framework for analyzing **CLIP and ARTR-seq** NGS data interogating transcriptome-wide RNA-protein interaction. The pipeline is built to work with a **Docker container** and is running a **test dataset on chromosome 3** by default unless user provides custom data and parameters.

> **Note:** This pipeline does **not** support multi-core processing yet (`params.nb_cores = 1`). A version using a **Singularity container** for HPC environments will be published in a separate repository (**TBA**).

## Reference
This pipeline is based on the original **omniCLIP** framework. Please cite the following resources:

- **omniCLIP GitHub**: [philippdre/omniCLIP](https://github.com/philippdre/omniCLIP)
- **omniCLIP Publication**: [Genome Biology (2023)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1521-2)
- **Original ARTR-seq Data**: [Nature Methods (2023)](https://www.nature.com/articles/s41592-023-02146-w)

## Container Information
This pipeline uses my **Docker container** available on Docker Hub:

- **DockerHub Repository**: [kalasnty/omniclip_container](https://hub.docker.com/r/kalasnty/omniclip_container/tags/)

### **Limitations**
- This pipeline is **not optimized for HPC environments** due to the use of Docker instead of Singularity.
- Multi-core support is currently **disabled** (`params.nb_cores = 1`).
- Results will be stored in the `Results` directory **by default** and will be overwritten if not saved separately.
- Bam files for your custom analysis need to be copied into the **Input** folder, not other folder naming is permitted.
- The pipeline is set-up to run a test analysis by default, to change that please adapt parameters and metadata.csv accordingly.

## Quick Start
1. **Clone the repository:**
   ```bash
   git clone https://github.com/kalasNTY/nf-omniCLIP_docker.git
   cd nf-omniCLIP_docker
   ```

2. **Run the pipeline with default parameters (chromosome 3 test dataset):**
   ```bash
   nextflow run omniclip_v4.nf
   ```

3. **To run on custom data, ensure that:**
   - You provide per-chromosome FASTA files in a **folder**.
   - The genome annotation matches the BAM file annotation (**UCSC vs. NCBI**).

## Output
- The default output directory is `Results/`.
- The contents of `Results/` will be **overwritten** in subsequent runs unless manually saved.
- The contents of Results/ need to be moved and saved to a different location prior to the next execution
- Please manually delete all unneccessary temp files (e.g. _dat and .db) unless you want them to be re-used.

## Contribution & Development
This pipeline was developed, debugged, and optimized by **Julian A. Zagalak**.

For future improvements, including **multi-core support and HPC compatibility**, check for updates on this repository.

---

If you encounter any issues, please report them via GitHub Issues.

