# AberrantTerminalExonDetection
Pipelines for detecting aberrant events in terminal exons of RNA-seq data using adapted versions of OUTRIDER [1] and FRASER [2]. Both of the adapted exon-level tools use as a starting point the HITindex [3] results.

##### References
[1] OUTRIDER: Brechtmann F, Mertes C, Matusevičiūtė A, et al. OUTRIDER: A Statistical Method for Detecting Aberrantly Expressed Genes in RNA Sequencing Data. Am J Hum Genet. 2018;103(6):907-917. https://doi.org/10.1016/j.ajhg.2018.10.025 | [GitHub](https://github.com/gagneurlab/drop/tree/master/drop/modules/aberrant-expression-pipeline)

[2] FRASER: Mertes, C., Scheller, I.F., Yépez, V.A. et al. Detection of aberrant splicing events in RNA-seq data using FRASER. Nat Commun 12, 529 (2021). https://doi.org/10.1038/s41467-020-20573-7 | [GitHub](https://github.com/gagneurlab/drop/tree/master/drop/modules/aberrant-splicing-pipeline)

[3] HITindex: Fiszbein A, McGurk M, Calvo Roitberg E, Kim GY, Burge CB, and Pai AA. (2021). Widespread occurrence of hybrid internal-terminal exons in human transcriptomes. (bioRxiv) doi: https://doi.org/10.1101/2021.05.27.446076 | [GitHub](https://github.com/thepailab/HITindex)





## 1. HITindex pipeline - Starting point for the detection of aberrant events at the exon-level

  HITindex results are the starting point of the exon-level tools implemented to detect aberrant events in terminal exons. A Snakemake-based workflow that employs the HITindex [3] pipeline was designed to process RNA-seq data to classify reads and calculate HITindex metrics, across samples. It is tailored for use with GTEx samples but can be adapted for other datasets.
  
  The HITindex is a pipeline to classify hybrid, internal, or terminal exons from RNA-seq data by modeling ratios of splice junction coverage. The pipeline involves two major steps:
   
  - HITindex_annotate: Annotate metaexons from a GTF file by collapsing overlapping constituent exons.

  - HITindex_classify: Calculate HIT index metrics and classify metaexons into one of 5 exon-types: first, first-internal, internal, internal-last, and last exons.
  
  
  #### Features
  - Classifies and quantifies hybrid, internal, or terminal exons from RNA-seq data by modeling ratios of splice junction coverage, across samples.
  - Calculates HITindex metrics such as the number of upstream and downstream reads.
  - Configurable for specific tissues or all tissues in a dataset.
  - Utilizes default HITindex parameters with the option for customization.
  
  #### Requirements
  - samtools (v1.3)
  - bedtools (v2.26.0)
  - python (v3.6.3)
  
  #### Dependencies
  - scipy (v1.5.2)
  - numpy (v1.19.2)
  - pysam (v0.16)
  - pybedtools (v0.8.1)
  - pandas (v0.25.3)
  - pymc3 (v3.9.3)
  
  #### Installation
  Clone the repository and install the required Python packages:
  
    git clone https://github.com/raquel-romao/AberrantTerminalExonDetection.git
    cd HITindex
    pip install -r requirements.txt
  
  #### Configuration
  
  Edit the config.yaml file to specify the paths to your input files and set the desired HITindex parameters for the analysis in the config.yaml file and the HIT_identity_parameters.txt file
  
  #### Usage
  
  To run the pipeline, use the following command:
  
    snakemake --configfile config.yaml --snakefile HITindex.smk --rerun-incomplete
  
  
  #### Input Data Format
  
  The pipeline expects the following input data:
  - A GTF file sorted by gene -> transcript -> exons. Use [AGAT](https://github.com/NBISweden/AGAT) (agat_convert_sp_gff2gtf.pl) if sorting is required.
  - A sample annotation file in TSV format (with a similar format to what is described [here](https://gagneurlab-drop.readthedocs.io/en/latest/prepare.html#creating-the-sample-annotation-table)).
  
  #### Output Description
  
  The pipeline generates the same outputs described in the [HITindex GitHub](https://github.com/thepailab/HITindex).


## 2. Exon-level OUTRIDER

This repository contains the exon-level adaptation of the OUTRIDER pipeline, which is used to detect aberrant expression events in terminal exons. The pipeline is an extension of the gene-level aberrant expression module from DROP, adapted to work with exon-level data derived from HITindex results.

### Overview

The exon-level OUTRIDER pipeline follows a two-part structure:

1. **Data Preparation**: Instead of using BAM files, this pipeline computes exon expression outliers from HITindex results. A count matrix is generated from the HIT result tables across samples, in particular the upstream and downstream read counts per exon, which is then used to create an OUTRIDER dataset.
2. **OUTRIDER Fitting**: The pipeline determines the optimal encoding dimension, fits the OUTRIDER model, computes p-values, and extracts results at the exon-level.


![adaptedOUTRIDER](https://github.com/raquel-romao/AberrantTerminalExonDetection/assets/92799656/ed7d9372-3aa9-4e36-a67d-1c1236757a3d)

*Figure 1: Adaptation of OUTRIDER for the detection of aberrant expression outliers at the exon-level.*

### Requirements

- R (version specified)
- OUTRIDER (version specified)
- Other dependencies listed in `requirements.txt`

### Installation

Clone the repository and install the required R packages:

    git clone https://github.com/raquel-romao/AberrantTerminalExonDetection.git
    cd exon_level_OUTRIDER
    Rscript -e 'install.packages(c("required", "packages"))'

### Configuration

Edit wbuild.yaml to adjust the parameters for the analysis. Paths must be set in the YAML headers of the scripts.

### Usage

To run the pipeline, execute the following command:

    snakemake --jobs 1 --cluster "sbatch -N 1 -n 1 -c {threads} --mem {resources.mem_mb}M --parsable"

  Note: Cluster parameters and job management may vary based on your computing environment.
  
### Input Data Format

The pipeline expects input data in the form of HITindex result tables, which should contain the variables 'nUP' and 'nDOWN' for each sample.

### Output Description

The pipeline outputs text files with the results of the aberrant expression analysis, similar to the [original tool](https://github.com/gagneurlab/outrider) but at the exon-level (instead of the gene-level).

### Customization

Parameters for the analysis can be adjusted in the wbuild.yaml file. This includes setting thresholds for FPKM filtering, FDR, and Z-score cutoffs.

## 3. Exon-level FRASER

