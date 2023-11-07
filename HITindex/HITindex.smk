#!/usr/bin/env python
"""
HITindex Pipeline Snakemake Script

This script defines the Snakemake rules for running the HITindex pipeline on a set of samples.
It processes RNA-seq data to classify reads and calculate HITindex metrics.

Usage:
    snakemake --configfile config.yaml --snakefile HITindex.smk --rerun-incomplete

Note:
    - The config file should define paths and parameters required for the analysis, including:
        - gtf_path: Path to the GTF file.
        - output_dir: Directory where output files will be saved.
        - annotation_file_path: Path to the sample annotation file, with the format described here: https://gagneurlab-drop.readthedocs.io/en/latest/prepare.html#creating-the-sample-annotation-table .
        - tissue: Specific tissue(s) to analyze or None to analyze all tissues.
    - The values provided in the 'config.yaml' and 'HIT_identity_parameters.txt' files represent the default/recommended settings for the HITindex pipeline. These were chosen based on the standard analysis requirements and can be modified to suit specific needs. 
    - The GTF file used must be sorted by gene -> transcript -> exons. If it's not sorted, use AGAT to sort it:
      [agat_convert_sp_gff2gtf.pl](https://github.com/NBISweden/AGAT)

Author: Raquel Romão
Date: November 2023
"""

import pandas as pd
import os


configfile: "config.yaml"

sample_annot = pd.read_csv(config['annotation_file_path'], sep='\t')


def get_output_files(tissues):
    """
    Generate a list of output file paths for the given tissues.
    If tissues is None, all unique tissues from the sample annotation datatable are used.
    """
    output_files = []
    # If tissues is a single string, convert it to a list
    if isinstance(tissues, str):
        tissues = [tissues]
    # If tissues is None, get all unique tissues from the sample annotation
    elif tissues is None:
        tissues = sample_annot["TISSUE"].unique()

    # Iterate over each tissue to get the output files
    for tissue in tissues:
        output_dir_tissue = config['output_dir'] + tissue + "/"
        filtered_samples = sample_annot[sample_annot["TISSUE"] == tissue]
        for _, row in filtered_samples.iterrows():
            sample_name = row["RNA_ID"]
            output_files.extend([
                output_dir_tissue + f"{sample_name}_HITindex.AFEPSI",
                output_dir_tissue + f"{sample_name}_HITindex.ALEPSI",
                output_dir_tissue + f"{sample_name}_HITindex.exon",
                output_dir_tissue + f"{sample_name}.sorted.junctions.bam",
                output_dir_tissue + f"{sample_name}.sorted.junctions.bam.bai",
                output_dir_tissue + f"{sample_name}.sorted.junctions.bam_header.txt",
                output_dir_tissue + f"{sample_name}.sorted.junctions_exonjuncs.bed.gz"
            ])
    return output_files


def get_input_bam(sample_name):
    """
    Retrieve the BAM file path for a given sample name.
    """
    for _, row in sample_annot.iterrows():
       if row["RNA_ID"] == sample_name:
           bam_path = row["RNA_BAM_FILE"]
           return bam_path
    print("Sample ID not found.")

def memory_mb(attempts):
    """
    Calculate the memory requirement in MB based on the number of snakemake attempts.
    """
    return 25000 + (5000 * attempts)



rule all:
    input:
        get_output_files(config.get('tissue', None))
         
rule run_HITindex_annotate:
    input:
        gtf_file = config['gtf_path']
    params:
        ss3buffer = config['ss3buffer'],
        ss5buffer = config['ss5buffer']
    output:
        out_file = config['output_dir'] + "metaexons.bed_ss3-50ss5-20.buffer"
    shell:
        "python HITindex_annotate.py --gtf {input.gtf_file} --ss3buffer {params.ss3buffer} --ss5buffer {params.ss5buffer} --outfile {output.out_file}"
         
rule run_HITindex_classify:
    input:
        bam_file = lambda wildcards: get_input_bam(wildcards.sample_name),
        bed = rules.run_HITindex_annotate.output.out_file
    params:
        hitindex_output_prefix = lambda wildcards: config['output_dir'] + wildcards.tissue + "/" + wildcards.sample_name + "_HITindex",
        readstrand = config['readstrand'],
        readtype = config['readtype'],
        overlap = config['overlap'],
        readnum = config['readnum'],
        bootstrap = config['bootstrap']
    output:
        junctions_bam = config['output_dir'] + "{tissue}/{sample_name}.sorted.junctions.bam",
        junctions_bam_bai = config['output_dir'] + "{tissue}/{sample_name}.sorted.junctions.bam.bai",
        junctions_bam_header = config['output_dir'] + "{tissue}/{sample_name}.sorted.junctions.bam_header.txt",
        junctions_exonjuncs_bed = config['output_dir'] + "{tissue}/{sample_name}.sorted.junctions_exonjuncs.bed.gz",
        hitindex_AFEPSI = config['output_dir'] + "{tissue}/{sample_name}_HITindex.AFEPSI",
        hitindex_ALEPSI = config['output_dir'] + "{tissue}/{sample_name}_HITindex.ALEPSI",
        hitindex_exon = config['output_dir'] + "{tissue}/{sample_name}_HITindex.exon"
    resources:
        mem_mb = lambda wildcards, attempt: memory_mb(attempt),
        threads = 1
    shell:
        """
        rm -f {output.junctions_bam}.*.bam
        python HITindex_classify.py --junctionReads --bam {input.bam_file} --juncbam {output.junctions_bam} --readstrand {params.readstrand} --readtype {params.readtype} --HITindex --bed {input.bed} --overlap {params.overlap} --readnum {params.readnum} --classify --calculatePSI --outname {params.hitindex_output_prefix} --parameters HIT_identity_parameters.txt --bootstrap {params.bootstrap}
        """
       #included cleanup step of eventual temporary files not deleted when a job fails
   