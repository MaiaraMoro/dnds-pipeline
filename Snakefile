    #######################################
    # Snakefile for Dn/Ds Analysis        #
    #                                     #
    # This Snakefile orchestrates the     #
    # workflow for calculating Dn/Ds      #
    # ratios in gene sequences. It        #
    # includes configuration, parameters, #
    # and rules for each step of the      #
    # analysis pipeline.                  #
    #######################################

import yaml, os, textwrap
import sys
from pathlib import Path
from utils.load_gene_file import load_gene_names

configfile: "config.yaml"

THREADS_MAFFT = int(config['threads_mafft'])
THREADS_TOTAL = int(config['threads_total'])
GENES_TSV = config['genes_file']
GENE_IDS = load_gene_names(GENES_TSV)


    ########################################
    #               RULES                  #
    ########################################


rule all:
    input:
        expand("results/dnds/{gene}.txt", gene=GENE_IDS),
        "results/dnds/summary.tsv"


#######################################################
# 1 - Fetch CDS sequences
#######################################################
rule fetch_cds:
    output: "data/raw_cds/{gene}.fasta"
    threads: 1
    conda: "envs/dnds.yaml"
    shell: "python src/01_fetch_cds.py {wildcards.gene} {output}"

#######################################################
# 2 - Translate CDS to protein sequences
#######################################################
rule translate_cds:
    input: "data/raw_cds/{gene}.fasta"
    output: "data/raw_prot/{gene}.fasta"
    conda: "envs/dnds.yaml"
    shell: "python src/02_translate_cds.py {input} {output}       "

######################################################
# 3 - Align protein sequences (MAFFT)
######################################################
rule align_protein: 
    input: "data/raw_prot/{gene}.fasta"
    output: "data/align_prot/{gene}.fasta"
    threads: THREADS_MAFFT
    conda: "envs/dnds.yaml"
    shell: """mafft --auto --thread {threads} {input} > {output}""".format(threads=THREADS_MAFFT, input=input, output=output)

######################################################
# 4 - Back-translate aligned protein sequences to CDS
######################################################
rule pal2nal:
    input:
        protein = "data/align_prot/{gene}.fasta",
        cds = "data/raw_cds/{gene}.fasta"
    output: "data/align_codon/{gene}.fasta"
    conda: "envs/dnds.yaml"
    shell: """pal2nal.pl {input.protein} {input.cds} -output fasta -nogap > {output}""".format(input=input, output=output)

######################################################
