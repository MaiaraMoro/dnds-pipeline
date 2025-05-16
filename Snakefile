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
        expand("results/paml/{gene}/mlc", gene=GENE_IDS)

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
    shell: """mafft --auto --thread {threads} {input} > {output}"""

######################################################
# 4 - Back-translate aligned protein sequences to CDS
######################################################
rule pal2nal:
    input:
        protein = "data/align_prot/{gene}.fasta",
        cds = "data/raw_cds/{gene}.fasta"
    output: "data/align_codon/{gene}.fasta"
    conda: "envs/dnds.yaml"
    shell: """pal2nal.pl {input.protein} {input.cds} -output fasta -nogap > {output}"""

######################################################
# 5 - Maximum likelihood tree IQTREE
######################################################
rule build_tree:
    input: "data/align_codon/{gene}.fasta"
    output: "results/trees/{gene}.treefile"
    threads: THREADS_TOTAL
    conda: "envs/dnds.yaml"
    shell: """iqtree2 -s {input} -st CODON -m MFP -nt 2 -pre data/trees/{wildcards.gene}"""

######################################################
# 6 - Convert FASTA to PHYLIP
######################################################
rule fasta_to_phylip:
    input: "data/align_codon/{gene}.fasta"
    output: "data/align_codon_phylip/{gene}.phy"
    conda: "envs/dnds.yaml"
    shell: """python src/03_fasta_to_phylip.py {input} {output}"""

#######################################################
# 7 - PAML (Nssites = 0,1,2,7,8) for every gene 
#######################################################
rule ctl_global:
    input: 
        aln = "data/align_codon_phylip/{gene}.phy",
        tree = "data/trees/{gene}.treefile"
    output: 
        ctl = "results/paml/{gene}/codeml.ctl"
    run:
        os.makedirs(os.path.dirname(output.ctl), exist_ok=True)
        ctl.write(textwrap.dedent(f"""
                seqfile = {os.path.abspath(input.aln)}
                treefile = {os.path.abspath(input.tree)}
                outfile = mlc

                noisy = 3
                verbose = 1
                runmode = 0
                cleandata = 0

                seqtype = 1
                ndata = 1 #one
                CodonFreq = 2
                icode = 0
                clock = 0
                fix_kappa = 0
                kappa = 2
                fix_omega = 0
                omega = 0.5

                model = 0
                NSsites = 0 1 2 7 8 
            """))

rule run_codeml_global:
    input: "results/paml/{gene}/codeml.ctl"
    output: "results/paml/{gene}/mlc"
    conda:
        "envs/paml.yaml"
    shell:
        """
        cd results/paml/{wildcards.gene} && codeml codeml.ctl
        """

