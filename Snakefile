import pathlib, sys, os

PROJECT_ROOT = pathlib.Path.cwd()      # diretório onde o Snakemake foi iniciado
SRC_DIR = PROJECT_ROOT / "src"

sys.path.insert(0, str(SRC_DIR))       # põe no topo do sys.path
os.environ["PYTHONPATH"] = (
    f"{SRC_DIR}{os.pathsep}{os.environ.get('PYTHONPATH', '')}"
)

from utils.load_gene_file import load_gene_names

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
        expand("results/paml/{gene}/mlc",      gene=GENE_IDS),
        expand("results/paml/{gene}/lrt.tsv",  gene=GENE_IDS),
        expand("results/paml/{gene}/beb.tsv",  gene=GENE_IDS)

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
    output: "data/trees/{gene}.treefile"
    threads: THREADS_TOTAL
    conda: "envs/dnds.yaml"
    shell: """iqtree2 -s {input} -st CODON -m MFP -nt {threads} -pre data/trees/{wildcards.gene}"""

rule add_header_tree:
    input:
        fasta = "data/align_codon/{gene}.fasta",
        tree  = "data/trees/{gene}.treefile"
    output: "data/trees_hdr/{gene}.tree"
    conda:  "envs/dnds.yaml"
    run:
        import os, re
        from Bio import SeqIO

        # Nº de táxons (= nº de sequências no alinhamento)
        ntaxa = len(list(SeqIO.parse(input.fasta, "fasta")))

        # Newick sem quebras/brancos
        with open(input.tree) as ih:
            newick = re.sub(r"\s+", "", ih.read().strip())

        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        with open(output[0], "w") as oh:
            oh.write(f"{ntaxa} 1\n{newick}\n") 

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
        tree = "data/trees_hdr/{gene}.tree"
    output: 
        ctl = "results/paml/{gene}/codeml.ctl"
    run:
        os.makedirs(os.path.dirname(output.ctl), exist_ok=True)
        with open(output.ctl, "w") as f:
            f.write(textwrap.dedent(f"""
                seqfile = {os.path.abspath(input.aln)}
                treefile = {os.path.abspath(input.tree)}
                outfile = mlc

                noisy = 3
                verbose = 1
                runmode = 0
                cleandata = 1
                fix_blength = 2

                seqtype = 1
                ndata = 1 
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
        "envs/dnds.yaml"
    shell:
        """
        cd results/paml/{wildcards.gene} && codeml codeml.ctl
        """

rule analyze_paml_results:
    input: "results/paml/{gene}/mlc"
    output:
        lrt = "results/paml/{gene}/lrt.tsv",
        beb = "results/paml/{gene}/beb.tsv"
    conda: "envs/dnds.yaml"
    shell: "python src/04_analyze_paml_output.py {wildcards.gene} {output.lrt} {output.beb}"
