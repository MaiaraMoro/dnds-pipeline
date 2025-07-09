import pathlib, sys, os

PROJECT_ROOT = pathlib.Path.cwd()    
SRC_DIR = PROJECT_ROOT / "src"

sys.path.insert(0, str(SRC_DIR))       
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
# 2 - Translate CDS to protein sequences and filter
#     to keep only valid CDS (good quality)
#######################################################
rule translate_cds:
    input: "data/raw_cds/{gene}.fasta"
    output: "data/raw_prot/{gene}.fasta"
    conda: "envs/dnds.yaml"
    shell: "python src/02_translate_cds.py {input} {output}       "

rule filter_valid_cds: 
    input:
        cds = "data/raw_cds/{gene}.fasta",
        prot = "data/raw_prot/{gene}.fasta"
    output: "data/filtered_cds/{gene}.fasta"
    conda: "envs/dnds.yaml"
    shell: """
    seqkit grep -f <(seqkit seq -n -i {input.prot}) {input.cds} -o {output}
    """
######################################################
# 3 - Align protein sequences (PRANK)
######################################################
rule align_protein: 
    input: "data/raw_prot/{gene}.fasta"
    output: "data/align_prot/{gene}.fasta"
    conda: "envs/dnds.yaml"
    shell: """
    prank -d={input} -o={output} -protein -quiet
    mv {output}.best.fas {output}
    """

######################################################
# 4 - Back-translate aligned protein sequences --> codons
######################################################
rule pal2nal:
    input:
        protein = "data/align_prot/{gene}.fasta",
        cds = "data/filtered_cds/{gene}.fasta"
    output: "data/align_codon/{gene}.fasta"
    conda: "envs/dnds.yaml"
    shell: """pal2nal.pl {input.protein} {input.cds} -output fasta > {output}"""

rule reorder_codon_alignment:
    input: "data/align_codon/{gene}.fasta"
    output: "data/align_codon_reordered/{gene}.fasta"
    conda: "envs/dnds.yaml"
    shell: """
    python src/03_reorder_codon_alignment.py {input} {output}
    """

######################################################
# 5 - Maximum likelihood tree IQTREE
######################################################
rule build_tree:
    input: "data/align_codon_reordered/{gene}.fasta"
    output: "data/trees/{gene}.treefile"
    threads: THREADS_TOTAL
    conda: "envs/dnds.yaml"
    shell: """iqtree2 -s {input} -st CODON -m MFP -nt {threads} -pre data/trees/{wildcards.gene}"""

rule add_header_tree:
    input:
        fasta = "data/align_codon_reordered/{gene}.fasta",
        tree  = "data/trees/{gene}.treefile"
    output: "data/trees_hdr/{gene}.tree"
    conda:  "envs/dnds.yaml"
    shell: "python src/04_add_header_tree.py {input.fasta} {input.tree} {output}" 
        
######################################################
# 6 - Convert FASTA to PHYLIP
######################################################
rule fasta_to_phylip:
    input: "data/align_codon_reordered/{gene}.fasta"
    output: "data/align_codon_phylip/{gene}.phy"
    conda: "envs/dnds.yaml"
    shell: """python src/05_fasta_to_phylip.py {input} {output}"""

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
    shell: "python src/06_analyze_paml_output.py {wildcards.gene} {output.lrt} {output.beb}"
