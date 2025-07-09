"""
Add a header to a tree file 
"""
import sys, os, re
from Bio import SeqIO

fasta_in, tree_in, tree_out = sys.argv[1:]

ntaxa = sum(1 for _ in SeqIO.parse(fasta_in, "fasta"))

with open(tree_in) as ih:
    newick = re.sub(r"\s+", "", ih.read().strip())

os.makedirs(os.path.dirname(tree_out), exist_ok=True)
with open(tree_out, "w") as oh:
    oh.write(f"{ntaxa} 1\n{newick}\n")