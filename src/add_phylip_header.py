import sys, re

treefile, ntaxa = sys.argv[1], int(sys.argv[2])
txt = open(treefile).read().strip()
with open(treefile, "w") as f:
    f.write(f"{ntaxa} 1\n{re.sub(r'\\s+', '', txt)}\n")