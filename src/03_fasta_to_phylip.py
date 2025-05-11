from Bio import AlignIO
import sys

input_file = sys.argv[1] 
output_file = sys.argv[2]

alignment = AlignIO.read(input_file, "fasta")
AlignIO.write(alignment, output_file, "phylip")
print(f"Converted {input_file} to {output_file} in PHYLIP format.")