from Bio import AlignIO
import sys

input_file = sys.argv[1] 
output_file = sys.argv[2]

alignment = AlignIO.read(input_file, "fasta")

nseq = len(alignment)
length = alignment.get_alignment_length()

with open(output_file, "w") as out:
    # Cabeçalho
    out.write(f" {nseq} {length}\n")
    for record in alignment:
        # Ajustar ID (máx 50, espaçamento fixo com padding)
        id_fixed = record.id[:50].ljust(52)
        seq = str(record.seq).replace("-", "N")
        seq_no_wrap = ''.join(seq)
        out.write(f"{id_fixed}{seq_no_wrap}\n")
