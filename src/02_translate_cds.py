import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def translate_fasta(input_path, output_path):
    """
    Translate a FASTA file containing CDS sequences to protein sequences.
    The input FASTA file should contain sequences in the format:
    >species_name
    <sequence>
    
    The output FASTA file will contain the translated protein sequences.
    """
    records = []
    
    for record in SeqIO.parse(input_path, "fasta"):
        # Translate the CDS sequence to protein
        seq = record.seq
        
        if len(seq) % 3 != 0:
            print(f"Warning: Sequence length is not a multiple of 3 for {record.id}. Skipping translation.", file=sys.stderr)
            continue

        if "N" in seq:
            print(f"Warning: Sequence contains 'N' for {record.id}.", file=sys.stderr)
        
        try: 
            protein = seq.translate(table=1, cds=True)
        
        except Exception as e:
            print(f"Error translating sequence {record.id}: {e}", file=sys.stderr)
            continue
    
        records.append(SeqRecord(protein, id=record.id, description=""))

    SeqIO.write(records, output_path, "fasta")

if __name__ == "__main__":
    input_path = sys.argv[1]  # Input FASTA file
    output_path = sys.argv[2]  # Output FASTA file
    translate_fasta(input_path, output_path)