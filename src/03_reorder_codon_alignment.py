"""
Script para reordenar um arquivo de alinhamento FASTA, garantindo que
a sequência de Homo sapiens seja a primeira da lista.
"""
import sys
from Bio import SeqIO

def reorder_fasta_for_paml(input_path: str, output_path: str, ref_species: str = "homo_sapiens"):
    """
    Lê um FASTA, encontra a sequência da espécie de referência e a escreve primeiro,
    seguida pelas outras.
    """
    try:
        records = list(SeqIO.parse(input_path, "fasta"))
    except FileNotFoundError:
        print(f"ERRO: Arquivo de entrada não encontrado: {input_path}", file=sys.stderr)
        sys.exit(1)

    ref_record = None
    other_records = []

    # Separa o registro de referência dos outros
    for record in records:
        # Procura pelo nome da espécie no cabeçalho (ID)
        if ref_species in record.id:
            ref_record = record
        else:
            other_records.append(record)

    with open(output_path, "w") as f:
        if ref_record:
            # Escreve o registro de referência primeiro
            SeqIO.write(ref_record, f, "fasta")
            # Escreve os registros restantes
            SeqIO.write(other_records, f, "fasta")
            print(f"Arquivo reordenado com '{ref_species}' como primeira sequência em '{output_path}'.")
        else:
            # Se a referência não for encontrada, apenas copia o arquivo original
            SeqIO.write(records, f, "fasta")
            print(f"AVISO: Espécie de referência '{ref_species}' não encontrada. A ordem original foi mantida.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(f"Uso: python {sys.argv[0]} <input.fasta> <output.fasta>")

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    reorder_fasta_for_paml(input_file, output_file)