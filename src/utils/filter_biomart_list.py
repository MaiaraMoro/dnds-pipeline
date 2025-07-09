import csv
import sys

# Caminhos fixos
BIOMART_FILE   = "/home/maiara/Documents/bioinfo-masters/analise_dn_ds/data/genes/genes.tsv"
GENE_LIST_FILE = "/home/maiara/Documents/bioinfo-masters/analise_dn_ds/pre-analysis/candidate_genes_expanded - selected_genes_analise_dn_ds_leading_edge.csv"
OUTPUT_FILE    = "/home/maiara/Documents/bioinfo-masters/analise_dn_ds/pre-analysis/filtered_genes.tsv"


def load_gene_list(path):
    """
    Lê lista de genes (símbolos) de um CSV: SYMBOL,...
    Retorna um set de símbolos.
    """
    genes = set()
    with open(path) as f:
        for line in f:
            parts = line.strip().split(',')
            if parts and parts[0]:
                genes.add(parts[0].strip())
    return genes


def filter_biomart(biomart_path, gene_list_path, output_path):
    """
    Mantém apenas linhas cujo Gene name esteja em gene_list e
    cuja coluna 'Ensembl Canonical' seja 1 ou 1.0.
    """
    genes = load_gene_list(gene_list_path)
    
    with open(biomart_path) as inp, open(output_path, 'w', newline='') as outp:
        reader = csv.DictReader(inp, delimiter='\t')
        writer = csv.DictWriter(outp, fieldnames=reader.fieldnames, delimiter='\t')
        writer.writeheader()
        for row in reader:
            symbol = row.get('Gene name', '').strip()
            canonical_str = row.get('Ensembl Canonical', '').strip()
            if not symbol or symbol not in genes:
                continue
            try:
                canonical_val = float(canonical_str)
            except ValueError:
                continue
            if canonical_val == 1.0:
                writer.writerow(row)

if __name__ == '__main__':
    filter_biomart(BIOMART_FILE, GENE_LIST_FILE, OUTPUT_FILE)
    print(f"Filtragem completa: {OUTPUT_FILE}")
