import csv

def load_genes_dict(genes_file:str) -> dict:
    """ 
    Load gene names from a file.
    Return a dictionary with gene names as keys and their Ensembl IDs as values
    """

    genes_dict = {}
    with open(genes_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene_name = row['Gene name'].strip()
            ensembl_id = row['Gene stable ID'].strip()
            genes_dict[gene_name] = ensembl_id 
    
    return genes_dict


def load_gene_names(genes_file:str) -> list:
    """ 
    Load gene names from a file.
    Return a list of gene names
    """

    names = []

    with open(genes_file, 'r', encoding='utf-8') as f:

        reader = csv.DictReader(f, delimiter='\t')

        for row in reader:
            gene_name = row['Gene name'].strip()
            names.append(gene_name)

    return names