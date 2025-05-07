import requests
import os
import csv
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent.parent

# This script fetches CDS sequences for a list of genes from Ensembl REST API
# and saves them in FASTA format. It also filters the sequences based on a
# list of species provided in a separate file. 

# Ensembl database version: 113 

def load_genes(genes_file):
    """ 
    Load genes from a file and return a list of tuples (ensembl_id, gene_name).

    The file should contain one gene per line in the format:
    <Gene name>	<Ensembl Canonical>	<Gene stable ID>

    Third column is human ensembl gene ID.
    """

    genes = []

    with open(genes_file, 'r', encoding='utf-8') as f:

        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene_name = row['Gene name'].strip()
            ensembl_id = row['Gene stable ID'].strip()
            genes.append((ensembl_id, gene_name))

    return genes

def load_species(species_file):
    """
    Read a file containing species names of interest and return a set to 
    facilitate filtering fastas after downloading.
    The file should contain one species name per line.
    
    """

    species_set = set()

    with open(species_file, 'r', encoding='utf-8') as f:
        for line in f:
            sp = line.strip()
            if sp:
                species_set.add(sp)
    return species_set

def get_all_orthologs_cds(symbol,
                      species="homo_sapiens",
                      target_taxon=40674,
                      remove_gaps=False):
    """
    Get all orthologs for a given human gene ID from Ensembl REST API.
    The focus is on the CDS sequences.
    Then, filter orthologs 1:1 
    Return a dictionary species: sequence
    """
    # base url 
    url = (f"https://rest.ensembl.org/homology/symbol/{species}/{symbol}"
           f"?type=orthologues;target_taxon={target_taxon};sequence=cds")
    headers = {"Accept": "application/json"}
    
    r = requests.get(url, headers=headers)

    if not r.ok:
        print("Erro HTTP:", r.status_code, r.text)
        return []
    
    data = r.json()
    if "data" not in data or not data['data']:
        print(f"Nenhum dado retornado para gene {symbol}")
        return []
    
    homologies = data['data'][0]['homologies']

    one2one = [h for h in homologies if h["type"] == "ortholog_one2one"]

    if not one2one:
        print(f"Nenhum ortólogo 1:1 para gene {symbol}")
        return {}
    
    results = {}

    source_info = homologies[0]["source"]

    if source_info.get("species") == "homo_sapiens":
        seq_hs_with_gaps = source_info["align_seq"]
        if remove_gaps:
            seq_hs = seq_hs_with_gaps.replace("-", "")
        else:
            seq_hs = seq_hs_with_gaps
    
        ensembl_id_hs = source_info.get("id")

        results["homo_sapiens"] = (seq_hs, ensembl_id_hs)

    for homology in one2one:

        target_info = homology['target']
        sp_name = target_info['species']
        seq_with_gaps = target_info["align_seq"]
        ensembl_id = target_info["id"]

        if remove_gaps: # remove gaps from the sequence
            seq = seq_with_gaps.replace("-", "")
        else:
            seq = seq_with_gaps

        results[sp_name] = (seq, ensembl_id)

    
    return results


def save_as_fasta(gene_name, species_seq_dict, outdir = "fasta_output"):
    """ 
    Save the sequences in FASTA format.
    The file name will be <gene_name>.fasta and
    the sequences will be in the format:
    >species_name
    <sequence>
    \t

    """

    os.makedirs(outdir, exist_ok=True)

    output_file = os.path.join(outdir, f"{gene_name}.fasta")

    with open(output_file, 'w', encoding='utf-8') as f:
        for species, (seq,ensembl_id) in species_seq_dict.items():
            header = f"{ensembl_id}|{species}"
            f.write(f">{header}\n")
            f.write(f"{seq}\n")
    
    print(f"Saved {len(species_seq_dict)} sequences for {gene_name} to {output_file}")
    
def main(): 
    genes_file = ROOT_DIR / "data" / "genes" / "genes.tsv" # Path to the file containing genes
    species_file = ROOT_DIR / "data" / "species" / "mammals.txt" # Path to the file containing species
    
    gene_list = load_genes(genes_file)

    if os.path.exists(species_file):
        species_set = load_species(species_file)
        print(f"Loaded {len(species_set)} species from {species_file}")
    else:
        species_set = None
        print(f"Species file {species_file} not found. Skipping species filtering.")
    
    outdir = ROOT_DIR / "data" / "raw_cds"
    os.makedirs(outdir, exist_ok=True)


    # for each gene, get the orthologs and save them in FASTA format
    for ensembl_id, gene_name in gene_list:
        print(f"Processing gene {gene_name} with Ensembl ID {ensembl_id}")
        orthologs = get_all_orthologs_cds(gene_name)
        
        if not orthologs:
            print(f"No orthologs found for {gene_name}")
            continue
        if species_set:
            # Filter orthologs by species
            orthologs = {sp: seq for sp, seq in orthologs.items() if sp in species_set}
            print(f"Filtered to {len(orthologs)} orthologs for {gene_name} based on species set")
        else:
            print(f"No desired species found for {gene_name}")
            continue

        # Save the orthologs in FASTA format
        save_as_fasta(gene_name, orthologs, outdir)

if __name__ == "__main__":
    main()







