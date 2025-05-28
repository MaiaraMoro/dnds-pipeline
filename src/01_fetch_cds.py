import requests
import os
import sys
from pathlib import Path
from utils.load_gene_file import load_genes_dict

ROOT_DIR = Path(__file__).resolve().parent.parent

# Ensembl database version: 113 

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


def save_as_fasta(gene_name, species_seq_dict, output_file):
    """ 
    Save the sequences in FASTA format.
    The file name will be <gene_name>.fasta and
    the sequences will be in the format:
    >species_name
    <sequence>
    \t

    """
    with open(output_file, 'w', encoding='utf-8') as f:
        for species, (seq,ensembl_id) in species_seq_dict.items():
            header = f"{ensembl_id}|{species}"
            f.write(f">{header}\n")
            f.write(f"{seq}\n")
    
if __name__ == "__main__":
    gene_symbol = sys.argv[1] # gene name
    out_path = sys.argv[2] # output path

    genes_file = ROOT_DIR / "data" / "genes" / "genes.tsv"
    species_file = ROOT_DIR / "data" / "species" / "mammals.txt"


    genes_dict = load_genes_dict(genes_file)
    if gene_symbol not in genes_dict:
        print(f"Gene {gene_symbol} not found in the gene list.")
        Path(out_path).write_text("")
        sys.exit(0)
    
    ensembl_id = genes_dict[gene_symbol]
    species_set = load_species(species_file)

    print(f"Fetching sequences for {gene_symbol}...")
    seqs = get_all_orthologs_cds(gene_symbol)

    if not seqs:
        print(f"Any 1:1 orthologs found for {gene_symbol}")
        Path(out_path).write_text("")
        sys.exit(0)
    
    if species_set:
        seqs = {sp: seq for sp, seq in seqs.items() if sp in species_set}
        if not seqs:
            print(f"Any orthologs found for {gene_symbol} in the species list")
            Path(out_path).write_text("")
            sys.exit(0)
            
    os.makedirs(Path(out_path).parent, exist_ok=True)
    save_as_fasta(gene_symbol, seqs, out_path)
    print(f"Saved {len(seqs)} sequences for {gene_symbol} to {out_path}")








