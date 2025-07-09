import requests
import os
import sys
import time
from pathlib import Path
from utils.load_gene_file import load_genes_dict
from Bio.Seq import Seq 

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

def get_all_orthologs_ids(symbol, species_filter_set,
                          species="homo_sapiens",
                          target_taxon=40674):
    """
    PASSO 1: Busca apenas os IDs estáveis dos ortólogos 1-para-1 no Ensembl.
    Retorna um dicionário {species_name: ensembl_id}.
    """
    print(f"Buscando IDs de ortólogos 1:1 para '{symbol}'...")
    url = (f"https://rest.ensembl.org/homology/symbol/{species}/{symbol}"
           f"?type=orthologues;target_taxon={target_taxon}")
    headers = {"Content-Type": "application/json"}
    
    try:
        r = requests.get(url, headers=headers, timeout=30)
        r.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Erro de rede ao buscar IDs: {e}")
        return None
    
    data = r.json()
    if "data" not in data or not data['data']:
        print(f"Nenhum dado de homologia retornado para o gene {symbol}")
        return None
    
    homologies = data['data'][0]['homologies']
    one2one = [h for h in homologies if h["type"] == "ortholog_one2one"]

    if not one2one:
        print(f"Nenhum ortólogo 1:1 encontrado para o gene {symbol}")
        return {}
    
    ortholog_ids = {}
    
    # Adiciona o ID da espécie de origem (humano)
    source_info = homologies[0].get("source")
    if source_info and source_info.get("id"):
        ortholog_ids[source_info.get("species")] = source_info.get("id")

    # Adiciona os IDs das espécies alvo
    for homology in one2one:
        target_info = homology['target']
        sp_name = target_info['species']

        if sp_name in species_filter_set:
            ensembl_id = target_info.get("id")
            if ensembl_id:
                ortholog_ids[sp_name] = ensembl_id

    print(f"Encontrados {len(ortholog_ids)} IDs de ortólogos.")
    return ortholog_ids
  
def get_canonical_transcript_id(gene_id, session):
    """
    Função auxiliar para encontrar o ID do transcrito canônico para um dado ID de gene.
    """
    url = f"https://rest.ensembl.org/lookup/id/{gene_id}?expand=1"
    headers = {"Content-Type": "application/json"}
    try:
        r = session.get(url, headers=headers, timeout=60)
        r.raise_for_status()
        data = r.json()

        if 'Transcript' not in data:
            return None

        # Procura pelo transcrito marcado como canônico
        for transcript in data['Transcript']:
            if transcript.get('is_canonical') == 1:
                return transcript.get('id')

        # Fallback: Se NENHUM for canônico, retorna o primeiro da lista por padrão.
        # Uma heurística melhor seria pegar o mais longo, mas para simplificar:
        if data['Transcript']:
             return data['Transcript'][0].get('id')

    except requests.exceptions.RequestException as e:
        print(f"  └─ AVISO: Falha ao buscar transcritos para {gene_id}. Erro: {e}")
        return None
    return None


def fetch_canonical_cds_from_genes(species_id_dict):
    """
    PASSO 2 e 3: Para cada ID de GENE, encontra seu TRANSCRITO canônico e
    então busca a sequência de CDS daquele transcrito específico.
    """
    if not species_id_dict:
        return {}

    print("Buscando e validando CDS dos transcritos canônicos...")
    session = requests.Session()
    final_results = {}
    
    for species, gene_id in species_id_dict.items():
        print(f"Processando {species} ({gene_id})...")
        
        # Passo 2: Encontrar o ID do transcrito canônico
        transcript_id = get_canonical_transcript_id(gene_id, session)
        time.sleep(0.1) # Pausa

        if not transcript_id:
            print(f"  └─ AVISO: Não foi possível encontrar um transcrito canônico para o gene {gene_id}. Pulando espécie.")
            continue
        
        print(f"  └─ Transcrito canônico encontrado: {transcript_id}")

        # Passo 3: Buscar o CDS para o ID do transcrito específico
        seq_url = f"https://rest.ensembl.org/sequence/id/{transcript_id}?type=cds"
        headers = {"Content-Type": "application/json"}
        
        try: 
            r = session.get(seq_url, headers=headers, timeout=60)
            r.raise_for_status()
            data = r.json()
            seq = data.get("seq")

            if not seq:
                print(f"    └─ AVISO: Nenhuma sequência encontrada para o transcrito {transcript_id}. Pulando.")
                continue

            # Validação (mesma de antes)
            if len(seq) % 3 != 0:
                print(f"    └─ AVISO: CDS do transcrito {transcript_id} não é múltiplo de 3. Pulando.")
                continue
            
            protein = Seq(seq).translate()
            if "*" in protein[:-1]:
                print(f"    └─ AVISO: CDS do transcrito {transcript_id} tem stop códon interno. Pulando.")
                continue
        
            final_results[species] = (seq, gene_id) # Salva com o ID do gene original para consistência

        except requests.exceptions.RequestException as e:
            print(f"    └─ AVISO: Falha ao buscar CDS para o transcrito {transcript_id}. Erro: {e}. Pulando.")
            continue
        
        time.sleep(0.1)

    return final_results


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
        if "homo_sapiens" in species_seq_dict:
            seq, ensembl_id = species_seq_dict["homo_sapiens"]
            header = f"{ensembl_id}|homo_sapiens"
            f.write(f">{header}\n")
            f.write(f"{seq}\n")
            
            # Remove a entrada humana do dicionário para não ser escrita novamente
            del species_seq_dict["homo_sapiens"]
            
        for species, (seq,ensembl_id) in species_seq_dict.items():
            header = f"{ensembl_id}|{species}"
            f.write(f">{header}\n")
            f.write(f"{seq}\n")
    
if __name__ == "__main__":
    # (bloco if __name__ == "__main__" ajustado para chamar a nova função)
    if len(sys.argv) < 3:
        print("Uso: python <nome_do_script>.py <GENE_SYMBOL> <caminho_de_saida.fasta>")
        sys.exit(1)

    gene_symbol = sys.argv[1]
    out_path = sys.argv[2]

    species_file = ROOT_DIR / "data" / "species" / "mammals.txt"
    species_set = load_species(species_file)

    print(f"--- Iniciando processo para o gene: {gene_symbol} ---")
    print(f"Filtro de espécies ativo para {len(species_set)} espécies.")

    ortholog_ids = get_all_orthologs_ids(gene_symbol, species_filter_set=species_set)

    if ortholog_ids is None or not ortholog_ids:
        print(f"Nenhum ortólogo 1:1 válido encontrado para {gene_symbol}. Encerrando.")
        Path(out_path).write_text("")
        sys.exit(0)
    
    # Chamada da nova função principal
    validated_seqs = fetch_canonical_cds_from_genes(ortholog_ids)
    
    if not validated_seqs:
        print(f"Nenhuma sequência de CDS válida pôde ser baixada para {gene_symbol}.")
        Path(out_path).write_text("")
        sys.exit(0)
                
    os.makedirs(Path(out_path).parent, exist_ok=True)
    save_as_fasta(gene_symbol, validated_seqs, out_path)
    print(f"\nSUCESSO: Salvas {len(validated_seqs)} sequências para {gene_symbol} em {out_path}")






