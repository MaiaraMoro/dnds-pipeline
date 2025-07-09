import re
import os
import sys
from scipy.stats import chi2

def parse_mlc_file(filepath):
    """
    Parses the MLC file to extract the log-likelihood values.
    """
    with open(filepath, 'r') as file:
        lines = file.readlines()

    models = {}
    current_model_tag = None

    model_map = {
        "Model 0": "M0",
        "Model 1": "M1a",
        "Model 2": "M2a",
        "Model 7": "M7",
        "Model 8": "M8"
    }

    for line in lines:
        for paml_tag, tag in model_map.items():
            if paml_tag in line and "NSsites" in line:
                current_model_tag = tag
                break

        if "lnL" in line and "np:" in line and current_model_tag:
            match = re.search(r"lnL.*np:\s*(\d+)\)\s*:\s*(-?\d+\.\d+)", line)
            if match:
                np = int(match.group(1))
                lnL = float(match.group(2))
                # Só adiciona se o modelo ainda não foi registrado (evita pegar o lnL de M0 para M1a, etc.)
                if current_model_tag not in models:
                    models[current_model_tag] = {"lnL": lnL, "np": np}
                current_model_tag = None  # Reseta para a próxima seção

    return models

def likelihood_ratio_test(lnL0, np0, lnL1, np1):
    """
    Performs the likelihood ratio test.
    """
    stat = 2 * (lnL1 - lnL0)
    df = np1 - np0
    p_value = 1 - chi2.cdf(stat, df)

    return stat, df, p_value

def extract_beb_sites(filepath, model_tag="M8"):
    with open(filepath, "r") as fh:
        content = fh.read()

    # Isola a seção do modelo de interesse (M2a ou M8)
    model_search_tag = "Model 2: Positive" if "M2a" in model_tag else "Model 8: beta&w>1"
    
    # Divide o arquivo pela tag do modelo para isolar a seção correta
    try:
        model_section = content.split(model_search_tag)[1]
    except IndexError:
        # Se o modelo não foi encontrado no arquivo, retorna lista vazia
        return []

    # Dentro da seção do modelo, procura pela tabela BEB
    try:
        beb_section = model_section.split("Bayes Empirical Bayes (BEB) analysis")[1]
    except IndexError:
        return []

    results = []
    # Itera sobre as linhas da seção BEB
    for line in beb_section.split('\n'):
        line = line.strip()
        
        # Procura por uma linha que comece com um número (a linha de dados)
        if re.match(r"^\d+", line):
            parts = line.split()
            try:
                site = int(parts[0])
                aa = parts[1]
                prob_token = parts[2]
                
                prob = float(prob_token.rstrip('*'))
                
                signif = ""
                if prob_token.endswith("**"):
                    signif = "**"
                elif prob_token.endswith("*"):
                    signif = "*"
                
                results.append((model_tag, site, aa, prob, signif))
            except (ValueError, IndexError):
                # Ignora linhas que não podem ser parseadas
                continue
    return results


def main(gene, output_lrt_path, output_beb_path):
    path = f"results/paml/{gene}/mlc"
    if not os.path.exists(path):
        print(f"{gene}\t NO MLC FILE")
        return
    
    models = parse_mlc_file(path)

    with open(output_lrt_path, 'w') as out_lrt:
        out_lrt.write("Gene\tTest\tLRT_statistic\tdf\tpval\n")
        for name, m0, m1 in [
            ("M1a_vs_M2a", "M1a", "M2a"),
            ("M7_vs_M8",  "M7",  "M8"),
        ]:
            if m0 in models and m1 in models:
                stat, df, pval = likelihood_ratio_test(
                    models[m0]["lnL"], models[m0]["np"],
                    models[m1]["lnL"], models[m1]["np"]
                )
                out_lrt.write(f"{gene}\t{name}\t{stat:.3f}\t{df}\t{pval:.4f}\n")
            else:
                out_lrt.write(f"{gene}\t{name}\tNA\tNA\tNA\n")

    with open(output_beb_path, 'w') as out_beb:
        out_beb.write("Gene\tModel\tSite\tAA\tProb\tSignif\n")
        for model in ["M2a", "M8"]:
            beb_sites = extract_beb_sites(path, model_tag=model)
            for tag, site, aa, prob, signif in beb_sites:
                out_beb.write(f"{gene}\t{tag}\t{site}\t{aa}\t{prob:.3f}\t{signif}\n")


if __name__ == "__main__":
    import sys
    gene = sys.argv[1]
    output_lrt = sys.argv[2]
    output_beb = sys.argv[3]
    main(gene, output_lrt, output_beb)
