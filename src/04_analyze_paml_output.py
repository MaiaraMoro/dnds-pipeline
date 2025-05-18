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
    current_model = None

    for i, line in enumerate(lines):
        # detecta o nome do modelo nas linhas que precedem o lnL
        if "Model 0" in line:
            current_model = "M0"
        elif "Model 1" in line:
            current_model = "M1a"
        elif "Model 2" in line:
            current_model = "M2a"
        elif "Model 7" in line:
            current_model = "M7"
        elif "Model 8" in line:
            current_model = "M8"

        # captura lnL e np
        if "lnL" in line and "np:" in line:
            match = re.search(r"lnL.*np:\s+(\d+)\)\s*:\s*(-?\d+\.\d+)", line)
            if match and current_model:
                np = int(match.group(1))
                lnL = float(match.group(2))
                models[current_model] = {"lnL": lnL, "np": np}
                current_model = None  # reseta para próxima leitura

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
    """
    Extracts BEB sites from the MLC file.
    """
    with open(filepath, 'r') as file:
        lines = file.readlines()

    results = []
    recording = False
    model_found = False

    for i, line in enumerate(lines):
        if f"Model {model_tag}" in line:
            model_found = True
            
        if model_found and "Bayes Empirical Bayes" in line:
            recording = True
            continue
        
        if recording:
            if re.match(f"^\s*\d+", line):
                parts = line.strip().split()
                if len(parts) >= 5:
                    try:
                        site = int(parts[0])
                        aa = parts[1]
                        prob = float(parts[2])
                        if prob > 0.95:
                            results.append((model_tag, site, aa, prob))
                    except ValueError:
                        continue
                    
            elif line.strip() == "":
                break 

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
        out_beb.write("Gene\tModel\tSite\tAA\tProb\n")
        for model in ["M2a", "M8"]:
            beb_sites = extract_beb_sites(path, model_tag=model)
            for tag, site, aa, prob in beb_sites:
                out_beb.write(f"{gene}\t{tag}\t{site}\t{aa}\t{prob:.3f}\n")


if __name__ == "__main__":
    import sys
    gene = sys.argv[1]
    output_lrt = sys.argv[2]
    output_beb = sys.argv[3]
    main(gene, output_lrt, output_beb)
