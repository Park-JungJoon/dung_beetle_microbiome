import pandas as pd
import os
import json
import numpy as np
from skbio.diversity import beta_diversity
from skbio.stats.distance import permdisp, mantel
from sklearn.preprocessing import StandardScaler
from skbio.stats.distance import DistanceMatrix

def calculate_permisp(file_path, host_species_name):
    df = pd.read_csv(file_path, sep="\t")

    # Start SVs from index 2 (3rd column)
    sv_start_idx = 2
    sv_cols = df.columns[sv_start_idx:]
    sv_counts = df[sv_cols]
    sv_rel_abundance = sv_counts.div(sv_counts.sum(axis=1), axis=0)

    # Metadata
    sample_ids = df["Sample_id"]
    site_type = df["Site_type"]
    site_type_series = pd.Series(site_type.values, index=sample_ids)

    # Bray-Curtis distance matrix
    bray_dm = beta_diversity("braycurtis", sv_rel_abundance, ids=sample_ids)

    # PERMDISP
    site_type_series = pd.Series(site_type.values, index=sample_ids, name="Site_type")
    permdisp_result = permdisp(bray_dm, site_type_series)

    # Mantel test: binary distance matrix for site type
    site_matrix = np.array([
        [0 if a == b else 1 for b in site_type]
        for a in site_type
    ])
    site_dm = DistanceMatrix(site_matrix, ids=sample_ids)

    # Mantel test
    mantel_result = mantel(bray_dm, site_dm, method='pearson', permutations=999)

    # Save results
    results = {
        "PERMDISP_F_value": permdisp_result["test statistic"],
        "PERMDISP_p_value": permdisp_result["p-value"],
        "PERMDISP_sample_size": permdisp_result["sample size"],
        "PERMDISP_n_groups": permdisp_result["number of groups"],
        "PERMDISP_n_permutations": permdisp_result["number of permutations"],
        "Mantel_r": mantel_result[0],
        "Mantel_p_value": mantel_result[1]
    }

    output_path = os.path.join(os.path.dirname(file_path), f"{host_species_name}_permdisp_mantel.json")
    with open(output_path, "w") as f:
        json.dump(results, f, indent=4)

# 기준 경로 설정
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(CURRENT_DIR, "output")

# Run for each filtered file
files = os.listdir(OUTPUT_DIR)
files = [file for file in files if "sv_filtered" in file]

for file in files:
    sample_id = file.split("_sv_filtered")[0]
    file_path = os.path.join(OUTPUT_DIR, file)
    calculate_permisp(file_path, sample_id)
