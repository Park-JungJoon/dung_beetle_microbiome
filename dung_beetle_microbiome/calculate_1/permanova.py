import os
import json
import pandas as pd
from skbio.diversity import beta_diversity
from skbio.stats.distance import permanova
from skbio import DistanceMatrix

def get_top_sv_column_names(file_path: str, top_ratio: float = 0.1) -> list:
    """
    TSV 파일에서 각 샘플의 SV 비율을 계산하고, 평균 비율 기준 상위 top_ratio에 해당하는 SV 컬럼 이름만 반환.
    """
    with open(file_path, "r") as f:
        header = f.readline().strip().split("\t")
        sv_start_idx = header.index("SV1")
        sv_columns = header[sv_start_idx:]

        sv_percentage_matrix = []
        for line in f:
            row = line.strip().split("\t")
            sv_values = list(map(float, row[sv_start_idx:]))
            total = sum(sv_values)
            if total == 0:
                continue
            sv_percent = [val / total for val in sv_values]
            sv_percentage_matrix.append(sv_percent)

    num_svs = len(sv_columns)
    mean_percents = [
        sum(sample[j] for sample in sv_percentage_matrix) / len(sv_percentage_matrix)
        for j in range(num_svs)
    ]

    n_top = max(1, int(num_svs * top_ratio))
    top_indices = sorted(
        sorted(range(num_svs), key=lambda i: mean_percents[i], reverse=True)[:n_top]
    )

    selected_sv_cols = [sv_columns[i] for i in top_indices]
    return selected_sv_cols


def calculate_permanova(species, input_dir, output_dir):
    input_path = os.path.join(input_dir, f"{species}_preprocessed_all_data.tsv")
    df = pd.read_csv(input_path, sep="\t")

    top_sv_cols = get_top_sv_column_names(input_path)
    sv_df = df[top_sv_cols]
    filtered_df = sv_df[top_sv_cols]

    sample_ids = df["Sample_id"].values
    site_type = df["Site_type"].values

    sv_rel = filtered_df.div(filtered_df.sum(axis=1), axis=0)
    bray_curtis_dm = beta_diversity("braycurtis", sv_rel, ids=sample_ids)

    result = permanova(distance_matrix=bray_curtis_dm, grouping=site_type, permutations=999)
    
    result_dict = {
        "test_statistic": result["test statistic"],
        "p_value": result["p-value"],
        "sample_size": result["sample size"],
        "num_groups": result["number of groups"],
        "stat_name": result["test statistic name"],
        "method": result["method name"],
        "num_permutations": result["number of permutations"]
    }

    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"{species}_permanova.json")
    with open(output_path, "w") as f:
        json.dump(result_dict, f, indent=4)


if __name__ == "__main__":
    CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
    INPUT_DIR = os.path.join(CURRENT_DIR, "output")
    RESULT_DIR = os.path.join(CURRENT_DIR, "calculate_1", "output")

    species_paths = [
        fname for fname in os.listdir(INPUT_DIR)
        if fname.endswith("_preprocessed_all_data.tsv")
    ]
    species_list = [fname.split("_preprocessed_all_data.tsv")[0] for fname in species_paths]

    for species in species_list:
        calculate_permanova(species, INPUT_DIR, RESULT_DIR)
