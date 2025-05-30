import os
from collections import defaultdict

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
    mean_percents = [sum(sample[j] for sample in sv_percentage_matrix) / len(sv_percentage_matrix) for j in range(num_svs)]

    n_top = max(1, int(num_svs * top_ratio))
    top_indices = sorted(
        sorted(range(num_svs), key=lambda i: mean_percents[i], reverse=True)[:n_top]
    )
    selected_sv_cols = [sv_columns[i] for i in top_indices]
    return selected_sv_cols


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

# 상대 경로 설정
meta_file = os.path.join(CURRENT_DIR, "data", "jsb_preprocessed_metadata_25_05.tsv")
full_meta_file = os.path.join(CURRENT_DIR, "data", "Final_metadata_jsb.tsv")
sv_data_file = os.path.join(CURRENT_DIR, "data", "Allsamples_dungbeetle_df.tsv")
merged_output_file = os.path.join(CURRENT_DIR, "output", "preprocessed_all_data.tsv")

# 샘플 필터링
using_samples = set()
with open(meta_file) as f:
    f.readline()
    for line in f:
        smp_id = line.split("\t")[0]
        using_samples.add(smp_id)

# 메타데이터 수집
using_species = ["gibbulus", "mopsus", "nuchicornis", "sordescens", "laticornis"]
smp2metadata = dict()
with open(full_meta_file) as f:
    header = ["smp_id"] + f.readline().rstrip("\n").split("\t")
    col2idx = {col: idx for idx, col in enumerate(header)}
    for line in f:
        row = line.rstrip("\n").split("\t")
        sample_id = row[col2idx["Sample_id"]]
        if sample_id not in using_samples:
            continue
        host_species = row[col2idx["Species"]]
        if host_species not in using_species:
            continue
        host_genus = row[col2idx["Genus"]]

        using_data = [
            sample_id,
            row[col2idx["Site_type"]],
            f"{host_genus}_{host_species}",
        ]
        using_data.extend(row[col2idx["bio1"]:col2idx["SG"]])
        smp2metadata[sample_id] = using_data

# SV 데이터 수집
smp2svdata = dict()
with open(sv_data_file) as f:
    header = f.readline().rstrip("\n").split("\t")
    col2idx = {col: idx for idx, col in enumerate(header)}
    for line in f:
        row = line.rstrip("\n").split("\t")
        sample_id = row[col2idx["Sample_id"]]
        if sample_id not in using_samples:
            continue
        smp2svdata[sample_id] = row[col2idx["SV1"]:]

# 전체 데이터 저장
meta_header = [
    "Sample_id", "Site_type", "Host_Genus_Species", "bio1", "bio2", "bio3", "bio4", "bio5", "bio6",
    "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17",
    "bio18", "bio19"
]
sv_header = header[1:]
all_header = meta_header + sv_header

os.makedirs(os.path.dirname(merged_output_file), exist_ok=True)
with open(merged_output_file, "w") as f:
    f.write("\t".join(all_header) + "\n")
    for sample in smp2metadata:
        if sample not in smp2svdata:
            continue
        all_data = smp2metadata[sample] + smp2svdata[sample]
        f.write("\t".join(all_data) + "\n")

# 종별로 분리
species2all_data = defaultdict(list)
with open(merged_output_file) as f:
    header = f.readline()
    col2idx = {col: idx for idx, col in enumerate(header.rstrip("\n").split("\t"))}
    for line in f:
        species = line.split("\t")[col2idx["Host_Genus_Species"]]
        species2all_data[species].append(line)

for species, lines in species2all_data.items():
    species_file = os.path.join(CURRENT_DIR, "output", f"{species}_preprocessed_all_data.tsv")
    with open(species_file, "w") as f:
        f.write(header)
        for line in lines:
            f.write(line)

# SV 상위 필터링
for species in species2all_data:
    input_file = os.path.join(CURRENT_DIR, "output", f"{species}_preprocessed_all_data.tsv")
    top_sv_col = get_top_sv_column_names(input_file)

    new_rows = []
    with open(input_file) as f:
        header = f.readline().rstrip("\n").split("\t")
        col2idx = {col: idx for idx, col in enumerate(header)}
        new_header = ["Sample_id", "Site_type"] + top_sv_col
        new_rows.append(new_header)
        for line in f:
            row = line.rstrip("\n").split("\t")
            new_line = [row[col2idx["Sample_id"]], row[col2idx["Site_type"]]]
            new_line += [row[col2idx[sv]] for sv in top_sv_col]
            new_rows.append(new_line)

    output_file = os.path.join(CURRENT_DIR, "output", f"{species}_sv_filtered_all_data.tsv")
    with open(output_file, "w") as g:
        for line in new_rows:
            g.write("\t".join(line) + "\n")
