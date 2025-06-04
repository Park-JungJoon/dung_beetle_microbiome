import pandas as pd
import os

def preprocess_to_cli(df: pd.DataFrame, file_name: str) -> None:
    sample_ids = df["SampleID"]
    class_labels = df["Class"]
    sv_cols = df.columns[2:]

    sv_df = df[sv_cols]
    sv_df.index = sample_ids
    sv_df = sv_df.T

    # 🔄 Class를 가장 위로 (먼저 삽입)
    class_row = pd.DataFrame([class_labels.values], index=["Class"], columns=sv_df.columns)
    sv_df = pd.concat([class_row, sv_df], axis=0)

    sv_df.to_csv(file_name, sep="\t")


files = os.listdir("calculate_3/input")
for file in files:
    sample_id = file.split("_sv_filtered_all_data.tsv")[0]
    df = pd.read_csv(f"calculate_3/input/{file}", sep="\t")
    preprocess_to_cli(df, f"calculate_3/output/{sample_id}_cli_formatted.tsv")
