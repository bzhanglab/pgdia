import argparse
import os
from pathlib import Path
import re

import pandas as pd
from Bio import SeqIO


sample_name = ""
data_path_prefix = ""
data_path = ""
isoform_annotation_path = None
novel_fasta_path = None
reference_fasta_path = None

seq_dict = {}
ensp_ensg_dict = {}
enst_ensg_dict = {}
ensg_gene_dict = {}
known_seqs = []

REFERENCE_PROTEIN_PREFIXES = ("ENSP", "ENSMUSP")
REFERENCE_TRANSCRIPT_RE = re.compile(r"(ENST\d+(?:\.\d+)?|ENSMUST\d+(?:\.\d+)?)")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Postprocess DIA-NN report parquet with global/run-specific FDR filtering."
    )
    parser.add_argument(
        "--sample-name",
        help="Optional sample name label. Defaults to value inferred from --report-parquet filename.",
    )
    parser.add_argument(
        "--report-parquet",
        required=True,
        help="Path to DIA-NN report parquet file.",
    )
    parser.add_argument(
        "--output-prefix",
        help="Output file prefix. Defaults to report basename without _report.parquet.",
    )
    parser.add_argument(
        "--novel-fasta",
        help="Novel/combined protein FASTA used for variant annotation.",
    )
    parser.add_argument(
        "--isoform-annotation",
        help="StringTie .tmap file for novel isoform gene mapping.",
    )
    parser.add_argument(
        "--reference-fasta",
        help="Reference protein FASTA used to map ENSP/ENST to gene IDs/symbols.",
    )
    return parser.parse_args()


def configure_context(args):
    global sample_name
    global data_path_prefix
    global data_path
    global isoform_annotation_path
    global novel_fasta_path
    global reference_fasta_path

    report_name = Path(args.report_parquet).name
    if report_name.endswith("_report.parquet"):
        inferred_prefix = report_name[: -len("_report.parquet")]
    elif report_name.endswith(".parquet"):
        inferred_prefix = report_name[: -len(".parquet")]
    else:
        inferred_prefix = report_name

    sample_name = args.sample_name or inferred_prefix
    data_path_prefix = args.output_prefix or inferred_prefix
    data_path = args.report_parquet
    isoform_annotation_path = args.isoform_annotation
    novel_fasta_path = args.novel_fasta
    reference_fasta_path = args.reference_fasta


def load_reference_maps(reference_fasta):
    seq_dict.clear()
    ensp_ensg_dict.clear()
    enst_ensg_dict.clear()
    ensg_gene_dict.clear()
    known_seqs.clear()

    if not reference_fasta:
        print("Reference FASTA was not provided; gene symbol annotation will be limited.")
        return

    if not os.path.exists(reference_fasta):
        print(f"Reference FASTA not found: {reference_fasta}. Gene symbol annotation will be limited.")
        return

    for record in SeqIO.parse(reference_fasta, "fasta"):
        known_seqs.append(str(record.seq))
        if record.id.startswith("Cont|"):
            continue

        parts = record.id.split("|")
        if len(parts) < 3:
            continue

        rid = parts[1]
        seq_dict[rid] = record

        rids = rid.split("_")
        if len(rids) >= 3:
            enst_ensg_dict[rids[1]] = rids[2]
            ensp_ensg_dict[rids[0]] = rids[2]
            enst_ensg_dict[rids[1].split(".")[0]] = rids[2]
            ensp_ensg_dict[rids[0].split(".")[0]] = rids[2]
            gene_symbol = parts[2].split()[0]
            ensg_gene_dict[rids[2]] = gene_symbol
            ensg_gene_dict[rids[2].split(".")[0]] = gene_symbol


def convert_parquet_to_tsv(file_path):
    df = pd.read_parquet(file_path, engine="pyarrow")

    return df


def reference_record_id(entry):
    token = str(entry).split(";")[0]
    parts = token.split("|")
    if len(parts) >= 3:
        return parts[1]
    return token


def is_reference_protein_id(entry):
    rid = reference_record_id(entry)
    return rid.startswith(REFERENCE_PROTEIN_PREFIXES)


def reference_gene_id(entry):
    rid = reference_record_id(entry)
    fields = rid.split("_")
    if len(fields) >= 3:
        return fields[2]
    return ensp_ensg_dict.get(rid, ensp_ensg_dict.get(rid.split(".")[0], "None"))


def prioritize_ensp(entry):
    if ";" in entry:
        parts = entry.split(";")
        first_part = parts[0]
        if not is_reference_protein_id(first_part):
            ensp_candidates = [p for p in parts if is_reference_protein_id(p)]
            return ensp_candidates[0] if ensp_candidates else first_part
        return first_part
    return entry


def get_id_num(df):
    # check if the dataframe has decoys and Cont|
    # df = df[(df["Decoy"] == 0) & (~df['Protein.Ids'].str.startswith("Cont|"))]
    num_peptides = df["Stripped.Sequence"].nunique()
    print(f"The number of peptides is {num_peptides}")


def filter_at_fdr_cutoff(df_subset, df_all, sort_column, is_novel):
    if df_subset.empty:
        return df_subset, df_subset, df_subset
    # Sort by Score Q.Value column(lower is better)
    df_subset = df_subset.sort_values(by=sort_column, ascending=True).reset_index(drop=True)

    # Compute cumulative counts
    df_subset["Cumulative Decoys"] = df_subset["Decoy"].cumsum()
    df_subset["Cumulative Targets"] = (1 - df_subset["Decoy"]).cumsum()

    df_subset["D+"] = df_subset["Cumulative Decoys"]
    df_subset["T+_V"] = df_subset["Cumulative Targets"]

    D_V = df_all[(df_all["Decoy"] == 1) & (df_all["Novel"] == is_novel)].shape[0]  # Total subtype decoy PSMs
    D = df_all[df_all["Decoy"] == 1].shape[0]  # Total decoy PSMs

    if sort_column == "Q.Value":
        fdr_column = "Cumulative FDR"
    else:
        fdr_column = "Cumulative global FDR"
    df_subset[fdr_column] = (df_subset["D+"] / df_subset["T+_V"]) * (D_V / D) if D > 0 else 0

    # Filter at 1% FDR
    df_filtered = df_subset[df_subset["Cumulative FDR"] <= 0.01]
    if "Cumulative global FDR" in df_filtered.columns:
        df_filtered = df_filtered[df_filtered["Cumulative global FDR"] <= 0.01]
    # print(df_filtered.shape[0])

    df_filtered_no_decoy = df_filtered[(df_filtered["Decoy"] == 0) & (~df_filtered['Protein.Ids'].str.startswith("Cont|"))]
    df_filtered_no_decoy = df_filtered_no_decoy.reset_index(drop=True)

    num_peptides = df_filtered_no_decoy["Stripped.Sequence"].nunique()
    print(f"The number of peptides is {num_peptides}")

    return df_subset, df_filtered, df_filtered_no_decoy


def compute_global_fdr(df_subset, sort_column):
    if df_subset.empty:
        return df_subset, df_subset, df_subset
    df_subset = df_subset.sort_values(by=sort_column, ascending=True).reset_index(drop=True)

    df_subset["Cumulative Decoys"] = df_subset["Decoy"].cumsum()
    df_subset["Cumulative Targets"] = (1 - df_subset["Decoy"]).cumsum()

    # fdr_column = " "
    if sort_column == "Q.Value":
        fdr_column = "Cumulative FDR"
    else:
        fdr_column = "Cumulative global FDR"

    df_subset[fdr_column] = df_subset["Cumulative Decoys"] / df_subset["Cumulative Targets"]

    df_filtered = df_subset[df_subset[fdr_column] <= 0.01]

    df_filtered_no_decoy = df_filtered[(df_filtered["Decoy"] == 0) & (~df_filtered['Protein.Ids'].str.startswith("Cont|"))]
    df_filtered_no_decoy = df_filtered_no_decoy.reset_index(drop=True)

    return df_subset, df_filtered, df_filtered_no_decoy


def compute_separate_fdr(df_all, fdr_level):
    # Apply the function to modify the "Protein.Ids" column
    # if there are multiple protein assignment for one peptide sequence, prioritize the reference one
    df_all["Processed_Protein.Ids"] = df_all["Protein.Ids"].apply(prioritize_ensp)
    global_score_cutoff = 0.01

    # df_all = df_all.sort_values(by="Q.Value", ascending=True).reset_index(drop=True)
    # Split into two groups, separate all IDs in the report
    df_other = df_all[~df_all["Processed_Protein.Ids"].apply(is_reference_protein_id)]  # Non-ref group

    if fdr_level == "psm":
        global_score_column = "Q.Value"
        df_filtered = df_all[df_all[global_score_column] < global_score_cutoff]

        D = df_all["Decoy"].sum()
        D_plus = df_filtered["Decoy"].sum()

        D_v = df_other["Decoy"].sum()
        T_plus_v = df_other[(~df_other["Decoy"]) & (df_other[global_score_column] < global_score_cutoff)].shape[0]
        fdr_v = (D_plus / T_plus_v) * (D_v / D) if T_plus_v > 0 and D > 0 else 0

        return fdr_v
    else:
        return -1


def split_by_run_index(df):
    """
    Splits a DataFrame into separate DataFrames based on the 'Run.index' column.

    Parameters:
    df (pd.DataFrame): Input DataFrame containing a 'Run.index' column.

    Returns:
    dict: A dictionary where keys are unique 'Run.index' values and values are the corresponding DataFrames.
    """
    return {r: df[df["Run.Index"] == r].reset_index(drop=True) for r in df["Run.Index"].unique()}


def save_unique_peptides(df, output_file_path, novel_flag=False):
    # Filter target peptides (Decoy == 0) and remove contaminants (Protein.Ids does not start with "Cont|")
    n_isoform_peptides = 0
    n_variant_peptides = 0
    if "Decoy" in df.columns:
        df = df[(df["Decoy"] == 0) & (~df['Protein.Ids'].str.startswith("Cont|", na=False))]

    else:
        df = df[~df['Protein.Ids'].str.startswith("Cont|", na=False)]
    unique_peptides = df["Stripped.Sequence"].unique()
    unique_precursors = df["Precursor.Id"].unique()
    unique_proteins = df["Genes"].unique()

    print(f"Precursors: {len(unique_precursors)}")
    print(f"Peptides: {len(unique_peptides)}")
    print(f"Proteins: {len(unique_proteins)}")

    if novel_flag:
        n_isoform_peptides = df.loc[df['Processed_Protein.Ids'].str.startswith("STRG"), "Stripped.Sequence"].nunique()
        n_variant_peptides = len(unique_peptides) - n_isoform_peptides
        print(f"Variant Peptides: {n_variant_peptides}")
        print(f"Isoform Peptides: {n_isoform_peptides}")

    # Save to a TXT file (one peptide per line)
    with open(output_file_path, "w") as f:
        for seq in unique_peptides:
            f.write(seq + "\n")

    print(f"Unique peptides saved to {output_file_path}")


def get_pr_matrix(df):
    sub_dfs = split_by_run_index(df)
    intensity_dfs = []
    for ri, sub_df in sub_dfs.items():
        # Rename intensity column using Run.index
        intensity_column_name = f"Run_{ri}_Intensity"
        sub_df_intensity = sub_df[["Precursor.Id", "Precursor.Normalised"]].rename(
            columns={"Precursor.Normalised": intensity_column_name})
        intensity_dfs.append(sub_df_intensity)
    intensity_df = intensity_dfs[0]
    for df_int in intensity_dfs[1:]:
        intensity_df = intensity_df.merge(df_int, on="Precursor.Id", how="outer")
    sorted_columns = sorted(intensity_df.columns[1:], key=lambda x: int(x.split('_')[1]))
    updated_columns = [intensity_df.columns[0]]
    updated_columns.extend(sorted_columns)
    sorted_intensity_df = intensity_df[updated_columns]

    columns_to_keep = [
        "Precursor.Id", "Protein.Group", "Protein.Ids", "Processed_Protein.Ids", "Protein.Names", "Genes",
        "Proteotypic", "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge", "Novel"
    ]  # Core metadata columns

    final_combined_df = (pd.concat(list(sub_dfs.values()), ignore_index=True)[columns_to_keep].
                         drop_duplicates(subset=["Precursor.Id"]))
    final_combined_df_int = final_combined_df.merge(sorted_intensity_df, on="Precursor.Id", how="left")

    return final_combined_df_int


def check_novel_seqs(novel_df):

    # reference_seqs = get_known_seqs("../rna-seq-results/combined_protein_db/GENCODE.V42.basic.CHR_uniprot_like.fa")
    not_novel_cnt = 0
    for i in range(novel_df.shape[0]):
        pep = novel_df.at[i, "Stripped.Sequence"]
        if pep in known_seqs:
            print(novel_df.at[i, "Protein.Ids"], i)
            not_novel_cnt += 1
    print(f"There are {not_novel_cnt} of peptides which are not novel")


def extract_var_annotation(fasta_path):
    if not fasta_path or not os.path.exists(fasta_path):
        return {}

    fasta_dict = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        header = rec.id
        if header.startswith("var_") and REFERENCE_TRANSCRIPT_RE.search(header):
            k = header.split(";")[0]
            m_transcript = re.search(r"_(ENST\d+|ENSMUST\d+)", k)
            if m_transcript:
                k = k[:m_transcript.start()]
            fasta_dict[k] = header

    return fasta_dict


def annotate_var_id(protein_id, annotate_dict):
    m_transcript = re.search(r"_(ENST\d+|ENSMUST\d+)", protein_id)
    if m_transcript:
        protein_id = protein_id[:m_transcript.start()]
    if protein_id in annotate_dict:
        return annotate_dict[protein_id]
    return protein_id


def extract_novel_isoform_annotation(file_path):
    if not file_path or not os.path.exists(file_path):
        return {}

    # file1_path = f'../rna-seq-results/stringtie/output/{sample_name}/t.{sample_name}.stringtie_output.gtf.tmap'
    df1 = pd.read_csv(file_path, sep="\t", header=0)  # StringTie annotation
    df1 = df1[df1["class_code"].isin(["j", "i", "m", "x", "u", "o"])]  # j|i|m|x|u|o
    df1 = df1[df1["TPM"] > 8]

    gene_strg_id_dict = df1.set_index('qry_id')['ref_gene_id'].to_dict()

    return gene_strg_id_dict


def symbol_from_gene_id(gid: str) -> str:
    if gid and gid not in ("None", "NotPrimary", "NA"):
        return ensg_gene_dict.get(gid, ensg_gene_dict.get(gid.split(".")[0], "NA"))
    return "NA"


def get_gene_info(input_df, var_fasta_path, isoform_file_path):
    var_fasta_dict = extract_var_annotation(var_fasta_path)
    isoform_gene_dict = extract_novel_isoform_annotation(isoform_file_path)
    # start_pos = []
    # end_pos = []
    genes = []
    gene_symbols = []

    for i in range(input_df.shape[0]):
        protein_id = input_df.at[i, "Processed_Protein.Ids"]
        gene_id = "None"
        gene_symbol = "NA"

        if is_reference_protein_id(protein_id):
            gene_id = reference_gene_id(protein_id)
            gene_symbol = symbol_from_gene_id(gene_id)
        elif "STRG" in protein_id:
            protein_id = protein_id.split(";")[0]
            gene_flag = protein_id.find("p")
            qry_gene_id = protein_id[:gene_flag-1]
            gene_id = isoform_gene_dict.get(qry_gene_id, qry_gene_id)
            gene_symbol = symbol_from_gene_id(gene_id)

        elif "var" in protein_id:
            new_protein_id = annotate_var_id(protein_id, var_fasta_dict)
            input_df.at[i, "Processed_Protein.Ids"] = new_protein_id
            # 1) Get gene_id (ENSG) via ENST -> ENSG if possible
            m_enst = REFERENCE_TRANSCRIPT_RE.search(new_protein_id)
            if m_enst:
                enst_id = m_enst.group(1)
                gene_id = enst_ensg_dict.get(enst_id, "NotPrimary")

            # 2) Prefer gene symbol from gene_id_dict (ensg_gene_dict)
            gene_symbol = symbol_from_gene_id(gene_id)

            # 3) Fallback: extract gene symbol from the string (e.g., _ENST..._MRPL15_)
            if gene_symbol == "NA":
                m_sym = re.search(r'(?:ENST\d+(?:\.\d+)?|ENSMUST\d+(?:\.\d+)?)_([A-Za-z0-9]+)', new_protein_id)
                if m_sym:
                    gene_symbol = m_sym.group(2)

        else:
            print("Not reference ID")
            gene_symbol = "NA"

        genes.append(gene_id)
        gene_symbols.append(gene_symbol)

    input_df["Genes"] = genes
    input_df["Gene_Symbols"] = gene_symbols

    return input_df



def process_report_with_separate_fdr():
    df_fdr = convert_parquet_to_tsv(data_path)

    df_fdr["Processed_Protein.Ids"] = df_fdr["Protein.Ids"].apply(prioritize_ensp)
    df_fdr["Novel"] = ~df_fdr["Processed_Protein.Ids"].apply(is_reference_protein_id)

    run_dfs = split_by_run_index(df_fdr)
    filtered_ensp = []
    filtered_novel = []
    ensp_sub_dfs_fdr = []  # run-specific ensp dfs with new FDR column
    novel_sub_dfs_fdr = []  # run-specific novel dfs with new FDR column

    for run_index, run in run_dfs.items():
        # PSM level FDR calculation
        df_ensp = run[~run["Novel"]]
        df_novel = run[run["Novel"]]

        df_sub_ensp, df_ensp_filtered, _ = filter_at_fdr_cutoff(df_ensp, run, "Q.Value", False)
        df_sub_novel, df_novel_filtered, _ = filter_at_fdr_cutoff(df_novel, run, "Q.Value", True)

        filtered_ensp.append(df_ensp_filtered)
        filtered_novel.append(df_novel_filtered)

        ensp_sub_dfs_fdr.append(df_sub_ensp)
        novel_sub_dfs_fdr.append(df_sub_novel)

    merged_ensp_df = pd.concat(ensp_sub_dfs_fdr, ignore_index=True)
    merged_novel_df = pd.concat(novel_sub_dfs_fdr, ignore_index=True)
    print(merged_ensp_df.shape[0], merged_novel_df.shape[0])

    _, _, final_ensp_df = filter_at_fdr_cutoff(merged_ensp_df, df_fdr, "Global.Q.Value", False)
    _, _, final_novel_df = filter_at_fdr_cutoff(merged_novel_df, df_fdr, "Global.Q.Value", True)

    final_ensp_df = final_ensp_df[(final_ensp_df["Cumulative FDR"] <= 0.01)]
    final_novel_df = final_novel_df[(final_novel_df["Cumulative FDR"] <= 0.01)]

    final_combined_ensp_int = get_pr_matrix(final_ensp_df)
    final_combined_novel_int = get_pr_matrix(final_novel_df)

    final_combined_ensp_int = final_combined_ensp_int.sort_values(by="Precursor.Id").reset_index(drop=True)
    final_combined_novel_int = final_combined_novel_int.sort_values(by="Precursor.Id").reset_index(drop=True)

    save_unique_peptides(final_combined_ensp_int, data_path_prefix + "_separate_reference_fdr_peptide.txt")
    save_unique_peptides(final_combined_novel_int, data_path_prefix + "_separate_novel_fdr_peptide.txt")

    get_id_num(final_combined_ensp_int)
    get_id_num(final_combined_novel_int)
    final_combined_ensp_int.to_csv(data_path_prefix + "_separate_reference_fdr.tsv", sep='\t', index=None)
    final_combined_novel_int.to_csv(data_path_prefix + "_separate_novel_fdr.tsv", sep='\t', index=None)


def recalculate_global_fdr():
    df_fdr = convert_parquet_to_tsv(data_path)

    df_fdr["Processed_Protein.Ids"] = df_fdr["Protein.Ids"].apply(prioritize_ensp)
    df_fdr["Novel"] = ~df_fdr["Processed_Protein.Ids"].apply(is_reference_protein_id)

    run_dfs = split_by_run_index(df_fdr)

    # Compute run-specific and global FDR together
    run_specific_filtered_dfs = []
    fdr_run_dfs = []
    for run_index, run_df in run_dfs.items():
        fdr_run_df, filtered_run_df, _ = compute_global_fdr(run_df, "Q.Value")
        fdr_run_dfs.append(fdr_run_df)
        run_specific_filtered_dfs.append(filtered_run_df)

    merged_fdr_run_df = pd.concat(fdr_run_dfs, ignore_index=True)
    # merged_run_df = pd.concat(run_specific_filtered_dfs, ignore_index=True)
    global_df, _, global_fdr_df = compute_global_fdr(merged_fdr_run_df, "Global.Q.Value")
    print(f"global FDR only: {global_fdr_df.shape[0]}")
    get_id_num(global_fdr_df)
    global_fdr_df = global_fdr_df[(global_fdr_df["Cumulative FDR"] <= 0.01)]
    print(f"global FDR and run-specific FDR: {global_fdr_df.shape[0]}")
    get_id_num(global_fdr_df)
    # global_fdr_df.to_csv(data_path_prefix + "_global_fdr_report.tsv", sep='\t', index=None)

    global_pr_matrix = get_pr_matrix(global_fdr_df)
    global_pr_matrix = global_pr_matrix.sort_values(by="Precursor.Id").reset_index(drop=True)

    if novel_fasta_path or isoform_annotation_path:
        global_pr_matrix = get_gene_info(global_pr_matrix, novel_fasta_path, isoform_annotation_path)
    else:
        print("Skipping gene annotation: --novel-fasta/--isoform-annotation were not provided.")

    global_pr_matrix_ensp = global_pr_matrix[~global_pr_matrix["Novel"]].reset_index(drop=True)
    global_pr_matrix_ensp.to_csv(data_path_prefix + "_ref_matrix.tsv", sep='\t', index=None)

    global_pr_matrix_novel = global_pr_matrix[global_pr_matrix["Novel"]].reset_index(drop=True)
    if known_seqs:
        check_novel_seqs(global_pr_matrix_novel)  # double check whether the novel sequences are included in reference DB
    else:
        print("Skipping novel sequence check: reference FASTA was not loaded.")
    global_pr_matrix_novel.to_csv(data_path_prefix + "_novel_matrix.tsv", sep='\t', index=None)


def main():
    args = parse_args()
    configure_context(args)
    load_reference_maps(reference_fasta_path)

    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Report parquet not found: {data_path}")

    recalculate_global_fdr()


if __name__ == "__main__":
    main()
