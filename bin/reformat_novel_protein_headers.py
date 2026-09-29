#!/usr/bin/env python3
import argparse
import re
from pathlib import Path


def parse_fasta(path):
    header = None
    seq_chunks = []

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                if header is None:
                    raise ValueError(f"Sequence found before first FASTA header in {path}")
                seq_chunks.append(line.replace(" ", "").replace("\t", ""))

    if header is not None:
        yield header, "".join(seq_chunks)


def write_fasta_record(handle, header, sequence, wrap=60):
    handle.write(f">{header}\n")
    for start in range(0, len(sequence), wrap):
        handle.write(sequence[start : start + wrap] + "\n")


def safe_token(value):
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    cleaned = cleaned.strip("_.-")
    return cleaned or "UNKNOWN"


def extract_gene_name(original_header, protein_id):
    match = re.search(r"\bGENE\.([^\s~]+)", original_header)
    if match:
        return safe_token(match.group(1))

    # TransDecoder protein IDs look like STRG.1006.1.p19; the transcript-like
    # portion is the best fallback for GN=.
    if ".p" in protein_id:
        return safe_token(protein_id.rsplit(".p", 1)[0])
    return safe_token(protein_id)


def make_uniprot_like_header(original_header, sample_id, tissue, organism, taxonomy_id, reviewed):
    protein_id = original_header.split(None, 1)[0]
    safe_protein = safe_token(protein_id)
    gene_name = extract_gene_name(original_header, protein_id)

    db_prefix = "sp" if reviewed else "st"
    accession = safe_protein
    entry_name = safe_protein

    tissue_field = f" Tissue={safe_token(tissue)}" if tissue else ""

    return (
        f"{db_prefix}|{accession}|{entry_name} "
        f"OS={organism} OX={taxonomy_id} GN={gene_name} PE=5 SV=1 Sample={safe_token(sample_id)}{tissue_field} "
        f"OriginalID={safe_protein} OriginalHeader={original_header}"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Reformat TransDecoder novel protein FASTA headers to UniProt-like headers."
    )
    parser.add_argument("--input", required=True, help="Input novel protein FASTA")
    parser.add_argument("--output", required=True, help="Output FASTA with UniProt-like headers")
    parser.add_argument("--sample-id", required=True, help="Sample/database identifier")
    parser.add_argument("--tissue", default=None, help="Tissue name for Tissue=")
    parser.add_argument("--organism", required=True, help="Organism name for OS=")
    parser.add_argument("--taxonomy-id", required=True, help="NCBI taxonomy id for OX=")
    parser.add_argument(
        "--reviewed",
        action="store_true",
        help="Use sp| instead of st|. Novel entries should normally remain st|.",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    count = 0
    with open(output_path, "w", encoding="utf-8") as out_handle:
        for original_header, sequence in parse_fasta(input_path):
            new_header = make_uniprot_like_header(
                original_header=original_header,
                sample_id=args.sample_id,
                tissue=args.tissue,
                organism=args.organism,
                taxonomy_id=args.taxonomy_id,
                reviewed=args.reviewed,
            )
            write_fasta_record(out_handle, new_header, sequence)
            count += 1

    print(f"Reformatted {count} novel protein headers: {input_path} -> {output_path}")


if __name__ == "__main__":
    main()
