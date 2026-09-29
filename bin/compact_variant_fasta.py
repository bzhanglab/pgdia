#!/usr/bin/env python3
"""Create a compact variant FASTA and a TSV holding the full provenance."""

import argparse
import csv
import re


def fasta_records(path):
    header, sequence = None, []
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header, sequence = line[1:], []
            else:
                sequence.append(line)
    if header is not None:
        yield header, "".join(sequence)


def field(header, name, default=""):
    match = re.search(rf"(?:^|\s){re.escape(name)}=([^\s]+)", header)
    return match.group(1) if match else default


def write_fasta_record(handle, header, sequence, width=60):
    handle.write(f">{header}\n")
    for start in range(0, len(sequence), width):
        handle.write(sequence[start : start + width] + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_fasta")
    parser.add_argument("output_fasta")
    parser.add_argument("provenance_tsv")
    args = parser.parse_args()

    seen_ids = set()
    with open(args.output_fasta, "w", encoding="utf-8") as fasta, open(
        args.provenance_tsv, "w", encoding="utf-8", newline=""
    ) as tsv:
        writer = csv.DictWriter(
            tsv,
            delimiter="\t",
            fieldnames=[
                "protein_id",
                "gene",
                "amino_acid_change",
                "protein_position",
                "observed_run_count",
                "observed_tissue_count",
                "observed_runs",
                "observed_tissues",
            ],
        )
        writer.writeheader()

        for header, sequence in fasta_records(args.input_fasta):
            protein_id = header.split(None, 1)[0]
            if protein_id in seen_ids:
                raise ValueError(f"Duplicate protein ID: {protein_id}")
            seen_ids.add(protein_id)

            gene = field(header, "GN", protein_id.rsplit("|", 1)[-1])
            change = field(header, "AminoAcidChange")
            position = field(header, "ProteinPosition")
            runs = field(header, "ObservedInRuns")
            tissues = field(header, "ObservedInTissues")
            run_count = len(runs.split(",")) if runs else 0
            tissue_count = len(tissues.split(",")) if tissues else 0

            compact_header = (
                f"{protein_id} GN={gene} AA={change} Pos={position} "
                f"Runs={run_count} Tissues={tissue_count}"
            )
            write_fasta_record(fasta, compact_header, sequence)
            writer.writerow(
                {
                    "protein_id": protein_id,
                    "gene": gene,
                    "amino_acid_change": change,
                    "protein_position": position,
                    "observed_run_count": run_count,
                    "observed_tissue_count": tissue_count,
                    "observed_runs": runs,
                    "observed_tissues": tissues,
                }
            )


if __name__ == "__main__":
    main()
