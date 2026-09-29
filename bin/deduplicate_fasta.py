#!/usr/bin/env python3
import argparse
import re
from pathlib import Path


def parse_fasta(path):
    header = None
    chunks = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks)


def write_record(handle, header, sequence, wrap):
    handle.write(f">{header}\n")
    for start in range(0, len(sequence), wrap):
        handle.write(sequence[start : start + wrap] + "\n")


def source_from_path(path):
    name = Path(path).name
    for suffix in [".fa", ".fasta", ".faa", ".pep"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    name = re.sub(r"_var_modified_peptides$", "", name)
    name = re.sub(r"_var_peptides$", "", name)
    return safe_token(name)


def safe_token(value):
    cleaned = re.sub(r"[^A-Za-z0-9_.:/,-]+", "_", str(value).strip())
    cleaned = cleaned.strip("_.:/,-")
    return cleaned or "UNKNOWN"


def value_from_header(header, key, fallback=None):
    match = re.search(rf"(?:^|\s){re.escape(key)}=([^\s]+)", header)
    if match:
        return safe_token(match.group(1))
    return fallback


def sample_from_header(header, fallback):
    return value_from_header(header, "Sample", fallback)


def tissue_from_header(header):
    return value_from_header(header, "Tissue")


def strip_aggregate_fields(header):
    header = re.sub(r"\s+ObservedIn=[^\s]+", "", header)
    header = re.sub(r"\s+ObservedInRuns=[^\s]+", "", header)
    header = re.sub(r"\s+ObservedInTissues=[^\s]+", "", header)
    return header.strip()


def add_aggregate_fields(header, samples, tissues):
    clean_header = strip_aggregate_fields(header)
    fields = []
    if samples:
        fields.append(f"ObservedInRuns={','.join(samples)}")
    if tissues:
        fields.append(f"ObservedInTissues={','.join(tissues)}")
    return f"{clean_header} {' '.join(fields)}".rstrip()


def _old_sample_from_header(header, fallback):
    match = re.search(r"(?:^|\s)Sample=([^\s]+)", header)
    return safe_token(match.group(1)) if match else fallback


def main():
    parser = argparse.ArgumentParser(description="Merge FASTA files and remove duplicate sequences.")
    parser.add_argument("--output", required=True)
    parser.add_argument("--dedup-by", choices=["sequence", "header"], default="sequence")
    parser.add_argument(
        "--min-samples",
        type=int,
        default=1,
        help="Keep only records observed in at least this many input samples/runs after deduplication.",
    )
    parser.add_argument("inputs", nargs="+")
    args = parser.parse_args()

    records = []
    index_by_key = {}
    skipped = 0
    for fasta in args.inputs:
        fallback_source = source_from_path(fasta)
        for header, sequence in parse_fasta(fasta):
            key = sequence if args.dedup_by == "sequence" else header.split(None, 1)[0]
            sample = sample_from_header(header, fallback_source)
            tissue = tissue_from_header(header)
            if key in index_by_key:
                idx = index_by_key[key]
                if sample not in records[idx]["samples"]:
                    records[idx]["samples"].append(sample)
                if tissue and tissue not in records[idx]["tissues"]:
                    records[idx]["tissues"].append(tissue)
                skipped += 1
                continue
            index_by_key[key] = len(records)
            records.append({
                "header": header,
                "sequence": sequence,
                "samples": [sample],
                "tissues": [tissue] if tissue else [],
            })

    with open(args.output, "w", encoding="utf-8") as out_handle:
        for record in records:
            if len(record["samples"]) < args.min_samples:
                continue
            header = add_aggregate_fields(record["header"], sorted(record["samples"]), sorted(record["tissues"]))
            write_record(out_handle, header, record["sequence"], 60)

    kept = sum(1 for record in records if len(record["samples"]) >= args.min_samples)
    filtered = len(records) - kept
    print(
        f"Kept {kept} FASTA records; removed {skipped} duplicate records; "
        f"filtered {filtered} records with fewer than {args.min_samples} samples/runs."
    )


if __name__ == "__main__":
    main()
