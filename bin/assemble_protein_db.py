#!/usr/bin/env python3
import sys


def parse_fasta(path: str):
    """
    Streaming FASTA parser.
    Yields (header_line_without_>, seq, description_after_first_space).
    - header_line_without_> is the full header (up to end of line) without the leading '>'.
    - seq is the concatenated sequence (no whitespace).
    - description_after_first_space is the part after the first space (or "" if none).
    """
    header = None
    desc = ""
    seq_chunks = []

    with open(path, "r") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                # flush previous record
                if header is not None:
                    yield header, "".join(seq_chunks), desc
                full_header = line[1:].strip()
                if not full_header:
                    raise ValueError(f"Empty FASTA header in {path}")
                parts = full_header.split(None, 1)
                header = parts[0] if parts else full_header
                desc = parts[1] if len(parts) == 2 else ""
                seq_chunks = []
            else:
                if header is None:
                    raise ValueError(f"FASTA sequence encountered before any header in {path}")
                seq_chunks.append(line.replace(" ", "").replace("\t", ""))

    # flush last record
    if header is not None:
        yield header, "".join(seq_chunks), desc


def write_fasta_record(out_fh, header: str, desc: str, seq: str, wrap: int = 60) -> None:
    if desc:
        out_fh.write(f">{header} {desc}\n")
    else:
        out_fh.write(f">{header}\n")

    if wrap and wrap > 0:
        for i in range(0, len(seq), wrap):
            out_fh.write(seq[i : i + wrap] + "\n")
    else:
        out_fh.write(seq + "\n")


def count_and_copy_fasta(in_path: str, out_fh) -> int:
    count = 0
    for header, seq, desc in parse_fasta(in_path):
        write_fasta_record(out_fh, header, desc, seq)
        count += 1
    return count


def main() -> int:
    # Keep the same argv layout as your original script:
    # sys.argv[1] Variant DB
    # sys.argv[2] novel isoforms DB
    # sys.argv[3] combined DB output
    # sys.argv[4] novel DB output (variant + isoforms)
    # sys.argv[5] reference protein DB
    if len(sys.argv) < 6:
        print(
            "Usage:\n"
            "  script.py <variant_fa> <isoform_fa> <combined_out_fa> <novel_out_fa> <ref_fa>",
            file=sys.stderr,
        )
        return 2

    fasta_file_2 = sys.argv[1]  # Variant DB
    fasta_file_3 = sys.argv[2]  # novel isoforms DB
    combined_fasta_path = sys.argv[3]  # combined DB
    novel_fasta_path = sys.argv[4]  # novel DB (variant + isoforms)
    fasta_file_1 = sys.argv[5]  # reference protein DB

    record_counts = {}
    total_records = 0

    # Write combined: ref + variant + isoforms
    with open(combined_fasta_path, "w") as out_fh:
        for fasta_file in [fasta_file_1, fasta_file_2, fasta_file_3]:
            cnt = count_and_copy_fasta(fasta_file, out_fh)
            record_counts[fasta_file] = cnt
            total_records += cnt

    for fasta_file, cnt in record_counts.items():
        print(f"{fasta_file}: {cnt} records")
    print(f"Total number of records: {total_records}")

    # Write novel: variant + isoforms
    novel_total = 0
    with open(novel_fasta_path, "w") as out_fh:
        for fasta_file in [fasta_file_2, fasta_file_3]:
            novel_total += count_and_copy_fasta(fasta_file, out_fh)

    print(f"Total number of novel records: {novel_total}")
    return 0
    

if __name__ == "__main__":
    raise SystemExit(main())