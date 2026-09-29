#!/usr/bin/env python3
"""Sample FASTQ records and report the maximum observed read length."""

import argparse
import gzip
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Detect FASTQ read length before STAR index generation and alignment."
    )
    parser.add_argument("reads", nargs="+", help="FASTQ or FASTQ.gz input files")
    parser.add_argument(
        "--max-reads",
        type=int,
        default=10000,
        help="Maximum records to inspect in each FASTQ (default: 10000)",
    )
    parser.add_argument(
        "--expected-length",
        type=int,
        default=None,
        help="Fail if the maximum observed length differs from this configured value",
    )
    parser.add_argument("--length-output", required=True)
    parser.add_argument("--report", required=True)
    return parser.parse_args()


def open_fastq(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def inspect_fastq(path, max_reads):
    lengths = []
    with open_fastq(path) as handle:
        while len(lengths) < max_reads:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not sequence or not plus or not quality:
                raise ValueError(f"Incomplete FASTQ record in {path}")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"Malformed FASTQ record in {path}")
            sequence = sequence.rstrip("\r\n")
            quality = quality.rstrip("\r\n")
            if len(sequence) != len(quality):
                raise ValueError(
                    f"Sequence/quality length mismatch in {path}: "
                    f"{len(sequence)} != {len(quality)}"
                )
            lengths.append(len(sequence))

    if not lengths:
        raise ValueError(f"No FASTQ records found in {path}")
    return len(lengths), min(lengths), max(lengths)


def main():
    args = parse_args()
    if args.max_reads < 1:
        raise ValueError("--max-reads must be a positive integer")
    if args.expected_length is not None and args.expected_length < 1:
        raise ValueError("--expected-length must be a positive integer")

    rows = []
    for read_path in args.reads:
        count, minimum, maximum = inspect_fastq(read_path, args.max_reads)
        rows.append((Path(read_path).name, count, minimum, maximum))

    detected = max(row[3] for row in rows)
    if args.expected_length is not None and detected != args.expected_length:
        raise ValueError(
            f"Configured read length ({args.expected_length}) does not match "
            f"the maximum observed FASTQ read length ({detected}). Remove "
            "--read_length to auto-detect it or provide the observed value."
        )

    with open(args.length_output, "w", encoding="utf-8") as handle:
        handle.write(f"{detected}\n")

    with open(args.report, "w", encoding="utf-8") as handle:
        handle.write("fastq\treads_checked\tminimum_length\tmaximum_length\n")
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")

    print(
        f"Detected maximum read length {detected} bp from "
        f"{len(rows)} FASTQ file(s)."
    )


if __name__ == "__main__":
    main()
