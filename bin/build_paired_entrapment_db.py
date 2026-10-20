#!/usr/bin/env python3
"""Build a per-sample 1:1 peptide-level paired shuffled-entrapment database.

Random-shuffle helpers ported from build_cell_line_foreign_species_entrapment.py.
Entrapments preserve composition and the final residue; I/L-equivalent matches
anywhere in the input proteins are forbidden. Search with empty cleavage rules,
zero missed cleavages and no methionine excision. Generate search decoys separately.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
import re
from collections import Counter, deque
from pathlib import Path

from Bio import SeqIO

STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")


def fasta_records(path):
    with Path(path).open() as handle:
        yield from SeqIO.parse(handle, 'fasta')


def normalize_peptide(peptide: str, collapse_i_to_l: bool) -> str:
    return peptide.replace("I", "L") if collapse_i_to_l else peptide


def canonical_segments(sequence: str) -> list[str]:
    """Split at stop codons, gaps, ambiguity codes, and non-standard residues."""
    segments: list[str] = []
    current: list[str] = []
    for residue in sequence.upper():
        if residue in STANDARD_AA:
            current.append(residue)
        elif current:
            segments.append("".join(current))
            current = []
    if current:
        segments.append("".join(current))
    return segments


def digest_segment(
    sequence: str,
    min_length: int,
    max_length: int,
    missed_cleavages: int,
    collapse_i_to_l: bool,
) -> set[str]:
    """Digest after K/R, including K/R followed by P (no proline suppression)."""
    tryptic_parts: list[str] = []
    start = 0
    for index, residue in enumerate(sequence):
        if residue in {"K", "R"}:
            tryptic_parts.append(sequence[start : index + 1])
            start = index + 1
    if start < len(sequence):
        tryptic_parts.append(sequence[start:])

    peptides: set[str] = set()
    for start_index in range(len(tryptic_parts)):
        for missed in range(missed_cleavages + 1):
            stop_index = start_index + missed + 1
            if stop_index > len(tryptic_parts):
                break
            peptide = "".join(tryptic_parts[start_index:stop_index])
            if min_length <= len(peptide) <= max_length:
                peptides.add(normalize_peptide(peptide, collapse_i_to_l))
    return peptides


def digest_fasta(
    fasta_path: Path,
    min_length: int,
    max_length: int,
    missed_cleavages: int,
    collapse_i_to_l: bool,
) -> set[str]:
    peptides: set[str] = set()
    for record in fasta_records(fasta_path):
        for segment in canonical_segments(str(record.seq)):
            peptides.update(
                digest_segment(
                    segment,
                    min_length,
                    max_length,
                    missed_cleavages,
                    collapse_i_to_l,
                )
            )
    if not peptides:
        raise ValueError(f"No peptides generated from {fasta_path}")
    return peptides


def digest_record(record, **digest_options) -> set[str]:
    peptides: set[str] = set()
    for segment in canonical_segments(str(record.seq)):
        peptides.update(digest_segment(segment, **digest_options))
    return peptides


def canonical_transcript_gene_map(canonical_fasta: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for record in fasta_records(canonical_fasta):
        parts = record.id.split("|")
        if len(parts) >= 3 and parts[0] == "gc":
            gene = parts[2]
            transcript_match = re.search(r"(ENST\d+(?:\.\d+)?)", parts[1])
            if transcript_match:
                mapping[transcript_match.group(1).split(".")[0]] = gene
    return mapping


def source_gene(record, transcript_gene: dict[str, str]) -> str:
    description = record.description
    uniprot_gene = re.search(r"(?:^|\s)GN=([^\s]+)", description)
    if uniprot_gene:
        return uniprot_gene.group(1)

    parts = record.id.split("|")
    if len(parts) >= 3 and parts[0] == "gc":
        return parts[2]
    if record.id.startswith("STRG."):
        match = re.match(r"(STRG\.\d+\.\d+)", record.id)
        if match:
            return match.group(1)

    transcript_match = re.search(r"(ENST\d+)", record.id)
    if transcript_match:
        transcript = transcript_match.group(1)
        explicit_gene = re.search(r"ENST\d+_([A-Za-z0-9.-]+)_[A-Z*]+\d", record.id)
        if explicit_gene:
            return explicit_gene.group(1)
        if transcript in transcript_gene:
            return transcript_gene[transcript]
        return transcript
    return record.id


def safe_header_token(value: str, maximum_length: int = 100) -> str:
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")
    return (token or "UNKNOWN")[:maximum_length]


def write_annotated_fasta_entries(
    handle,
    class_prefix: str,
    selected_peptides: set[str],
    source_fasta: Path,
    transcript_gene: dict[str, str],
    digest_options: dict[str, object],
) -> tuple[int, int]:
    """Write one entry per peptide-source association with UniProt-like headers."""
    entry_index = 0
    observed_peptides: set[str] = set()
    for record in fasta_records(source_fasta):
        record_peptides = digest_record(record, **digest_options) & selected_peptides
        if not record_peptides:
            continue
        gene = safe_header_token(source_gene(record, transcript_gene))
        source_id = record.id
        source_token = safe_header_token(source_id)
        for peptide in sorted(record_peptides):
            entry_index += 1
            accession = f"{class_prefix}_{entry_index:08d}"
            entry_name = f"{source_token}_{class_prefix}"
            handle.write(
                f">sp|{accession}|{entry_name} GN={gene}\n"
                f"{peptide}\n"
            )
            observed_peptides.add(peptide)

    missing = selected_peptides - observed_peptides
    if missing:
        examples = ", ".join(sorted(missing)[:5])
        raise ValueError(
            f"{class_prefix}: {len(missing)} selected peptides could not be mapped "
            f"back to {source_fasta}; examples: {examples}"
        )
    return entry_index, len(observed_peptides)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def cell_seed(seed: int, cell_line: str) -> int:
    value = hashlib.sha256(f"{seed}:{cell_line}".encode()).digest()
    return int.from_bytes(value[:8], byteorder="big", signed=False)


def database_substring_hits(peptides: set[str], fasta_paths: list[Path]) -> set[str]:
    """Find I/L-equivalent matches anywhere in proteins using Aho-Corasick.

    Includes non-enzymatic occurrences, independently of digestion settings.
    Nonstandard residues break matches; proteins are never joined together.
    """
    edges: list[dict[str, int]] = [{}]
    failure = [0]
    terminal: list[str | None] = [None]
    output_link = [0]
    for peptide in sorted(peptides):
        node = 0
        for residue in peptide.replace("I", "L"):
            if residue not in edges[node]:
                edges[node][residue] = len(edges)
                edges.append({})
                failure.append(0)
                terminal.append(None)
                output_link.append(0)
            node = edges[node][residue]
        terminal[node] = peptide.replace("I", "L")
    queue = deque(edges[0].values())
    while queue:
        parent = queue.popleft()
        for residue, child in edges[parent].items():
            queue.append(child)
            fallback = failure[parent]
            while fallback and residue not in edges[fallback]:
                fallback = failure[fallback]
            failure[child] = edges[fallback].get(residue, 0)
            suffix = failure[child]
            output_link[child] = suffix if terminal[suffix] else output_link[suffix]
    hits: set[str] = set()
    for path in fasta_paths:
        for record in fasta_records(path):
            node = 0
            for residue in str(record.seq).upper().replace("I", "L"):
                while node and residue not in edges[node]:
                    node = failure[node]
                node = edges[node].get(residue, 0)
                match = node
                while match:
                    if terminal[match]:
                        hits.add(terminal[match])
                    match = output_link[match]
    return hits


def available_shuffle_exhaustive(target: str, blocked: set[str], used: set[str]) -> str | None:
    """Resolve small permutation spaces exactly; never call budget exhaustion proof.

    I/L-equivalent permutations are enumerated once, then literal I content is
    restored so the candidate retains the original amino-acid composition.
    """
    counts = Counter(target[:-1].replace("I", "L"))
    permutations = math.factorial(len(target) - 1)
    for count in counts.values():
        permutations //= math.factorial(count)
    if permutations > 10000:
        raise ValueError(f"{target}: random search exhausted but {permutations} permutations "
                         "remain possible; cannot prove unshufflability. Increase --max-shuffle-attempts.")

    def heads(prefix: str):
        if len(prefix) == len(target) - 1:
            yield prefix
            return
        for residue in sorted(counts):
            if counts[residue]:
                counts[residue] -= 1
                yield from heads(prefix + residue)
                counts[residue] += 1

    for head in heads(""):
        key = head + target[-1].replace("I", "L")
        if key not in blocked and key not in used:
            restored = head.replace("L", "I", target[:-1].count("I"))
            return restored + target[-1]
    return None


def shuffle_pairs(targets, forbidden, fasta_paths, seed, max_attempts):
    """Use random shuffles; report exclusions only after exhausting small spaces."""
    rng = random.Random(seed)
    blocked = {peptide.replace('I', 'L') for peptide in forbidden}
    pairs, exclusions, used = {}, {}, set()
    attempts = Counter()
    pending = sorted(targets)
    while pending:
        for target in pending:
            if len(set(target[:-1].replace('I', 'L'))) <= 1:
                exclusions[target] = 'no_distinct_permutation_after_I_L_equivalence'
                continue
            candidate = None
            while attempts[target] < max_attempts:
                attempts[target] += 1
                residues = list(target[:-1])
                rng.shuffle(residues)
                draw = ''.join(residues) + target[-1]
                key = draw.replace('I', 'L')
                if key not in blocked and key not in used:
                    candidate = draw
                    break
            if candidate is None:
                candidate = available_shuffle_exhaustive(target, blocked, used)
            if candidate is None:
                exclusions[target] = 'all_permutations_match_input_databases_or_assigned_entrapments'
                continue
            pairs[target] = candidate
            used.add(candidate.replace('I', 'L'))
        collisions = database_substring_hits({pairs[target] for target in pending if target in pairs}, fasta_paths)
        blocked.update(collisions)
        pending = sorted(target for target in pairs if pairs[target].replace('I', 'L') in collisions)
        for target in pending:
            used.remove(pairs.pop(target).replace('I', 'L'))
    return pairs, exclusions


def build(reference, novel, sample, prefix, seed=20260817, missed_cleavages=1,
          min_length=7, max_length=30, max_attempts=10000):
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', sample):
        raise ValueError('Sample must contain only letters, numbers, underscores, dots or hyphens')
    if not 0 <= missed_cleavages <= 2 or not 1 <= min_length <= max_length or max_attempts < 1:
        raise ValueError('Invalid digestion settings or shuffle attempt limit')
    paths = [Path(str(prefix) + suffix) for suffix in ('.fasta', '.pairs.tsv', '.excluded.tsv', '.manifest.json')]
    for path in paths:
        if path.exists():
            raise FileExistsError(f'Refusing to overwrite {path}')
        if not path.parent.is_dir():
            raise FileNotFoundError(path.parent)
    options = dict(min_length=min_length, max_length=max_length,
                   missed_cleavages=missed_cleavages, collapse_i_to_l=False)
    canonical = digest_fasta(reference, **options)
    novel_peptides = digest_fasta(novel, **options)
    canonical_keys = {peptide.replace('I', 'L') for peptide in canonical}
    representatives = {}
    for peptide in sorted(novel_peptides):
        key = peptide.replace('I', 'L')
        if key not in canonical_keys:
            representatives.setdefault(key, peptide)
    targets = set(representatives.values())
    if not targets:
        raise ValueError(f'{sample}: no non-reference target peptides')
    pairs, exclusions = shuffle_pairs(targets, canonical | novel_peptides,
                                      [reference, novel], cell_seed(seed, sample), max_attempts)
    if not pairs:
        raise ValueError(f'{sample}: no pairable non-reference targets remain')
    with paths[0].open('x') as fasta, paths[1].open('x') as table:
        write_annotated_fasta_entries(fasta, 'T', canonical, reference,
                                     canonical_transcript_gene_map(reference), options)
        table.write('sequence\tdecoy\tproteins\tpeptide_type\tpeptide_pair_index\n')
        for index, target in enumerate(sorted(pairs)):
            for peptide, kind in ((target, 'target'), (pairs[target], 'p_target')):
                accession = f'NOVEL_TARGET_{sample}_{index:08d}'
                if kind == 'p_target':
                    accession += '_p_target'
                protein = f'sp|{accession}|{accession}'
                fasta.write(f'>{protein}\n{peptide}\n')
                table.write(f'{peptide}\tNo\t{protein}\t{kind}\t{index}\n')
    with paths[2].open('x') as table:
        table.write('sequence\treason\n')
        for peptide, reason in sorted(exclusions.items()):
            table.write(f'{peptide}\t{reason}\n')
    manifest = dict(
        method='paired peptide shuffling', sample=sample,
        inputs={str(path): file_sha256(path) for path in (reference, novel)},
        seed=seed, sample_seed=cell_seed(seed, sample), digestion=options,
        cleavage='after K/R without proline suppression', shuffle_selection='random',
        max_shuffle_attempts=max_attempts, nr_targets=len(pairs), entrapments=len(pairs), ratio=1,
        nr_targets_before_exclusion=len(targets), excluded_nr_targets=len(exclusions),
        constraints=['exact amino-acid composition preserved', 'last residue fixed',
                     'I/L-equivalent uniqueness', 'entrapments absent anywhere in both input FASTAs'],
        diann_options=['--cut ""', '--missed-cleavages 0', 'no --met-excision', '--report-decoys'],
        scope='paired novel-target/entrapment list; reference peptides are search background',
        outputs={name: str(path) for name, path in zip(('fasta', 'pairs', 'excluded'), paths)},
    )
    with paths[3].open('x') as handle:
        json.dump(manifest, handle, indent=2)
        handle.write('\n')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-fasta', type=Path, required=True)
    parser.add_argument('--novel-fasta', type=Path, required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--output-prefix', type=Path, required=True)
    parser.add_argument('--seed', type=int, default=20260817)
    parser.add_argument('--missed-cleavages', type=int, choices=[0, 1, 2], default=1)
    parser.add_argument('--min-length', type=int, default=7)
    parser.add_argument('--max-length', type=int, default=30)
    parser.add_argument('--max-shuffle-attempts', type=int, default=10000)
    args = parser.parse_args()
    summary = build(args.reference_fasta, args.novel_fasta, args.sample, args.output_prefix,
                    args.seed, args.missed_cleavages, args.min_length, args.max_length,
                    args.max_shuffle_attempts)
    print(f"{summary['nr_targets']} novel targets; {summary['entrapments']} entrapments; "
          f"{summary['excluded_nr_targets']} excluded; ratio=1")


if __name__ == '__main__':
    main()
