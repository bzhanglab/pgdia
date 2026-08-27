#!/usr/bin/env python3
"""Generate the small RNA-seq and reference inputs used by the test profile."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DATA = Path(__file__).resolve().parent / "data"
TRANSCRIPT_START = 501

PEPTIDES = (
    "VLGHENR",
    "PEPTIDER",
    "SDFGHLEK",
    "VAVVQAEK",
    "FQGQVLDSLR",
    "LSEPAELTDAVK",
    "LGEYGFQNALIVR",
    "AGLCQTFVYGGCR",
    "TMTDIDGDGTITTK",
    "YLQEIYNSNNQK",
)

CODON = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT",
    "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTG",
    "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT",
    "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
}


def run(command: list[str]) -> None:
    subprocess.run(command, check=True)


def write_fasta(path: Path, name: str, sequence: str) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f">{name}\n")
        for offset in range(0, len(sequence), 80):
            handle.write(sequence[offset:offset + 80] + "\n")


def main() -> None:
    if shutil.which("samtools") is None:
        raise SystemExit("samtools is required to regenerate the BAM and reference indexes")

    DATA.mkdir(parents=True, exist_ok=True)
    protein = "M" + "K".join(PEPTIDES) + "K" + ("AGST" * 15)
    coding_sequence = "".join(CODON[aa] for aa in protein) + "TAA"
    transcript_end = TRANSCRIPT_START + len(coding_sequence) - 1
    genome = "A" * (TRANSCRIPT_START - 1) + coding_sequence + "A" * 500

    genome_fasta = DATA / "simulated_genome.fasta"
    write_fasta(genome_fasta, "chrSim", genome)
    run(["samtools", "faidx", str(genome_fasta)])
    with (DATA / "simulated_genome.dict").open("w", encoding="utf-8") as handle:
        handle.write("@HD\tVN:1.6\tSO:unsorted\n")
        handle.write(f"@SQ\tSN:chrSim\tLN:{len(genome)}\n")

    gtf = DATA / "simulated_annotation.gtf"
    attributes = 'gene_id "SIM_GENE001"; transcript_id "SIM_TX001"; gene_name "SIM_GENE001";'
    with gtf.open("w", encoding="utf-8") as handle:
        handle.write(
            f"chrSim\tPGDIA\ttranscript\t{TRANSCRIPT_START}\t{transcript_end}\t.\t+\t.\t{attributes}\n"
        )
        handle.write(
            f"chrSim\tPGDIA\texon\t{TRANSCRIPT_START}\t{transcript_end}\t.\t+\t.\t{attributes}\n"
        )
        handle.write(
            f"chrSim\tPGDIA\tCDS\t{TRANSCRIPT_START}\t{transcript_end - 3}\t.\t+\t0\t{attributes}\n"
        )

    bam = DATA / "simulated_rna.bam"
    with tempfile.TemporaryDirectory(prefix="pgdia_example_") as temp_dir_name:
        temp_dir = Path(temp_dir_name)
        sam = temp_dir / "simulated_rna.sam"
        unsorted_bam = temp_dir / "simulated_rna.unsorted.bam"
        with sam.open("w", encoding="utf-8") as handle:
            handle.write("@HD\tVN:1.6\tSO:coordinate\n")
            handle.write(f"@SQ\tSN:chrSim\tLN:{len(genome)}\n")
            handle.write("@RG\tID:SIM\tSM:SIM_DIA\tPL:ILLUMINA\n")
            read_number = 1
            for start in range(TRANSCRIPT_START, transcript_end - 74, 20):
                sequence = genome[start - 1:start - 1 + 75]
                handle.write(
                    f"SIM_READ_{read_number:04d}\t0\tchrSim\t{start}\t60\t75M\t*\t0\t0\t"
                    f"{sequence}\t{'I' * 75}\tRG:Z:SIM\n"
                )
                read_number += 1
        run(["samtools", "view", "--no-PG", "-b", "-o", str(unsorted_bam), str(sam)])
        run(["samtools", "sort", "--no-PG", "-o", str(bam), str(unsorted_bam)])
    run(["samtools", "index", str(bam)])

    variant_position = TRANSCRIPT_START + 2
    vcf = DATA / "simulated_vep_annotated.vcf"
    csq_format = (
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|"
        "EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids"
    )
    csq_value = (
        "A|missense_variant|MODERATE|SIM_GENE001|SIM_GENE001|Transcript|SIM_TX001|"
        "protein_coding|1/1||c.3G>A|p.Met1Ile|3|3|1/"
        f"{len(protein)}|M/I"
    )
    with vcf.open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write(f"##contig=<ID=chrSim,length={len(genome)}>\n")
        handle.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n')
        handle.write(f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: {csq_format}">\n')
        handle.write('##INFO=<ID=Source,Number=1,Type=String,Description="Variant source">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSIM_DIA\n")
        handle.write(
            f"chrSim\t{variant_position}\tSIM_VAR_001\tG\tA\t60\tPASS\t"
            f"AF=0.50;CSQ={csq_value};Source=simulated\tGT\t0/1\n"
        )

    relative_data = "examples/runnable_example/data"
    samplesheet = Path(__file__).resolve().parent / "samplesheet.csv"
    with samplesheet.open("w", encoding="utf-8") as handle:
        handle.write("sample,fastq_1,fastq_2,bam,bai,cram,crai,vcf,tbi,dia_raw,strandedness\n")
        handle.write(
            f"SIM_DIA,,,{relative_data}/simulated_rna.bam,"
            f"{relative_data}/simulated_rna.bam.bai,,,,,"
            f"{relative_data}/simulated_dia.mzML,unstranded\n"
        )

    print(f"Wrote workflow inputs under {DATA.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
