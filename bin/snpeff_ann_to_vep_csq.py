#!/usr/bin/env python3
import argparse
import re


CSQ_FIELDS = [
    "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type",
    "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp",
    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
    "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS",
    "VARIANT_CLASS", "SYMBOL_SOURCE", "HGNC_ID", "CANONICAL", "MANE",
    "MANE_SELECT", "MANE_PLUS_CLINICAL", "TSL", "APPRIS", "CCDS", "ENSP",
    "SWISSPROT", "TREMBL", "UNIPARC", "UNIPROT_ISOFORM", "GENE_PHENO",
    "DOMAINS", "miRNA", "AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF",
    "SAS_AF", "gnomADe_AF", "gnomADe_AFR_AF", "gnomADe_AMR_AF",
    "gnomADe_ASJ_AF", "gnomADe_EAS_AF", "gnomADe_FIN_AF",
    "gnomADe_MID_AF", "gnomADe_NFE_AF", "gnomADe_REMAINING_AF",
    "gnomADe_SAS_AF", "gnomADg_AF", "gnomADg_AFR_AF", "gnomADg_AMI_AF",
    "gnomADg_AMR_AF", "gnomADg_ASJ_AF", "gnomADg_EAS_AF",
    "gnomADg_FIN_AF", "gnomADg_MID_AF", "gnomADg_NFE_AF",
    "gnomADg_REMAINING_AF", "gnomADg_SAS_AF", "MAX_AF", "MAX_AF_POPS",
    "FREQS", "CLIN_SIG", "SOMATIC", "PHENO", "PUBMED", "MOTIF_NAME",
    "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE",
    "TRANSCRIPTION_FACTORS",
]

AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",
}


def parse_info(info):
    fields = []
    values = {}
    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            fields.append((key, value))
            values[key] = value
        else:
            fields.append((item, None))
            values[item] = None
    return fields, values


def aa_change_from_hgvsp(hgvsp):
    hgvsp = hgvsp.replace("p.", "")
    match = re.match(r"([A-Z][a-z]{2}|Ter)\d+([A-Z][a-z]{2}|Ter)$", hgvsp)
    if not match:
        return ""
    return f"{AA3_TO_1.get(match.group(1), match.group(1))}/{AA3_TO_1.get(match.group(2), match.group(2))}"


def ann_to_csq(ann_entry):
    parts = ann_entry.split("|")
    parts += [""] * (16 - len(parts))

    csq = [""] * len(CSQ_FIELDS)
    csq[0] = parts[0]   # Allele
    csq[1] = parts[1]   # Consequence
    csq[2] = parts[2]   # IMPACT
    csq[3] = parts[3]   # SYMBOL
    csq[4] = parts[4]   # Gene
    csq[5] = parts[5]   # Feature_type
    csq[6] = parts[6]   # Feature
    csq[7] = parts[7]   # BIOTYPE
    csq[8] = parts[8]   # EXON/Rank
    csq[10] = parts[9]  # HGVSc
    csq[11] = parts[10] # HGVSp
    csq[12] = parts[11] # cDNA_position
    csq[13] = parts[12] # CDS_position
    csq[14] = parts[13] # Protein_position
    csq[15] = aa_change_from_hgvsp(parts[10])
    csq[18] = parts[14] # DISTANCE
    return "|".join(csq)


def convert_vcf(input_path, output_path):
    csq_header = (
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations '
        'converted from snpEff ANN for pypgatk compatibility. Format: '
        + "|".join(CSQ_FIELDS)
        + '">'
    )
    saw_csq_header = False
    inserted_csq_header = False

    with open(input_path, "r", encoding="utf-8") as inp, open(output_path, "w", encoding="utf-8") as out:
        for raw_line in inp:
            line = raw_line.rstrip("\n")
            if line.startswith("##INFO=<ID=CSQ,"):
                saw_csq_header = True
            if line.startswith("#CHROM") and not saw_csq_header and not inserted_csq_header:
                out.write(csq_header + "\n")
                inserted_csq_header = True

            if line.startswith("#"):
                out.write(line + "\n")
                continue

            cols = line.split("\t")
            if len(cols) < 8:
                out.write(line + "\n")
                continue

            info_fields, info_map = parse_info(cols[7])
            ann = info_map.get("ANN")
            if ann and "CSQ" not in info_map:
                csq_value = ",".join(ann_to_csq(entry) for entry in ann.split(",") if entry)
                info_fields.append(("CSQ", csq_value))
                cols[7] = ";".join(
                    key if value is None else f"{key}={value}"
                    for key, value in info_fields
                )

            out.write("\t".join(cols) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Convert snpEff ANN annotations to VEP-like CSQ annotations.")
    parser.add_argument("--input", required=True, help="Input snpEff-annotated VCF")
    parser.add_argument("--output", required=True, help="Output VCF with CSQ annotations")
    args = parser.parse_args()
    convert_vcf(args.input, args.output)


if __name__ == "__main__":
    main()
