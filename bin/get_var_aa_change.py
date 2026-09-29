import argparse
import re
import sys

TARGET_CONSEQUENCES = {
    "missense_variant",
    "frameshift_variant",
    "inframe_insertion",
    "inframe_deletion",
}

AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Annotate variant peptide FASTA headers with AA-change metadata.")
    parser.add_argument("vcf")
    parser.add_argument("translated_fasta")
    parser.add_argument("reference_fasta")
    parser.add_argument("output_fasta")
    parser.add_argument("--sample-id", default=None)
    parser.add_argument("--tissue", default=None)
    parser.add_argument("--organism", required=True)
    parser.add_argument("--taxonomy-id", required=True)
    return parser.parse_args()


def safe_token(value):
    cleaned = re.sub(r"[^A-Za-z0-9_.:/-]+", "_", str(value).strip())
    cleaned = cleaned.strip("_.:/-")
    return cleaned or "UNKNOWN"


def parse_fasta(path):
    record_id = None
    description = None
    sequence_chunks = []

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if record_id is not None:
                    yield record_id, description, "".join(sequence_chunks)
                description = line[1:].strip()
                if not description:
                    raise ValueError(f"Empty FASTA header in {path}")
                record_id = description.split(None, 1)[0]
                sequence_chunks = []
            else:
                if record_id is None:
                    raise ValueError(f"Sequence found before the first FASTA header in {path}")
                sequence_chunks.append(line.replace(" ", "").replace("\t", ""))

    if record_id is not None:
        yield record_id, description, "".join(sequence_chunks)


def write_fasta(records, path, wrap=60):
    with open(path, "w", encoding="utf-8") as handle:
        for description, sequence in records:
            handle.write(f">{description}\n")
            for start in range(0, len(sequence), wrap):
                handle.write(sequence[start : start + wrap] + "\n")


def reference_identifiers(record_id, description):
    identifiers = {record_id}
    identifiers.update(part for part in record_id.split("|") if part)
    for pattern in (
        r"(?:^|\s)Transcript=([^\s;]+)",
        r"(?:^|\s)transcript_id=([^\s;]+)",
        r"(?:^|\s)RefSeq=([^\s;]+)",
    ):
        identifiers.update(re.findall(pattern, description))
    return identifiers


def parse_info_field(info):
    parsed = {}
    for field in info.split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            parsed[key] = value
    return parsed


def aa_from_vep(amino_acids):
    if amino_acids and "/" in amino_acids:
        return amino_acids.split("/", 1)
    if amino_acids and len(amino_acids) >= 2:
        return amino_acids[0], amino_acids[1:]
    return "-", "-"


def aa_from_snpeff_hgvsp(hgvsp):
    hgvsp = hgvsp.replace("p.", "")
    if not hgvsp or hgvsp in {".", "?"}:
        return "-", "-"

    match = re.match(r"([A-Z][a-z]{2}|Ter)(\d+)([A-Z][a-z]{2}|Ter)$", hgvsp)
    if match:
        return AA3_TO_1.get(match.group(1), match.group(1)), AA3_TO_1.get(match.group(3), match.group(3))

    match = re.match(r"([A-Z][a-z]{2}|Ter)(\d+).*(fs|del|ins|dup)", hgvsp)
    if match:
        return AA3_TO_1.get(match.group(1), match.group(1)), match.group(3)

    return "-", "-"


def add_aa_change(aa_changes, chrom, pos, var_id, ref, alt, gene_name, transcript_id,
                  protein_position, ref_aa, alt_aa, source, consequence):
    aa_changes.append({
        'Gene Name': gene_name,
        'Transcript ID': transcript_id,
        'Protein Position': protein_position,
        'Reference AA': ref_aa,
        'Alternate AA': alt_aa,
        'Chromosome': chrom,
        'Position': pos,
        'Reference': ref,
        'Alternate': alt,
        'Variant ID': var_id,
        'Source': source,
        'Consequence': consequence
    })


def correct_missense_sequence(change, reference_sequence):
    """Reconstruct only a verified single-residue missense substitution."""
    consequences = set(change['Consequence'].split('&')) & TARGET_CONSEQUENCES
    if consequences != {'missense_variant'} or not reference_sequence:
        return None
    ref_aa, alt_aa = change['Reference AA'], change['Alternate AA']
    amino_acids = set(AA3_TO_1.values()) - {'*'}
    if ref_aa not in amino_acids or alt_aa not in amino_acids:
        return None
    position = change['Protein Position'].split('/')[0]
    if not position.isdigit():
        return None
    index = int(position) - 1
    if not 0 <= index < len(reference_sequence) or reference_sequence[index] != ref_aa:
        return None
    return reference_sequence[:index] + alt_aa + reference_sequence[index + 1:]


args = parse_args()
vcf_file_path = args.vcf
# # Read the VCF file to extract information
# vcf_entries = []
# with open(vcf_file_path, 'r') as file:
#     for line in file:
#         if not line.startswith('#'):
#             vcf_entries.append(line.strip().split('\t'))
#
# # Extracting the CSQ field from the VCF entries and parsing for amino acid changes
# aa_changes = []
#
# for entry in vcf_entries:
#     info_field = entry[7]
#     # Extract the CSQ field
#     csq_field = [field for field in info_field.split(';') if field.startswith('CSQ=')]
#     source_field = [x for x in info.split(';') if x.startswith('Source=')]
#
#     source = source_field[0].replace('Source=', '') if source_field else "NA"
#
#     if csq_field:
#         csq_field = csq_field[0].replace('CSQ=', '')
#         fields = csq_field.split('|')
#
#         # We are interested in missense variants
#         consequence = fields[1]
#         if 'missense_variant' in consequence or 'frameshift_variant' in consequence or \
#                 'inframe_insertion' in consequence or 'inframe_deletion' in consequence:
#             gene_name = fields[3]
#             transcript_id = fields[6]
#             protein_position = fields[14]
#             amino_acids = fields[15]
#
#             if amino_acids and '/' in amino_acids:
#                 # If we have a '/' delimiter, split as before
#                 ref_aa, alt_aa = amino_acids.split('/')
#             else:
#                 # Handle cases where no '/' is present
#                 # For frameshifts, you might have just one amino acid before shift
#                 # You can define logic here, for example:
#                 ref_aa = amino_acids[0] if amino_acids else '-'
#                 alt_aa = amino_acids[2:] if len(amino_acids) > 2 else (amino_acids if amino_acids else '-')
#
#             aa_changes.append({
#                     'Gene Name': gene_name,
#                     'Transcript ID': transcript_id,
#                     'Protein Position': protein_position,
#                     'Reference AA': ref_aa,
#                     'Alternate AA': alt_aa,
#                     'Chromosome': entry[0],
#                     'Position': entry[1],
#                     'Reference': entry[3],
#                     'Alternate': entry[4],
#                     'Variant ID': entry[2]
#                 })


# Read the VCF file to extract variant info + Source annotation.
aa_changes = []
with open(vcf_file_path, 'r') as file:
    for line in file:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 8:
            continue

        chrom, pos, var_id, ref, alt, qual, fltr, info = fields[:8]
        info_map = parse_info_field(info)
        source = info_map.get("Source", "NA")

        for csq_entry in info_map.get("CSQ", "").split(","):
            if not csq_entry:
                continue
            csq_parts = csq_entry.split('|')
            if len(csq_parts) < 16:
                continue

            consequence = csq_parts[1]
            if not any(keyword in consequence for keyword in TARGET_CONSEQUENCES):
                continue

            ref_aa, alt_aa = aa_from_vep(csq_parts[15])
            add_aa_change(
                aa_changes,
                chrom,
                pos,
                var_id,
                ref,
                alt,
                csq_parts[3],
                csq_parts[6],
                csq_parts[14],
                ref_aa,
                alt_aa,
                source,
                consequence,
            )

        for ann_entry in info_map.get("ANN", "").split(","):
            if not ann_entry:
                continue
            ann_parts = ann_entry.split('|')
            if len(ann_parts) < 14:
                continue

            consequence = ann_parts[1]
            if not any(keyword in consequence for keyword in TARGET_CONSEQUENCES):
                continue

            ref_aa, alt_aa = aa_from_snpeff_hgvsp(ann_parts[10])
            add_aa_change(
                aa_changes,
                chrom,
                pos,
                var_id,
                ref,
                ann_parts[0] or alt,
                ann_parts[3],
                ann_parts[6],
                ann_parts[13],
                ref_aa,
                alt_aa,
                source,
                consequence,
            )

translated_fasta_path = args.translated_fasta
reference_fasta_path = args.reference_fasta
reference_proteins = {}

for record_id, description, protein_sequence in parse_fasta(reference_fasta_path):
    for identifier in reference_identifiers(record_id, description):
        reference_proteins.setdefault(identifier, protein_sequence)

modified_protein_records = []

# Creating a dictionary from the extracted aa_changes to make the lookup easier
aa_changes_dict = {}
for change in aa_changes:
    key = f"var_{change['Variant ID']}_{change['Chromosome']}.{change['Position']}.{change['Reference']}.{change['Alternate']}_{change['Transcript ID']}"
    aa_changes_dict[key] = change
# print(aa_changes_dict["var_._chr1.970887.C.T_ENST00000379410"])

corrected_non_m_cnt = 0
uncorrected_non_m_cnt = 0
# Read the original translated FASTA file
for record_id, _original_description, record_sequence in parse_fasta(translated_fasta_path):
    # Match the entry name with our extracted amino acid change information
    if record_id in aa_changes_dict:
        change_info = aa_changes_dict[record_id]
        source = change_info.get("Source", "NA")
        # Adding amino acid change information to the FASTA entry name
        gene_name = safe_token(change_info['Gene Name'])
        variant_label = safe_token(
            f"{record_id}_{gene_name}_{change_info['Reference AA']}{change_info['Protein Position']}{change_info['Alternate AA']}"
        )
        sample_field = f" Sample={safe_token(args.sample_id)}" if args.sample_id else ""
        tissue_field = f" Tissue={safe_token(args.tissue)}" if args.tissue else ""
        species_fields = (
            f"OS={args.organism} OX={safe_token(args.taxonomy_id)} "
            f"GN={gene_name} PE=5 SV=1{sample_field}{tissue_field}"
        )
        if source != "NA":
            modified_id = f"sv|{variant_label}_Source={safe_token(source)}|{gene_name}"
            modified_description = (
                f"{modified_id} {species_fields} OriginalID={safe_token(record_id)} "
                f"AminoAcidChange={change_info['Reference AA']}>{change_info['Alternate AA']} "
                f"ProteinPosition={change_info['Protein Position']} Source={safe_token(source)}"
            )
        else:
            modified_id = f"sv|{variant_label}|{gene_name}"
            modified_description = (
                f"{modified_id} {species_fields} OriginalID={safe_token(record_id)} "
                f"AminoAcidChange={change_info['Reference AA']}>{change_info['Alternate AA']} "
                f"ProteinPosition={change_info['Protein Position']}"
            )

        new_protein_seq = record_sequence
        if not record_sequence.startswith("M"):
            protein_seq = reference_proteins.get(change_info["Transcript ID"])
            corrected = correct_missense_sequence(change_info, protein_seq)
            if corrected is None:
                # Preserve translated indels/frameshifts and any substitution
                # that cannot be verified against the reference protein.
                uncorrected_non_m_cnt += 1
            else:
                corrected_non_m_cnt += 1
                new_protein_seq = corrected
        modified_protein_records.append((modified_description, new_protein_seq))
    else:
        original_id = safe_token(record_id)
        sample_field = f" Sample={safe_token(args.sample_id)}" if args.sample_id else ""
        tissue_field = f" Tissue={safe_token(args.tissue)}" if args.tissue else ""
        modified_id = f"sv|{original_id}|UNKNOWN"
        modified_description = (
            f"{modified_id} OS={args.organism} OX={safe_token(args.taxonomy_id)} "
            f"GN=UNKNOWN PE=5 SV=1{sample_field}{tissue_field} OriginalID={original_id}"
        )
        modified_protein_records.append((modified_description, record_sequence))

# Write the modified protein sequences to a new FASTA file
# modified_fasta_path = "../rna-seq-results/annotation/s_74YA27_var_modified_peptides.fa"
modified_fasta_path = args.output_fasta
write_fasta(modified_protein_records, modified_fasta_path)

print(f"Corrected {corrected_non_m_cnt} non-M variant protein records from a matching reference.")
if uncorrected_non_m_cnt:
    print(
        f"WARNING: Preserved {uncorrected_non_m_cnt} non-M variant protein records because "
        "a single-residue missense correction could not be verified against the reference.",
        file=sys.stderr,
    )
