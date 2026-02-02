from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys

# Load the VEP annotated VCF file
# vcf_file_path = "../rna-seq-results/annotation/haplotypecaller/s_74YA27/s_74YA27.haplotypecaller.filtered_VEP.ann.vcf"
vcf_file_path = sys.argv[1]
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


# Read the VCF file to extract variant info + Source annotation
vcf_entries = []
aa_changes = []
with open(vcf_file_path, 'r') as file:
    for line in file:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 8:
            continue

        chrom, pos, var_id, ref, alt, qual, fltr, info = fields[:8]

        # Parse CSQ and Source fields from INFO
        info_fields = info.split(';')
        csq_field = [x for x in info_fields if x.startswith('CSQ=')]
        source_field = [x for x in info_fields if x.startswith('Source=')]
        source = source_field[0].replace('Source=', '') if source_field else "NA"

        if csq_field:
            csq_entries = csq_field[0].replace('CSQ=', '').split(',')

            for csq_entry in csq_entries:
                csq_parts = csq_entry.split('|')
                if len(csq_parts) < 16:
                    continue  # skip malformed or incomplete entries

                consequence = csq_parts[1]
                if var_id == "rs3829740":
                    print(consequence)
                if any(keyword in consequence for keyword in
                       ["missense_variant", "frameshift_variant", "inframe_insertion", "inframe_deletion"]):
                    gene_name = csq_parts[3]
                    transcript_id = csq_parts[6]
                    protein_position = csq_parts[14]
                    amino_acids = csq_parts[15]

                    if amino_acids and '/' in amino_acids:
                        ref_aa, alt_aa = amino_acids.split('/')
                    elif amino_acids and len(amino_acids) >= 2:
                        ref_aa = amino_acids[0]
                        alt_aa = amino_acids[1:]
                    else:
                        ref_aa, alt_aa = '-', '-'

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
                        'Source': source
                    })

# translated_fasta_path = "../rna-seq-results/annotation/s_74YA27_var_peptides.fa"
translated_fasta_path = sys.argv[2]
reference_fasta_path = "./rna-seq-results/combined_protein_db/GENCODE.V42.basic.CHR.combined_contaminants.fa"
reference_records = SeqIO.parse(reference_fasta_path, "fasta")
reference_proteins = {}

for record in reference_records:
    transcript_id = record.id.split("|")[1]

    # Convert Seq object to a string
    protein_sequence = str(record.seq)

    reference_proteins[transcript_id] = protein_sequence

modified_protein_records = []

# Creating a dictionary from the extracted aa_changes to make the lookup easier
aa_changes_dict = {}
for change in aa_changes:
    key = f"var_{change['Variant ID']}_{change['Chromosome']}.{change['Position']}.{change['Reference']}.{change['Alternate']}_{change['Transcript ID']}"
    aa_changes_dict[key] = change
# print(aa_changes_dict["var_._chr1.970887.C.T_ENST00000379410"])

mis_translation_cnt = 0
# Read the original translated FASTA file
with open(translated_fasta_path, "r") as fasta_file:
    protein_records = SeqIO.parse(fasta_file, "fasta")

    for record in protein_records:
        # Match the entry name with our extracted amino acid change information
        if record.id in aa_changes_dict:
            change_info = aa_changes_dict[record.id]
            source = change_info.get("Source", "NA")  # ✅ fallback to NA if not found
            # Adding amino acid change information to the FASTA entry name
            if source != "NA":
                modified_id = f"{record.id}_{change_info['Gene Name']}_{change_info['Reference AA']}{change_info['Protein Position']}{change_info['Alternate AA']}_Source={source}"
                modified_description = f"{record.description} | Amino Acid Change: {change_info['Reference AA']} -> {change_info['Alternate AA']} at position {change_info['Protein Position']} | Source: {source}"
            else:
                modified_id = f"{record.id}_{change_info['Gene Name']}_{change_info['Reference AA']}{change_info['Protein Position']}{change_info['Alternate AA']}"
                modified_description = f"{record.description} | Amino Acid Change: {change_info['Reference AA']} -> {change_info['Alternate AA']} at position {change_info['Protein Position']}"

            # Creating a new SeqRecord with modified id and description
            new_protein_seq = str(record.seq)
            if not str(record.seq).startswith("M"):
                if change_info["Transcript ID"] not in reference_proteins.keys():
                    continue
                mis_translation_cnt += 1
                protein_seq = reference_proteins[change_info["Transcript ID"]]
                protein_pos = change_info['Protein Position'].split("/")[0].split("-")
                if len(protein_pos) == 1:
                    new_protein_seq = protein_seq[:int(protein_pos[0]) - 1] + change_info['Alternate AA'] + \
                                      protein_seq[int(protein_pos[0]):]
                else:
                    new_protein_seq = protein_seq[:int(protein_pos[0]) - 1] + change_info['Alternate AA'] + \
                                      protein_seq[int(protein_pos[1]):]
            modified_record = SeqRecord(seq=new_protein_seq, id=modified_id, description=modified_description)
            modified_protein_records.append(modified_record)
        else:
            modified_protein_records.append(record)

# Write the modified protein sequences to a new FASTA file
# modified_fasta_path = "../rna-seq-results/annotation/s_74YA27_var_modified_peptides.fa"
modified_fasta_path = sys.argv[3]
with open(modified_fasta_path, "w") as output_handle:
    SeqIO.write(modified_protein_records, output_handle, "fasta")

print(f"The number of missed translation records is {mis_translation_cnt}")
