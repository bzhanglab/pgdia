from Bio import SeqIO
import sys

# Paths to the three different FASTA files
# fasta_file_1 = "../rna-seq-results/combined_protein_db/GENCODE.V42.basic.CHR.combined_contaminants.fa"
# fasta_file_2 = "../rna-seq-results/annotation/s_74YA27_var_modified_peptides.fa"
# fasta_file_3 = "../rna-seq-results/stringtie/s_74YA27_novel_transcripts.fasta.transdecoder.pep"

# fasta_file_1 = "../rna-seq-results/combined_protein_db/GENCODE.V42.basic.CHR.combined_contaminants.fa"
fasta_file_1 = sys.argv[5]  # reference protein DB: GENCODE.V42.basic.CHR.combined_contaminants.fa 
fasta_file_2 = sys.argv[1]  # Variant DB
fasta_file_3 = sys.argv[2]  # novel isoforms DB
combined_fasta_path = sys.argv[3]  # combined DB
novel_fasta_path = sys.argv[4]

# Output path for the combined FASTA file
# combined_fasta_path = "../rna-seq-results/combined_protein_db/s_74YA27_combined_protein_db.fa"

# List to store all sequences from the three files
combined_records = []

record_counts = {}
total_records = 0
# Process each file
for fasta_file in [fasta_file_1, fasta_file_2, fasta_file_3]:
    count = 0
    with open(fasta_file, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            combined_records.append(record)
            count += 1
    record_counts[fasta_file] = count
    total_records += count


# Write all combined records to the new output FASTA file
with open(combined_fasta_path, "w") as output_handle:
    SeqIO.write(combined_records, output_handle, "fasta")

for fasta_file, count in record_counts.items():
    print(f"{fasta_file}: {count} records")
print(f"Total number of records: {total_records}")

novel_records = []
for fasta_file in [fasta_file_2, fasta_file_3]:
    with open(fasta_file, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            novel_records.append(record)

# Write all combined records to the new output FASTA file
with open(novel_fasta_path, "w") as output_handle:
    SeqIO.write(novel_records, output_handle, "fasta")
print(f"Total number of novel records: {len(novel_records)}")