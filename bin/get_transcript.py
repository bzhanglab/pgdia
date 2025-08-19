import sys

id_file_path = sys.argv[1]
gtf_file_path = sys.argv[2]  # the GTF contains all isoforms info
output_path = sys.argv[3]  # The GTF file contains novel isoforms info only
# Load IDs into a set
# with open("../rna-seq-results/stringtie/ids_list.txt", "r") as f:
with open(id_file_path, "r") as f:
    ids = set(line.strip() for line in f)

# Open GTF file and filter by transcript_id
# with open("../rna-seq-results/stringtie/s_74YA27.stringtie_output.gtf", "r") as gtf,
# open("../rna-seq-results/stringtie/s_74YA27_novel_isoforms.gtf", "w") as output:
with open(gtf_file_path, "r") as gtf, open(output_path, "w") as output:
    for line in gtf:
        if line.startswith("#"):
            continue
        if 'transcript_id' in line:
            tid = line.split('transcript_id "')[1].split('"')[0]
            if tid in ids:
                output.write(line)