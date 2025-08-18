nextflow.enable.dsl=2

process combine_protein_dbs {
    tag "${sample_id}"

    input:
        val sample_id
        path var_modified_peptides
        path novel_pep

    output:
        path "${sample_id}_combined_protein_db.fa"
        path "r${sample_id}_novel_protein_db.fa"

    script:
    """
    python3 ./src/assemble_protein_db.py \
        $var_modified_peptides \
        $novel_pep \
        ${sample_id}_combined_protein_db.fa \
        ${sample_id}_novel_protein_db.fa
    """
}
