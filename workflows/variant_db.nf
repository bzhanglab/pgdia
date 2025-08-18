workflow generate_variant_db {
    take:
    sample_id
    annotated_vcf

    main:
    var_peptides = gen_var_db(sample_id, annotated_vcf)
    mod_var_peptides(sample_id, annotated_vcf, var_peptides)

    emit:
    out = mod_var_peptides.out
}

process gen_var_db {
    input:
    val sample_id
    path annotated_vcf

    output:
    path "results/variant_db/${sample_id}_var_peptides.fa"

    script:
    """
    # run the CLI tool to generate protein DB from VCF
    python3 -m pypgatk.pypgatk_cli vcf-to-proteindb \\
      --vcf ${annotated_vcf} \\
      --input_fasta GENCODE_V42_reference_v1.1.1/genome/input_fasta.fa \\
      --gene_annotations_gtf GENCODE_V42_reference_v1.1.1/gene_annotation/GENCODE.V42.basic.CHR.primary.selection.gtf \\
      --include_consequences missense_variant,frameshift_variant,inframe_insertion,inframe_deletion \\
      --af_field AF \\
      --af_threshold 0.05 \\
      -o results/variant_db/${sample_id}_var_peptides.fa

    """
}

process mod_var_peptides {
    input:
    val sample_id
    path annotated_vcf
    path var_peptides

    output:
    path "results/variant_db/${sample_id}_var_modified_peptides.fa"

    script:
    """
    python3 local/get_var_aa_change.py \\
      results/variant_annotation/${sample_id}.haplotypecaller.filtered_VEP.ann.vcf \\
      results/variant_db/${sample_id}_var_peptides.fa \\
      results/variant_db/${sample_id}_var_modified_peptides.fa
    """
}
