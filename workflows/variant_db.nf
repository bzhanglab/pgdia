workflow generate_variant_db {
    take:
    sample_id
    annotated_vcf

    main:
    def ref_fasta = params.fasta ?: (params.genomes && params.genome && params.genomes.containsKey(params.genome) ? params.genomes[params.genome].fasta : null)
    def ref_gtf   = params.gtf   ?: (params.genomes && params.genome && params.genomes.containsKey(params.genome) ? params.genomes[params.genome].gtf   : null)

    if (!ref_fasta) {
        error("Missing required reference fasta: set --fasta or provide a genome entry with a fasta path.")
    }
    if (!ref_gtf) {
        error("Missing required reference GTF: set --gtf or provide a genome entry with a gtf path.")
    }

    var_peptides_ch = gen_var_db(sample_id, annotated_vcf)
    mod = mod_var_peptides(var_peptides_ch)

    emit:
    variant_db = mod.out
}

process gen_var_db {
  tag "${sample_id}"

    input:
    val sample_id 
    path annotated_vcf

    output:
    tuple val(sample_id),
          path(annotated_vcf),
          path("results/variant_db/${sample_id}_var_peptides.fa")

    script:
    """
    set -euo pipefail
    mkdir -p results/variant_db

    # run the CLI tool to generate protein DB from VCF
    python3 -m pypgatk.pypgatk_cli vcf-to-proteindb \
      --vcf "${annotated_vcf}" \
      --input_fasta "${ref_fasta}" \
      --gene_annotations_gtf "${ref_gtf}" \
      --include_consequences missense_variant,frameshift_variant,inframe_insertion,inframe_deletion \
      --af_field AF \
      --af_threshold 0.05 \
      -o results/variant_db/${sample_id}_var_peptides.fa

    """
}
 
process mod_var_peptides {
    tag "${sample_id}"

    input:
    tuple val(sample_id),
          path(annotated_vcf),
          path(var_peptides)

    output:
    path("results/variant_db/${sample_id}_var_modified_peptides.fa")

    script:
    """
    set -euo pipefail
    mkdir -p results/variant_db

    get_var_aa_change.py \
      "${annotated_vcf}" \
      "${var_peptides}" \
      "results/variant_db/${sample_id}_var_modified_peptides.fa"
    """
}
