nextflow.enable.dsl=2

 include { GFFCOMPARE } from '../modules/nf-core/gffcompare/main'
 include { GFFREAD } from '../modules/nf-core/gffread/main' 


process FILTER_ISOFORMS {
  tag "${meta.id}"
  input:
    tuple val(meta), path(tmap), path(stringtie_gtf)
  output:
    tuple val(meta), path('filtered_novel_isoforms_ids.txt'), path(stringtie_gtf)
  script:
    """
    set -euo pipefail
    awk '\$3 ~ /^(j|i|m|x|u|o)\$/ && \$8 > 8' $tmap > filtered_novel_isoforms_ids.txt
    """
}

process EXTRACT_IDS {
  tag "${meta.id}"
  input:
    tuple val(meta), path(filtered_ids), path(stringtie_gtf)
  output:
    tuple val(meta), path('ids_list.txt'), path(stringtie_gtf)
  script:
    """
    set -euo pipefail
    cut -f 5 $filtered_ids > ids_list.txt
    """
}

process GET_NOVEL_TRANSCRIPT {
  tag "${meta.id}"

  conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

  
  input:
    tuple val(meta), path(ids_list), path(stringtie_gtf)
  output:
    tuple val(meta), path("${meta.id}_novel_isoforms.gtf")
  script:
    """
    set -euo pipefail
    python3 ${projectDir}/bin/get_transcript.py "$ids_list" "$stringtie_gtf" "${meta.id}_novel_isoforms.gtf"
    
    """
}

process FILTER_GTF_BY_FAI {
  tag "${meta.id}"
  input:
    tuple val(meta), path(novel_gtf)
    path(fai)
  output:
    tuple val(meta), path("${meta.id}_novel_isoforms.filtered.gtf")
  script:
    """
    set -euo pipefail
    awk 'BEGIN{FS="\\t"} FNR==NR{contigs[\$1]=1; next} /^#/ || contigs[\$1] { print }' "$fai" "$novel_gtf" > "${meta.id}_novel_isoforms.filtered.gtf"
    """
}



process TRANSDECODER_LONGORFS_PREDICT {
  tag "${meta.id}"

  input:
    tuple val(meta), path(fasta)
  output:
    tuple val(meta), path("${meta.id}.pep.fasta"), emit: pep_fasta

  script:
    """
    set -euo pipefail
    shopt -s nullglob

    rm -rf ./*.transdecoder_dir
    
    TransDecoder.LongOrfs -t ${fasta} -m 30 -O .

    TransDecoder.Predict -t "${fasta}" --retain_long_orfs_length 30 -O .

    # Copy output to a stable name
    pep=\$(ls -1 *.transdecoder.pep | head -n 1)
    cp "\$pep" "${meta.id}.pep.fasta"

    # Remove intermediate TransDecoder outputs (keeps only ${meta.id}.pep.fasta as declared output)
    rm -rf ./*.transdecoder_dir ./*.transdecoder*

    """
}


workflow GENERATE_NOVEL_ISOFORM_DB {

  take:
    annotated_gtf   // channel: [ val(meta), path(gtf) ]  (gtf = per-sample stringtie gtf)

  main:

    /*
     * Build the 3 inputs expected by GFFCOMPARE:
     *   1) per-sample tuple(meta), path(gtf)
     *   2) singleton tuple(meta2), path(fasta), path(fai)
     *   3) singleton tuple(meta3), path(reference_gtf)
     */
    def fasta_path = params.fasta ?: (params.genomes && params.genome && params.genomes.containsKey(params.genome) ? params.genomes[params.genome].fasta : null)
    def fai_path   = params.fasta_fai ?: (fasta_path ? fasta_path + '.fai' : null)
    def ref_gtf    = params.gtf ?: (params.genomes && params.genome && params.genomes.containsKey(params.genome) ? params.genomes[params.genome].gtf : null)

    if (!fasta_path) error("Missing --fasta (and no genomes[].fasta).")
    if (!ref_gtf)    error("Missing --gtf for gffcompare reference GTF.")

    def ch_genome = Channel.value(tuple([id:'genome'], file(fasta_path, checkIfExists:true), file(fai_path, checkIfExists:true)))
    def ch_refgtf = Channel.value(tuple([id:'refgtf'], file(ref_gtf, checkIfExists:true)))
    
    // 1) Run gffcompare
    GFFCOMPARE(annotated_gtf, ch_genome, ch_refgtf)

    // 2) Pair each sample's .tmap with its original stringtie gtf
    //    gffcompare_results.out.tmap: tuple(meta), path(tmap)
    def tmap_ch = GFFCOMPARE.out.tmap

    tmap_and_gtf_ch = tmap_ch
      .map { meta, tmap -> tuple(meta, tmap) }
      .join( annotated_gtf.map { meta, gtf -> tuple(meta, gtf) }, by: 0, failOnMismatch: true )
      .map { meta, tmap, gtf -> tuple(meta, tmap, gtf) }

    // 3) Filter .tmap file for novel isoforms
    def filtered_ids_ch = FILTER_ISOFORMS(tmap_and_gtf_ch)

    // 4) transcript ids list
    def ids_list_ch = EXTRACT_IDS(filtered_ids_ch)

    // 5) extract novel isoforms gtf
    def novel_gtf_ch = GET_NOVEL_TRANSCRIPT(ids_list_ch)
    // novel_gtf_ch: tuple( id, path("${id}_novel_isoforms.gtf") )

    // 6) filter out records on contigs absent from the genome fasta
    def ch_fai = Channel.value(file(fai_path, checkIfExists: true))
    def filtered_novel_gtf_ch = FILTER_GTF_BY_FAI(novel_gtf_ch, ch_fai)

    // 7) gffread to fasta (kept consistent with earlier wiring)
    gffread_in_ch = filtered_novel_gtf_ch.map { meta, novel_gtf ->
      tuple(meta, novel_gtf)
    }

    GFFREAD(
      gffread_in_ch,
      Channel.value(file(fasta_path, checkIfExists: true))
    )

    novel_fasta_ch = GFFREAD.out.gffread_fasta
    novel_fasta_by_id_ch = novel_fasta_ch.map { meta, fa -> tuple(meta, fa) }

    // 9) TransDecoder.LongOrfs TransDecoder.Predict
    def predict_pep_ch =  TRANSDECODER_LONGORFS_PREDICT(novel_fasta_by_id_ch).pep_fasta

  emit:
    isoform_db = predict_pep_ch
}
