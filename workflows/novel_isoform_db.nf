nextflow.enable.dsl=2

 include { GFFCOMPARE } from '../modules/nf-core/gffcompare/main'
 include { GFFREAD } from '../modules/nf-core/gffread/main' 


process filter_isoforms {
  tag "${id}"
  input:
    tuple val(id), path(tmap), path(stringtie_gtf)
  output:
    tuple val(id), path('filtered_novel_isoforms_ids.txt'), path(stringtie_gtf)
  script:
    """
    set -euo pipefail
    awk '\$3 ~ /^(j|i|m|x|u|o)\$/ && \$8 > 8' $tmap > filtered_novel_isoforms_ids.txt
    """
}

process extract_ids {
  tag "${id}"
  input:
    tuple val(id), path(filtered_ids), path(stringtie_gtf)
  output:
    tuple val(id), path('ids_list.txt'), path(stringtie_gtf)
  script:
    """
    set -euo pipefail
    cut -f 5 $filtered_ids > ids_list.txt
    """
}

process get_novel_transcripts {
  tag "${id}"
  input:
    tuple val(id), path(ids_list), path(stringtie_gtf)
  output:
    tuple val(id), path("${id}_novel_isoforms.gtf")
  script:
    """
    set -euo pipefail
    python3 ${params.get_transcript_py ?: "./bin/get_transcript.py"} \
      $ids_list $stringtie_gtf ${id}_novel_isoforms.gtf
    """
}

process transdecoder_longorfs {
  tag "${id}"
  input:
    tuple val(id), path(fasta)
  output:
    tuple val(id), path('novel_transcripts.fasta.transdecoder_dir')
  script:
    """
    set -euo pipefail
    TransDecoder.LongOrfs -t $fasta -m 30 -O .
    mkdir -p novel_transcripts.fasta.transdecoder_dir
    mv novel_transcripts.fasta.transdecoder* novel_transcripts.fasta.transdecoder_dir/
    """
}

process transdecoder_predict {
  tag "${id}"
  input:
    tuple val(id), path(longorfs_dir)
  output:
    path('*')
  script:
    """
    set -euo pipefail
    fasta=\$(ls -1 \$longorfs_dir/../novel_transcripts.fasta | head -n 1)
    TransDecoder.Predict -t \$fasta --retain_long_orfs_length 30 -O .
    """
}

workflow generate_novel_isoform_db {

  take:
    annotated_gtf   // channel: [ val(meta), path(gtf) ]  (gtf = per-sample stringtie gtf)

  main:

    /*
     * Build the 3 inputs expected by GFFCOMPARE:
     *   1) per-sample tuple(meta), path(gtf)
     *   2) singleton tuple(meta2), path(fasta), path(fai)
     *   3) singleton tuple(meta3), path(reference_gtf)
     */
    if (!params.fasta) {
        error("Missing required parameter --fasta for novel isoform generation.")
    }

    ch_gtf = annotated_gtf

    ch_genome = Channel
      .value(
          tuple(
              [id: 'genome'],
              file(params.fasta,      checkIfExists: true),
              file(params.fasta_fai ?: params.fasta + '.fai', checkIfExists: true)
          )
      )

    ch_refgtf = Channel
      .value( tuple([id: 'refgtf'], file(params.gtf, checkIfExists: true)) )

    // 1) Run gffcompare
    gffcompare_results = GFFCOMPARE(ch_gtf, ch_genome, ch_refgtf)

    // 2) Pair each sample's .tmap with its original stringtie gtf
    //    gffcompare_results.out.tmap: tuple(meta), path(tmap)
    tmap_and_gtf_ch = gffcompare_results.out.tmap
      .map { meta, tmap -> tuple(meta.id, tmap) }
      .join( annotated_gtf.map { meta, gtf -> tuple(meta.id, gtf) } )
      .map { id, tmap, gtf -> tuple(id, tmap, gtf) }

    // 3) Filter .tmap file for novel isoforms
    filtered_ids_ch = filter_isoforms(tmap_and_gtf_ch)

    // 4) transcript ids list
    ids_list_ch = extract_ids(filtered_ids_ch)

    // 5) extract novel isoforms gtf
    novel_gtf_ch = get_novel_transcripts(ids_list_ch)

    // 6) gffread to fasta (kept consistent with earlier wiring)
    novel_fasta_ch = novel_gtf_ch
      .map { id, novel_gtf -> tuple(id, novel_gtf, file(params.fasta, checkIfExists: true)) } \
      | GFFREAD

    // 8) TransDecoder.LongOrfs
    longorfs_ch = transdecoder_longorfs(novel_fasta_ch)

    // 9) TransDecoder.Predict
    predict_ch = transdecoder_predict(longorfs_ch)

  emit:
    isoform_db = predict_ch.out
}
