nextflow.enable.dsl=2

 include { GFFCOMPARE } from '../modules/nf-core/gffcompare/main'
 include { GFFREAD } from '../modules/nf-core/gffread/main' 


workflow generate_novel_isoform_db {
    take:
    annotated_gtf  // channel: [ val(meta), path(gtf) ]

    main:

    gffcompare_results = GFFCOMPARE(
        annotated_gtf.map { meta, gtf -> gtf },
        Channel.value(file(params.gtf)),
        annotated_gtf.map { meta, gtf -> meta.id }
    )

    // 2. Extract .tmap and pass sample id and stringtie GTF
    tmap_and_gtf_ch = gffcompare_results.out.combine(annotated_gtf).map { files, pair ->
        def tmap = files.find { it.name.endsWith('.tmap') }
        def meta = pair[0]
        def stringtie_gtf = pair[1]
        tuple(meta.id, tmap, stringtie_gtf)
    }


    // 3. Filter .tmap file for novel isoforms (awk step)
    filtered_ids_ch = filter_isoforms(tmap_and_gtf_ch)

    // 4. Cut -f 5 to get list of transcript ids
    ids_list_ch = extract_ids(filtered_ids_ch)

    // 5. Extract GTF entries by id
    novel_gtf_ch = get_novel_transcripts(ids_list_ch)

    // 6. GTF to fasta with gffread
    novel_fasta_ch = novel_gtf_ch.map { novel_gtf ->
        tuple(
            novel_gtf,
            file(params.reference_fasta),
            "novel_transcripts.fasta"
        )
    } | GFFREAD

    // 7. Remove old transdecoder results (can be a process or just ignore/overwrite)
    // Skipping rm for reproducibility; just overwrite in output dir

    // 8. TransDecoder.LongOrfs
    longorfs_ch = transdecoder_longorfs(novel_fasta_ch)

    // 9. TransDecoder.Predict
    predict_ch = transdecoder_predict(longorfs_ch)

    emit:
    out = predict_ch.out

}




process filter_isoforms {
    tag "${id}"
    input:
      tuple val(id), path(tmap), path(stringtie_gtf)
    output:
      tuple val(id), path('filtered_novel_isoforms_ids.txt'), path(stringtie_gtf)
    script:
      """
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
      python3 ./src/get_transcript.py $ids_list $stringtie_gtf ${id}_novel_isoforms.gtf
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
    fasta=\$(ls \$longorfs_dir/../novel_transcripts.fasta)
    TransDecoder.Predict -t \$fasta --retain_long_orfs_length 30 -O .
    """
}