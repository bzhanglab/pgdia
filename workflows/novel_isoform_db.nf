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

  conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

  
  input:
    tuple val(id), path(ids_list), path(stringtie_gtf)
  output:
    tuple val(id), path("${id}_novel_isoforms.gtf")
  script:
    """
    set -euo pipefail
    python3 ${projectDir}/bin/get_transcript.py "$ids_list" "$stringtie_gtf" "${id}_novel_isoforms.gtf"
    
    """
}

process filter_gtf_by_fai {
  tag "${id}"
  input:
    tuple val(id), path(novel_gtf)
    path(fai)
  output:
    tuple val(id), path("${id}_novel_isoforms.filtered.gtf")
  script:
    """
    set -euo pipefail
    awk 'BEGIN{FS="\\t"} FNR==NR{contigs[\$1]=1; next} /^#/ || contigs[\$1] { print }' "$fai" "$novel_gtf" > "${id}_novel_isoforms.filtered.gtf"
    """
}



process transdecoder_longorfs {
  tag "${id}"

  input:
    tuple val(id), path(fasta)
  output:
    tuple val(id), path(fasta), path("*.transdecoder_dir")
  script:
    """
    set -euo pipefail

    rm -rf "*.transdecoder_dir"
    
    TransDecoder.LongOrfs -t ${fasta} -m 30 -O .

    """
}

process transdecoder_predict {
  tag "${id}"

  input:
    tuple val(id), path(fasta), path(td_dir)

  output:
    tuple val(id), path("${id}.pep.fasta")
  
  script:
    """
    set -euo pipefail

    # Predict expects: <basename(fasta)>.transdecoder_dir
    expected_dir="\$(basename "$fasta").transdecoder_dir"

    # Nextflow may stage td_dir under a different name, so link it
    if [[ ! -d "\$expected_dir" ]]; then
      ln -s "$td_dir" "\$expected_dir"
    fi

    fasta=\$(ls -1 *.fasta | head -n 1)
    
    TransDecoder.Predict -t "\$fasta" --retain_long_orfs_length 30 -O .

    # TransDecoder outputs a peptide fasta like: <fasta>.transdecoder.pep
    pep=\$(ls -1 *.transdecoder.pep | head -n 1)
    cp "\$pep" "${id}.pep.fasta"
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
      .map { meta, tmap -> tuple(meta.id, tmap) }
      .join( annotated_gtf.map { meta, gtf -> tuple(meta.id, gtf) }, by: 0, failOnMismatch: true )
      .map { id, tmap, gtf -> tuple(id, tmap, gtf) }

    // 3) Filter .tmap file for novel isoforms
    def filtered_ids_ch = filter_isoforms(tmap_and_gtf_ch)

    // 4) transcript ids list
    def ids_list_ch = extract_ids(filtered_ids_ch)

    // 5) extract novel isoforms gtf
    def novel_gtf_ch = get_novel_transcripts(ids_list_ch)
    // novel_gtf_ch: tuple( id, path("${id}_novel_isoforms.gtf") )

    // 6) filter out records on contigs absent from the genome fasta
    def ch_fai = Channel.value(file(fai_path, checkIfExists: true))
    def filtered_novel_gtf_ch = filter_gtf_by_fai(novel_gtf_ch, ch_fai)

    // 7) gffread to fasta (kept consistent with earlier wiring)
    gffread_in_ch = filtered_novel_gtf_ch.map { id, novel_gtf ->
      tuple([id: id], novel_gtf)
    }

    GFFREAD(
      gffread_in_ch,
      Channel.value(file(fasta_path, checkIfExists: true))
    )

    novel_fasta_ch = GFFREAD.out.gffread_fasta
    novel_fasta_by_id_ch = novel_fasta_ch.map { meta, fa -> tuple(meta.id, fa) }

    // 9) TransDecoder.LongOrfs
    longorfs_ch = transdecoder_longorfs(novel_fasta_by_id_ch)

    // 10) TransDecoder.Predict
    predict_pep_ch = transdecoder_predict(longorfs_ch)

  emit:
    isoform_db = predict_pep_ch
}
