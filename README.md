[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.10.5-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

# PG-DIA (Zhang Lab)

PG-DIA is a Nextflow DSL2 workflow for building per-sample customized protein databases from RNA-seq and searching matched DIA-MS data against those databases. This Zhang Lab implementation combines RNA variant calling, transcript assembly, novel isoform ORF prediction, protein database assembly, and DIA-NN reporting in one pipeline.

The repository follows an nf-core-style layout, but the root `README.md` is the best high-level guide for the current Zhang Lab workflow.

![PG-DIA workflow overview](docs/images/pgdia_workflow.png)

## What the pipeline does

For each sample, the workflow:

1. Accepts RNA-seq input as FASTQ, BAM, or CRAM together with one DIA raw file path.
2. Runs the RNA variant branch from the bundled `rnavar` workflow to produce duplicate-marked BAMs and annotated VCFs.
3. Converts annotated variants into variant peptide FASTA entries with `pypgatk`, then annotates amino acid changes.
4. Runs StringTie on the RNA-seq alignment output.
5. Runs `gffcompare`, filters novel transcript models, extracts transcript FASTA with `gffread`, and predicts ORFs with TransDecoder.
6. Combines the reference proteome, variant peptides, and novel isoform peptides into per-sample protein databases.
7. Runs DIA-NN against the per-sample combined FASTA.
8. Post-processes the DIA-NN parquet report into reference and novel peptide/protein matrices.

This is a sample-matched workflow: each RNA-seq sample is paired to one DIA raw input and yields its own customized database and DIA-NN output.

## Inputs

### Samplesheet

The full PG-DIA workflow expects one row per sample with:

- one RNA-seq input mode: `fastq_1`/`fastq_2`, or `bam`/`bai`, or `cram`/`crai`
- one DIA input path in `dia_raw`

An example samplesheet for this branch looks like this:

```csv
sample,fastq_1,fastq_2,bam,bai,cram,crai,dia_raw,strandedness
SAMPLE_A,/data/rna/SAMPLE_A_R1.fastq.gz,/data/rna/SAMPLE_A_R2.fastq.gz,,,,,/data/dia/SAMPLE_A.d,forward
SAMPLE_B,,,/data/rna/SAMPLE_B.markdup.bam,/data/rna/SAMPLE_B.markdup.bam.bai,,,/data/dia/SAMPLE_B.RAW,unstranded
SAMPLE_C,,,,,/data/rna/SAMPLE_C.cram,/data/rna/SAMPLE_C.cram.crai,/data/dia/SAMPLE_C.d,reverse
```

Column notes:

| Column | Required | Description |
| --- | --- | --- |
| `sample` | Yes | Sample identifier. Multiple FASTQ lanes for the same sample may reuse the same name. |
| `fastq_1`, `fastq_2` | Conditional | RNA-seq FASTQs. Use these if starting from raw reads. |
| `bam`, `bai` | Conditional | Coordinate-sorted BAM and index if alignment is already done. |
| `cram`, `crai` | Conditional | CRAM and CRAI if starting from CRAM instead of BAM. |
| `dia_raw` | Yes | DIA-MS raw path for the sample, for example a timsTOF `.d` directory or vendor raw file. |
| `strandedness` | Optional | Passed to StringTie. Supported values are `forward`, `reverse`, and `unstranded`. |

Behavior to keep in mind:

- Multiple FASTQ lanes are supported and are concatenated by sample.
- Only one `dia_raw` entry is supported per sample.
- VCF input is still present in the inherited schema, but VCF-only rows are not sufficient for the full PG-DIA workflow because novel isoform discovery still requires RNA-seq alignment data for StringTie.

### Required reference and workflow parameters

At minimum, plan to provide:

- `--input`
- `--outdir`
- `--protein_reference_db`
- either `--genome` with a configured reference bundle, or explicit `--fasta` and `--gtf`
- `--read_length` when aligning from FASTQ

For the RNA variant branch:

- provide `--dbsnp` and/or `--known_indels` if you want base recalibration
- otherwise set `--skip_baserecalibration`

Important implementation detail:

- `--protein_reference_db` is required by the Zhang Lab workflow even though it is not yet surfaced in the top-level schema-generated docs.

## Quick start

Run from the repository root:

```bash
nextflow run . \
  -profile docker \
  --input samplesheet.csv \
  --outdir results \
  --genome all \
  --read_length 151 \
  --protein_reference_db /path/to/reference_proteome.fa \
  --skip_baserecalibration
```

If you have known sites for GATK RNA recalibration, use them instead of skipping BQSR:

```bash
nextflow run . \
  -profile docker \
  --input samplesheet.csv \
  --outdir results \
  --genome all \
  --read_length 151 \
  --protein_reference_db /path/to/reference_proteome.fa \
  --dbsnp /path/to/dbsnp.vcf.gz \
  --dbsnp_tbi /path/to/dbsnp.vcf.gz.tbi \
  --known_indels /path/to/known_indels.vcf.gz \
  --known_indels_tbi /path/to/known_indels.vcf.gz.tbi
```

Useful runtime options:

- `-profile docker`, `-profile singularity`, `-profile conda`, or `-profile mamba`
- `-resume` to continue from a previous run
- `--diann_image` to point to a DIA-NN container image name or tarball
- `--diann_bin` if the DIA-NN executable path inside the image differs from the default
- `--diann_cpus` to control DIA-NN threads

## Reference configuration

The bundled `conf/igenomes.config` defines a Zhang Lab-specific `--genome all` entry pointing to:

- genome FASTA
- annotation GTF
- STAR index
- VEP cache metadata

If you do not want to use that bundle, provide explicit reference files with parameters such as:

- `--fasta`
- `--fasta_fai`
- `--dict`
- `--gtf`
- `--star_index`
- `--dbsnp`
- `--known_indels`

`--read_length` matters for STAR index generation and alignment. Set it to the actual read length of the RNA-seq data, for example `151` for 2x151 bp libraries.

## Main outputs

Published results are written under `--outdir` and typically include:

| Path | Contents |
| --- | --- |
| `reports/` | MultiQC output and annotation summary reports |
| `pipeline_info/` | Nextflow execution metadata, params, and software versions |
| `annotation/<sample>/` | Annotated and decompressed VCF outputs from the RNA variant branch |
| `stringtie/` | Per-sample StringTie GTF assemblies |
| `isoform_db/<sample>/` | Predicted peptide FASTA from novel isoforms |
| `protein_db/<sample>/` | `<sample>_combined_protein_db.fa` and `<sample>_novel_protein_db.fa` |
| `diann_output/<sample>/` | DIA-NN parquet report, matrices parquet, and postprocessed novel/reference TSV matrices |

Key per-sample deliverables are:

- `<sample>_combined_protein_db.fa`
- `<sample>_novel_protein_db.fa`
- `<sample>_report.parquet`
- `<sample>_novel_matrix.tsv`
- `<sample>_ref_matrix.tsv`

## Repository structure

```text
.
├── main.nf
├── nextflow.config
├── conf/
├── workflows/
│   ├── rnavar_mini.nf
│   ├── variant_db.nf
│   ├── novel_isoform_db.nf
│   ├── combine_db.nf
│   └── diann_pipeline.nf
├── workflows/stringtie/
├── bin/
│   ├── assemble_protein_db.py
│   ├── get_transcript.py
│   ├── get_var_aa_change.py
│   └── process_parquet_report.py
├── assets/
└── docs/
```

## Current limitations and branch-specific notes

- The inherited nf-core template docs and test configs are not yet fully aligned with this Zhang Lab implementation. Use this `README.md` as the primary run guide for now.
- The top-level samplesheet schema still exposes some inherited `rnavar` entry modes. In practice, the full PG-DIA workflow needs RNA-seq alignment information plus `dia_raw` for each sample.
- Only one DIA raw path is supported per sample in the current grouping logic.
- DIA-NN runs against per-sample databases, not against one cohort-wide merged database.
- If `--diann_image` points to a tar archive, the helper step loads it with Docker, so that process needs Docker CLI access on the execution host.

## Credits

PG-DIA in this repository was written and adapted by Wenrong Chen in the Zhang Lab, building on nf-core and bundled `rnavar` components.

## Citation

See [CITATIONS.md](CITATIONS.md) for software and workflow citation details.
