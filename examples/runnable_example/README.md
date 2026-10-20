# PG-DIA runnable example

This directory supplies the data used by the repository's `test` profile. It is
a small, single-sample example for checking pipeline setup and execution. All
RNA-seq, variant, reference, and DIA-MS content is simulated and must not be
used to evaluate identification sensitivity, quantitative accuracy, or
biological performance.

## What the example exercises

The test run reads a coordinate-sorted RNA-seq BAM, performs the RNA variant
calling and filtering steps, assembles transcripts, builds variant and isoform
protein entries, combines them with the reference proteins, searches the
simulated DIA data with DIA-NN, and post-processes the DIA-NN report.

To keep the example small, variant annotation is not recomputed because a VEP
cache is much larger than the test data. The bundled VEP-style annotated VCF is
therefore used by the protein-database branch. Variant calling and filtering
still run and can be inspected in the Nextflow work directory.

## Included inputs

| File                                               | Purpose                                                                   |
| -------------------------------------------------- | ------------------------------------------------------------------------- |
| `samplesheet.csv`                                  | One-row PG-DIA input sheet linking the RNA-seq and DIA-MS files           |
| `data/simulated_rna.bam` and `.bai`                | Coordinate-sorted RNA-seq alignment and index                             |
| `data/simulated_genome.fasta`, `.fai`, and `.dict` | Small reference genome and indexes                                        |
| `data/simulated_annotation.gtf`                    | Matching gene and transcript annotation                                   |
| `data/simulated_vep_annotated.vcf`                 | One pre-annotated missense variant for database construction              |
| `data/simulated_dia.mzML`                          | 71 MS1 and 284 MS2 centroided scans in four DIA windows                   |
| `data/simulated_dia_ground_truth.tsv`              | Ten peptide sequences and their simulated precursor/retention-time values |
| `data/simulated_proteins.fasta`                    | Reference proteins containing the simulated DIA peptides                  |

## Step 1: install the runtime

The example requires:

- Nextflow 24.10.5 or later
- Java 17 or later
- Docker, with permission to run containers
- the project-specific DIA-NN Docker image used by PG-DIA

Check the first three requirements:

```bash
nextflow -version
java -version
docker version
```

## Step 2: configure DIA-NN

The default profile expects a local Docker image named `diann-2.0.2` and the
executable `/diann-2.0.2/diann-linux` inside that image. Confirm the image and
binary:

```bash
docker image inspect diann-2.0.2 >/dev/null
docker run --rm diann-2.0.2 /diann-2.0.2/diann-linux --help
```

If the image tag or executable differs, set both values before running:

```bash
export PGDIA_DIANN_IMAGE="your-local-diann-image:tag"
export PGDIA_DIANN_TAG="your-local-diann-image:tag"
export PGDIA_DIANN_BIN="/path/inside/container/diann-linux"
```

For an image tar file, set `PGDIA_DIANN_IMAGE` to the tar path and
`PGDIA_DIANN_TAG` to the image tag that should be assigned after loading it.
These variables affect only the test profile and avoid editing workflow files
for machine-specific paths.

## Step 3: inspect and validate the configuration

Run all commands from the repository root. The first command resolves the test
profile without launching tasks. The second validates the sample sheet and
constructs the execution graph without running the full workflow.

```bash
nextflow config -profile test,docker
nextflow run main.nf \
  -profile test,docker \
  --outdir results/test \
  -preview
```

## Step 4: run the example

The complete run can take several minutes and downloads workflow containers if
they are not already cached.

```bash
set -o pipefail
mkdir -p results/test/pipeline_info
nextflow run main.nf \
  -profile test,docker \
  --outdir results/test \
  -with-report results/test/pipeline_info/execution_report.html \
  -with-trace results/test/pipeline_info/execution_trace.txt \
  -with-timeline results/test/pipeline_info/execution_timeline.html \
  2>&1 | tee results/test/pipeline_info/nextflow.log
```

The profile uses Docker, limits each task to four CPUs, 15 GB memory, and one
hour, and assigns two threads to DIA-NN. To resume after an interrupted run,
append `-resume` before the log redirection. Monitor progress in the terminal or
inspect `results/test/pipeline_info/nextflow.log` and the trace file.

## Step 5: check the outputs

The key expected paths are listed in
[`expected_outputs/README.md`](expected_outputs/README.md). A successful run
should include the per-sample combined protein database and DIA-NN report under
`results/test`.

## Regenerate the fixtures

Regeneration is optional because all inputs are included. Python uses only its
standard library; `samtools` is additionally required for BAM and reference
indexes.

```bash
python3 examples/runnable_example/generate_simulated_dia.py
python3 examples/runnable_example/generate_workflow_inputs.py
xmllint --noout examples/runnable_example/data/simulated_dia.mzML
samtools quickcheck examples/runnable_example/data/simulated_rna.bam
```

The generators are deterministic and reproduce the committed fixtures.
