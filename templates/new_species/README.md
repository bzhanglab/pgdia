# New species / cohort

Copy this directory for each new dataset. The launcher targets the repository
containing this template. It does not switch branches or fetch updates.

```bash
cp -a templates/new_species /absolute/path/to/my_species_cohort
cd /absolute/path/to/my_species_cohort
```

When the copied directory is outside the repository, set `PGDIA_PIPELINE` to the
pipeline checkout when launching, or edit its default in `run_pgdia.sh`.

1. Edit `species.config`: fill in every marked species, taxonomy, cohort and
   reference field. Use a FASTA, GTF and protein reference for the same species
   and assembly. Use a short cohort ID without spaces. The `new_species` alias
   can stay as-is; if you rename it, change both `genome` and its `genomes` key.
2. Provide an **uncompressed** genome FASTA and GTF, and create the FASTA index
   before launching with `samtools faidx /absolute/path/to/genome.fa`. Set
   `fasta_fai` to the resulting index. The pipeline builds STAR and a local
   snpEff cache from these references; a prebuilt compatible STAR index is optional.
3. Edit `samplesheet.csv`: replace the example rows with your actual data and
   absolute file paths. Set `cohort_expected_runs` in the config to the number
   of unique run IDs (six in the example). Keep this value an integer in config.
4. Set `strandedness` to the actual library protocol: `forward`, `reverse` or
   `unstranded`. This cohort workflow requires one setting across the input runs.
   Review process memory/CPU settings for the run host, particularly the branch's
   180 GB default request for building a snpEff cache.

| Column | How to fill it |
| --- | --- |
| `sample` | Biological sample ID; reuse it for multiple runs of that sample. |
| `run` | Run ID unique across the entire cohort. Repeated lanes of the same run use the same sample/run pair. |
| `tissue` | Explicit tissue name, independent of how samples are named. Spaces become underscores in FASTA headers. |
| `fastq_1`, `fastq_2` | Gzipped FASTQs. For single-end input, leave `fastq_2` empty. |
| `bam`, `bai` | For aligned input, leave both FASTQ columns empty and supply a coordinate-sorted, duplicate-marked BAM and its BAI. |
| `strandedness` | The library's strand protocol; use one consistent setting for this cohort. |

Use one input mode per row. Leave unused columns empty; do not enter `NA` or
`null`. For BAM input the pipeline uses the supplied alignments without running
Picard duplicate marking again. Separate cohorts are appropriate for incompatible
references or strand protocols. The first FASTQ supplies the cohort read length;
use consistent read lengths across FASTQ inputs.

The template skips base recalibration because it assumes no suitable known-sites
resource. If one exists, configure the commented top-level `dbsnp` and `dbsnp_tbi`
paths and set `skip_baserecalibration = false`. Contig names must agree with the
genome. Variant annotation stays enabled. The template keeps variant sequences
seen in at least one run; change `variant_peptide_min_runs` to require more support.

After editing, print the launch command, then run:

```bash
DRY_RUN=1 bash run_pgdia.sh
bash run_pgdia.sh
```

If the templates were copied outside the repository:

```bash
PGDIA_PIPELINE=/absolute/path/to/pgdia DRY_RUN=1 bash run_pgdia.sh
PGDIA_PIPELINE=/absolute/path/to/pgdia bash run_pgdia.sh
```

`DRY_RUN=1` only prints the command; it does not validate reference contents or
run the pipeline. Unfilled placeholders cause the launcher to stop. Nextflow
validates input paths when it starts. The launcher resumes the latest session
in this cohort's run directory, or the session specified by `RESUME_ID`.

Outputs default to `run/results` and work files to `run/work` beside the copied
templates. Override `RUN_ROOT`, `OUTDIR`, `WORK_DIR`, `NEXTFLOW`, `PGDIA_PIPELINE`,
`CUSTOM_CONFIG`, `INPUT_SAMPLESHEET` or `PROFILE` with environment variables.
Changing `RUN_ROOT` does not rewrite reference or FASTQ paths. Docker and the
required images must be available for the default profile. Set `NXF_OFFLINE=false`
if Nextflow needs to fetch its plugins on a new host. Notifications are disabled
in this template.

This branch produces one combined protein database for the whole cohort,
run-level variant provenance, and StringTie outputs. Novel isoforms assembled
across multiple tissues have `Tissue=pooled`. It does not currently run DIA-NN,
so the samplesheet does not require a DIA raw-file column.
