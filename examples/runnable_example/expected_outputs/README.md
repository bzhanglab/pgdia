# Expected outputs

With `--outdir results/test`, the important output locations are:

| Path | Expected content |
| --- | --- |
| `results/test/stringtie/SIM_DIA.stringtie_output.gtf` | StringTie transcript assembly for `SIM_DIA` |
| `results/test/isoform_db/SIM_DIA/` | TransDecoder peptide FASTA and associated novel-isoform outputs |
| `results/test/var_db/SIM_DIA/SIM_DIA_var_modified_peptides.fa` | Variant-derived protein records with amino-acid-change annotation |
| `results/test/protein_db/SIM_DIA/SIM_DIA_combined_protein_db.fa` | Combined reference, variant, and isoform search database |
| `results/test/protein_db/SIM_DIA/SIM_DIA_novel_protein_db.fa` | Non-reference portion of the search database |
| `results/test/diann_output/SIM_DIA/SIM_DIA_report.parquet` | DIA-NN precursor-level report |
| `results/test/diann_output/SIM_DIA/SIM_DIA_matrices.parquet` | DIA-NN quantitative matrices, when emitted by the installed DIA-NN version |
| `results/test/diann_output/SIM_DIA/SIM_DIA_report_xic/` | Transition-level XIC parquet files, when emitted by the installed DIA-NN version |
| `results/test/diann_output/SIM_DIA/SIM_DIA_novel_matrix.tsv` | Postprocessed non-reference peptide/protein matrix |
| `results/test/diann_output/SIM_DIA/SIM_DIA_ref_matrix.tsv` | Postprocessed reference peptide/protein matrix |
| `results/test/pipeline_info/` | Execution report, trace, timeline, parameters, and software-version metadata |

The simulated ground truth describes signals inserted into the mzML file. Exact
DIA-NN identifications can vary with the installed DIA-NN version and settings;
the fixture is intended to test execution and data flow, not identification
performance.
