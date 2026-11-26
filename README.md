# hlb-rt-paper

Replication timing of the replication-dependent histone genes: 50 ENCODE cell lines (hg19) and a drug-treatment Repli-seq study.

![Pipeline](docs/pipeline.png)

## Run

Requirements: [Nextflow](https://www.nextflow.io) ≥ 24.10, [uv](https://docs.astral.sh/uv/), `curl`, `unzip` and internet access.

```bash
uv sync                                   # pinned Python environment in .venv
mkdir run && cd run
nextflow run ../pipeline/main.nf -profile local
```

The pipeline downloads all inputs from public sources, then executes the notebook.

## Analysis

[`analysis.ipynb`](analysis.ipynb) makes all the analysis plots. To open the notebook interactively, set `RT_DATA=run/results/notebook_inputs` (the pipeline run's `results/notebook_inputs`). The figures are in [`figures/`](figures/).

## Run options

Internet access is needed to encodeproject.org, ftp.ebi.ac.uk, hgdownload.soe.ucsc.edu and zenodo.org.

On a Slurm cluster, run Nextflow itself inside a batch job and use the `slurm` profile with a params file:

```yaml
# params.yaml
slurm_queue: my-partition               # partition
slurm_options: "--exclude=node1,node2"  # any extra sbatch options
```

```bash
nextflow run ../pipeline/main.nf -profile slurm -params-file params.yaml
```

Put `slurm_options` in the params file, not on the command line: Nextflow reads a value that starts with `--` as a new flag. The `slurm` profile runs the downloads (`FETCH_*`) on the node that runs Nextflow, at most 4 at a time. All other steps run as Slurm jobs.

| Parameter | Default | Meaning |
|---|---|---|
| `--outdir` | `results` | Output folder |
| `--cache_dir` | `cache` | Downloaded inputs, reused by later runs |
| `--cell_lines` | `all` | Comma-separated ENCODE labels, for a quick test |
| `--qc` | `true` | Run the UW track QC |
| `--python` | `../.venv/bin/python` | Python of the uv environment |

A full run downloads about 7 GB (108 ENCODE tracks, 104 RNA-seq tables, 1 GTF, 96 UW fraction bigWigs, and the 3.4 GB drug archive). (~1 hour on Slurm)

To rerun only the notebook on existing pipeline outputs, from the repository folder:

```bash
RT_DATA=run/results/notebook_inputs RT_FIGDIR=figures .venv/bin/jupyter nbconvert --to notebook --execute --inplace analysis.ipynb
```

| Step | Output (`results/`) | Content |
|---|---|---|
| `FETCH_*` | `cache/` | Downloaded inputs, md5-checked |
| `GENESETS` | `genesets/` | Protein-coding genes (gene-record span and strand from GENCODE) and the 66 RD histone genes |
| `GENE_RT` | `gene_rt/<track>.parquet` | Gene-body RT score of every protein-coding gene in one track (ENCODE accession or drug condition) |
| `RNASEQ_TOP5` | `rnaseq_top5.parquet` | Top 5% expressed protein-coding genes per label |
| `EXTRACT_DRUG` | `cache/drug_fractions/` | The 24 drug-treatment fraction bigWigs from the Zenodo archive |
| `WAVELET_DRUG` | `drug/tracks/` | One RT track per drug condition |
| `ASSEMBLE` | `notebook_inputs/` | The four tables the notebook reads, including `drug/drug_per_gene.csv` |
| `NOTEBOOK` | `figures/`, `analysis.executed.ipynb` | The figures, their plotted values and the executed notebook |
| `FETCH_UW`, `REPROCESS_UW` | `qc/reprocessed_uw/` | QC: the 16 UW tracks rebuilt from their six fractions |
| `QC_TRACKS` | `qc/track_qc/` | QC figures and correlation tables |

