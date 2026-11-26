# Methods

## Data sources

All input data are public. The sample sheets in `data/` list every file with its accession, experiment, md5 checksum and download URL.

| Sheet | Content | Files | Experiments | Cell-line labels |
|---|---|---|---|---|
| `data/samples_hg19_wavelet.csv` | ENCODE replication-timing signal tracks, wavelet-smoothed, hg19 (bigWig): Repli-chip (86 files), UW Repli-seq (16), Repli-seq (6) | 108 | 51 ENCODE experiments and 1 ENCODE UCSC-composite file set (all ENCSR) | 50 |
| `data/samples_rnaseq.csv` | ENCODE RNA-seq gene quantifications (RSEM, TSV) | 104 | 24 ENCODE experiments | 24 |
| `data/samples_annotation.csv` | GENCODE v44 basic annotation lifted to GRCh37 (`gencode.v44lift37.basic.annotation.gtf.gz`, EBI FTP) | 1 | — | — |
| `data/samples_uw_fractions.csv` | QC only: UW Repli-seq cell-cycle fraction signals, six per UW track (UCSC ENCODE archive `wgEncodeUwRepliSeq`, bigWig, hg19) | 96 | 16 UW tracks | 15 |
| `data/samples_drug_fractions.csv` | Drug-treatment Repli-seq, six cell-cycle fraction bigWigs per condition (hg19, 1 kb bins), members of the Zenodo archive of Koundrioukoff et al. 2026 | 24 | 4 conditions | JEFF lymphoblastoid cells |

The GENCODE file has md5 74b40c9bdb697cb085b005555eeef1a7. The drug-treatment archive is `Koundrioukoff_et_al_NAR_2026.zip` from Zenodo record 19093559 (doi:10.5281/zenodo.19093559, CC BY 4.0, md5 c50336d2b1f887bdd80456208713b92c, 3.4 GB); the pipeline downloads it once and extracts the 24 bigWigs, each md5-checked. The raw reads of that study are in GEO under GSE302744. Each cell-line label is the ENCODE biosample name. One label can pool several tracks (for example, replicates or both assays). The 50 labels fall in three groups: cancer (11), non-cancer (29) and stem (10).

## Gene sets

Protein-coding genes came from the GENCODE v44 lift37 basic annotation (20,176 genes; Frankish et al. 2021). Each gene has one transcript: the Ensembl canonical transcript, else the MANE Select transcript, else the longest. The TSS of that transcript is the gene's TSS. The replication-dependent histone genes (66) are the protein-coding genes of the four histone clusters HIST1 (6p22), HIST2 (1q21), HIST3 (1q42) and HIST4 (12p13), selected by HGNC symbol. Replication-independent variants (for example H3-3A, H2AZ1, CENPA) are not included.

## Replication-timing tracks

We used the ENCODE processed replication-timing files as released, with no reprocessing: all 108 human files of output type *wavelet-smoothed signal* (bigWig, hg19) on the ENCODE portal (ENCODE Project Consortium 2020). They come from three assays: Repli-chip (86 files, Gilbert lab; Ryba et al. 2010), UW Repli-seq with six cell-cycle fractions (16 files, Stamatoyannopoulos lab; Thurman et al. 2007; Hansen et al. 2010) and early/late Repli-seq (6 files, Gilbert lab, ENCODE pipeline ENCPL734EDH; Marchal et al. 2018). The labs produced these signals with their own processing, as recorded on the ENCODE portal; ENCODE has no uniform replication-timing pipeline. 

| Assay | Files | Labels | Lab | Processing recorded by ENCODE | Assay reference |
|---|---|---|---|---|---|
| Repli-chip | 86 | 33 | David Gilbert, FSU | Lab processing of the early/late microarray signal; 16 files from the ENCODE2 UCSC composite `wgEncodeFsuRepliChip` | Ryba et al. 2010 |
| Repli-seq (UW, 6 fractions) | 16 | 15 | John Stamatoyannopoulos, UW | ENCODE2 UCSC composite `wgEncodeUwRepliSeq`: weighted average of the cell-cycle fractions, wavelet-smoothed | Hansen et al. 2010 |
| Repli-seq (early/late) | 6 | 6 | David Gilbert, FSU | ENCODE pipeline ENCPL734EDH ("Gilbert Repli-seq Pipeline"), step "Repli-seq coverage step" (bedtools 2.26.0) | Marchal et al. 2018 |

The 6 early/late Repli-seq tracks are the only tracks of 6 cancer labels (A549, Caki2, G401, NCI-H460, SK-N-MC, T47D). Five of these labels have RNA-seq. So the figures need them: without them, section 1 would have 19 labels instead of 24 and section 2 would have 44 instead of 50.

On 2026-11-25, all 108 files were still released, none was superseded, every md5 checksum was unchanged, and the ENCODE search for released wavelet-smoothed signal returned exactly these 108 human files.

## Track QC

The analysis uses the released tracks. One check compares them with our own processing; it does not change the analysis.

We rebuilt each of the 16 UW Repli-seq tracks from its six cell-cycle fraction signals (UCSC ENCODE archive `wgEncodeUwRepliSeq`, `*PctSignalRep1.bigWig`, `data/samples_uw_fractions.csv`) with `scripts/wavelet_rt.py`. This script is a port of `15_wavelet_smooth.py` of the original analysis and follows the ENCODE/UW construction (Thurman et al. 2007; Hansen et al. 2010): per 1 kb bin, percentage-normalise the six fractions, take the weighted average with weights 0.917 (G1b), 0.750 (S1), 0.583 (S2), 0.417 (S3), 0.250 (S4) and 0 (G2), then smooth with a Daubechies-4 discrete wavelet transform at level 7 with all detail coefficients set to zero (PyWavelets). We compared each rebuilt track with each released track by Pearson r over 10 kb bin means (pyBigWig), chromosomes 1-22 and X.

| QC measure | Value |
|---|---|
| GM12878, r of released vs rebuilt track (284,341 bins) | 0.998 |
| r of each released UW track with its own rebuilt track (16 tracks) | 0.996 to 0.998 |
| Tracks whose own rebuilt track has the highest r | 16 of 16 |
| Highest r with another rebuilt track | 0.989 (the other BJ replicate) |

The input fraction signals are UW's own percent-normalised signals, the same inputs as the released tracks. So this check tests only the weighting and smoothing steps, not read alignment or normalisation. The 86 Repli-chip and 6 early/late tracks have no fraction signals and are not checked. The QC figures are in `figures/qc/`.

## Replication-timing score

For each track and each protein-coding gene, we took the mean signal over the TSS ± 5 kb (pyBigWig zoom-level mean). We checked the signal orientation with 11 early housekeeping genes against 6 late anchor genes. No track needed a sign flip. The genome reference of a track is its mean signal in consecutive 10 kb bins over chromosomes 1-22 and X. The replication-timing (RT) score of a gene is 1 minus the mid-rank percentile of its signal among these bins. A score of 0 is the earliest and 1 the latest replication in that track.

## Top 5% expressed genes

For each label with RNA-seq data, we averaged TPM over all its gene-quantification files per Ensembl gene. A protein-coding gene is "top 5% expressed" if its mean TPM is at or above the 95th percentile of that label's protein-coding genes. Four ENCODE files have no TPM column and were not used: three A549 featureCounts files (ENCFF056GXX, ENCFF720ZBH, ENCFF910WIQ) and the only T-cell file (ENCFF603KRU). So 24 of the 50 labels have a top-5% set.

The RNA-seq experiment of each label was chosen once, by the first ENCODE search result with gene-quantification files (polyA plus RNA-seq first, then total RNA-seq); `data/samples_rnaseq.csv` freezes that choice. On 2026-11-25, all 104 files were still released, none was superseded and every md5 checksum was unchanged. A new run of the same search chose a different experiment for 7 labels (HeLa-S3, HepG2, K562, NCI-H460, SK-N-SH, mesenchymal stem cell, T-cell). All 14 experiments involved were released before the original analysis, so the cause is the order of ENCODE search results, not new data. The frozen experiment is the more recent release for 5 of the 7 labels, the same date for 1 and older for 1 (SK-N-SH: 2012-02 vs 2012-06). We kept the frozen choice.

## Drug-treatment Repli-seq

The drug data are from Koundrioukoff et al. 2026: Repli-seq of lymphoblastoid JEFF cells, untreated (NT) or treated with an ATR inhibitor (ATRi), aphidicolin (Aph06) or both (ATRi+Aph06). The authors' Zenodo archive holds 1 kb fraction bigWigs (folder `1kb_adj`). Only four conditions have all six fractions there: NT replicate 3 and replicate 1 of the three treatments. The other samples have only the G1 fraction and are not used. These four tracks are not replicates of each other, so sections 3-5 report no replicate statistics.

For each condition, `scripts/wavelet_rt.py` builds one track with the same method as the UW QC above. The archive names the fractions G1, S1-S4 and G2M; we read them as G1b, S1-S4 and G2. `scripts/drug_gene_rt.py` then gives each protein-coding gene the mean signal over its exact gene body (1 kb bins that overlap [start, end]). The RT score is 1 minus the mid-rank percentile of that signal among the 10 kb genome bins (chromosomes 1-22 and X) of the same track. The output table keeps the 66 RD histone genes. A histone gene's class (H1, H2A, H2B, H3, H4) comes from its HGNC symbol.

Note that the ENCODE sections use the TSS ± 5 kb and the drug sections use the exact gene body, as in the original analysis.

## Statistics and figures

The notebook `analysis.ipynb` makes all figures. Its sections 1-5 are sections 1, 2, 6, 7 and 8 of the original `analysis_slim.ipynb`; the plotting code is unchanged except for input paths and figure saving. Sections 1-2 keep (track, gene) values with an RT score in [0, 1]. For each label, it takes the median RT score over all tracks of that label and all genes of a set. Each point in a figure is one label. Boxes show the median and interquartile range; whiskers extend to 1.5 times the interquartile range; outliers are not drawn separately.

| Figure | Sets | Labels | Reference lines |
|---|---|---|---|
| Section 1 | all protein-coding genes (the histone and top-5% genes included), top 5% expressed genes of that label, the 66 RD histone genes | 24 (labels with a top-5% set); grey lines join the three medians of one label | 0.5 |
| Section 2 | each histone cluster HIST1-HIST4; left panel pooled, right panel by group (cancer, non-cancer, stem) | 50 | 0.5, and the median over labels of the all-gene median (0.23) |

| Section 3 | H1 versus core (H2A, H2B, H3, H4) histone genes, and each class, per drug condition | 1 track per condition | — |
| Section 4 | per-gene RT score change of each treatment against NT (heatmap) | 1 track per condition | — |
| Section 5 | mean change over the three treatments, per histone repeat element (genes with the same cluster and copy number, for example H2AC4 and H3C4) and class | 1 track per condition | — |

Each figure is saved as PDF (editable text) and 300 dpi PNG in `figures/`, with its plotted values as `*_values.csv`.

Sections 2 and 3 draw strip-plot points with random horizontal jitter; the notebook fixes the seed, so the figure is the same on every run.

## Caveat

In these data, the cell-line group is confounded with the assay. All 11 cancer labels have Repli-seq data. Of the stem labels, 9 of 10 have Repli-chip data only; of the non-cancer labels, 20 of 29. The right panel of section 2 compares the groups; that comparison is confounded with the assay and would need an assay-matched design.

## Running the pipeline

`uv sync` creates the pinned Python environment (`.venv`). `nextflow run ../pipeline/main.nf -profile local` runs everything from a new folder. Internet access is needed to encodeproject.org, ftp.ebi.ac.uk, hgdownload.soe.ucsc.edu and zenodo.org.

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

A full run downloads about 7 GB (108 ENCODE tracks, 104 RNA-seq tables, 1 GTF, 96 UW fraction bigWigs, and the 3.4 GB drug archive). It takes about 1 hour on Slurm.

To rerun only the notebook on existing pipeline outputs, from the repository folder:

```bash
RT_DATA=run/results/notebook_inputs RT_FIGDIR=figures .venv/bin/jupyter nbconvert --to notebook --execute --inplace analysis.ipynb
```

| Step | Output (`results/`) | Content |
|---|---|---|
| `FETCH_*` | `cache/` | Downloaded inputs, md5-checked |
| `GENESETS` | `genesets/` | Protein-coding genes (one canonical TSS each) and the 66 RD histone genes |
| `GENE_RT` | `gene_rt/<accession>.parquet` | RT score of every protein-coding gene in one ENCODE track |
| `RNASEQ_TOP5` | `rnaseq_top5.parquet` | Top 5% expressed protein-coding genes per label |
| `EXTRACT_DRUG` | `cache/drug_fractions/` | The 24 drug-treatment fraction bigWigs from the Zenodo archive |
| `WAVELET_DRUG` | `drug/tracks/` | One RT track per drug condition |
| `DRUG_GENE_RT` | `drug/drug_per_gene.csv` | RT score of each histone gene in each drug track |
| `ASSEMBLE` | `notebook_inputs/` | The four tables the notebook reads |
| `NOTEBOOK` | `figures/`, `analysis.executed.ipynb` | The figures, their plotted values and the executed notebook |
| `FETCH_UW`, `REPROCESS_UW` | `qc/reprocessed_uw/` | QC: the 16 UW tracks rebuilt from their six fractions |
| `QC_TRACKS` | `qc/track_qc/` | QC figures and correlation tables |

| Repository path | Content |
|---|---|
| `analysis.ipynb` | The analysis, executed, with its plots |
| `pipeline/` | Nextflow pipeline (`main.nf`, `nextflow.config`, `modules/`) |
| `scripts/` | Python steps called by the pipeline |
| `data/` | Sample sheets: accession or archive path, md5 and URL of every input file |
| `figures/` | Figures (PDF and 300 dpi PNG), plotted values, and QC figures (`figures/qc/`) |
| `docs/pipeline.png`, `docs/pipeline.dot` | Pipeline graph (Graphviz; render with `dot -Tpng -Gdpi=150 docs/pipeline.dot -o docs/pipeline.png`) |
| `pyproject.toml`, `uv.lock` | Pinned Python environment |

## Reproduction check

A full run from the public inputs reproduces the original analysis. The original tables are the ones the original notebook read.

| Output | Rows compared | Result |
|---|---|---|
| Protein-coding gene set | 20,176 | identical |
| RD histone gene set | 66 | identical |
| Gene RT score, 108 ENCODE tracks | 1,959,932 | max abs difference 0.0 |
| Top 5% expressed genes, 24 labels | 23,917 | identical gene sets; TPM max abs difference 0.0 |
| Section 1 plotted values (label x gene-set medians) | 50 labels x 3 sets | max abs difference 1e-16 |
| Section 2 plotted values (label x cluster medians) | 200 | max abs difference 1e-16 |
| Drug gene RT table, 4 conditions | 249 (gene x condition) | rt_score and raw_signal max abs difference 0.0 |
| Sections 3-5 printed tables | 3 sections | identical text to the original notebook |

The 1e-16 differences come from the CSV round trip of the plotted values. No ENCODE track needed a sign flip.

## Software

Nextflow 25.10.2 (Di Tommaso et al. 2017); Python 3.12 with numpy 2.4.4, pandas 3.0.2, pyarrow 20.0.0, pyBigWig 0.3.25, matplotlib 3.10.9, seaborn 0.13.2, PyWavelets 1.10.0 and nbconvert (notebook execution). Exact versions are in `uv.lock`.

## References

- Di Tommaso P, et al. Nextflow enables reproducible computational workflows. Nat Biotechnol 35, 316-319 (2017).
- ENCODE Project Consortium. Expanded encyclopaedias of DNA elements in the human and mouse genomes. Nature 583, 699-710 (2020).
- Frankish A, et al. GENCODE 2021. Nucleic Acids Res 49, D916-D923 (2021).
- Hansen RS, et al. Sequencing newly replicated DNA reveals widespread plasticity in human replication timing. Proc Natl Acad Sci USA 107, 139-144 (2010).
- Koundrioukoff S, et al. ATR and TopBP1 oppose to control dormant origin activity and global replication dynamics, providing a first defense against replication stress. Nucleic Acids Res 54, gkag315 (2026). Data: Zenodo doi:10.5281/zenodo.19093559; GEO GSE302744.
- Marchal C, et al. Genome-wide analysis of replication timing by next-generation sequencing with E/L Repli-seq. Nat Protoc 13, 819-839 (2018).
- Ryba T, et al. Evolutionarily conserved replication timing profiles predict long-range chromatin interactions and distinguish closely related cell types. Genome Res 20, 761-770 (2010).
- Thurman RE, Day N, Noble WS, Stamatoyannopoulos JA. Identification of higher-order functional domains in the human ENCODE regions. Genome Res 17, 917-927 (2007).
