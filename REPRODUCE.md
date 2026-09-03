# Reproducing the analyses

This guide explains how to reproduce the population-genetic analyses and figures using the
**Supplementary Data** deposited on Zenodo (see the Data Availability section of the manuscript
for the most updated DOI) together with the code in this repository.

## Pipeline overview

The analysis is a chain rooted in the SNV catalogs:

```
SNV catalogs ──┬──► pairwise genome comparisons (dRep) ──► strain sharing   (Fig 1c, Fig 2, Ext. Data Fig 6)
               ├──► identical tracts ──► L99 ──► isolation years            (Fig 3)
               ├──► CP-HMM recombination ──► clonal divergence              (Ext. Data Fig 7)
               └──► site-frequency spectra ──► moments demography           (Fig 4)
```

The deposit provides the derived table at **each** stage, so there are two entry points:

- **(a) From the SNV catalogs** — the full path. Download the SNV-catalog archive
  (`Carter_et_al_2025_SNV_catalogs.tar.gz`, on the same Zenodo record) and recompute everything.
- **(b) From a deposited intermediate** — start midway and run only the downstream steps
  (e.g. compute L99 from `max_run_per_pair.csv` without re-deriving the tracts). This avoids the
  large, compute-heavy upstream steps.

## Setup

1. Install the Python/R dependencies (see [`README.md`](README.md)) and, for CP-HMM, the
   `cphmm` package (`pip install -e` a clone of [close_pair_hmm](https://github.com/zhiru-liu/close_pair_hmm)).
2. Download the Supplementary Data deposit and point `config.py` at it (set `data_base_path`
   and `databatch`).
3. The deposit uses publication-friendly filenames, while the code addresses files by the
   internal `{databatch}` naming. Place/rename the deposited files accordingly:

| Deposited file | Code expects (under the configured paths) |
|---|---|
| `01_metadata/mag_metadata.tsv`, `species_metadata.tsv` | `snv_catalog_path/{databatch}_snv_catalog_mag_metadata.tsv`, `..._species_metadata.tsv` |
| `02_pairwise_genome_comparisons/pairwise_genome_comparisons.csv` | `drep_res_path/{databatch}_drep.csv` (raw dRep columns `name/querry/reference/ani/perc_id` — see `utils/pairwise_utils.py::load_drep_res`) |
| `04_identical_tracts_L99/max_run_per_pair.csv` | `run_path/{databatch}_annotated_max_runs.csv` |
| `05_moments_demographic_inference/sfs/*.snps.txt` | `sfs_path/{databatch}_full/` |
| `03_clonal_divergence_recombination/cphmm_recombination_summary.tsv` | `cphmm_res_path/241022_inference_summary_full.tsv` |

## Running each analysis

Each module's `README.md` lists its scripts and the order to run them; the `Snakefile`s encode
the exact commands.

- **SNV catalogs** — `snv_catalog_pipeline/` builds the per-species SNV catalogs from MAGs. Included for completeness.
- **Identical tracts / L99 → Fig 3** — `identical_tract_analysis/`. The `Snakefile` computes
  tracts and `max_run_per_pair.csv` from the catalogs; then `simulate_mutation_accumulation.py`
  and `infer_isolation_years.py` convert between-population L99 into isolation years.
- **CP-HMM recombination → Ext. Data Fig 7** — `cphmm_analysis/`. Full path:
  `prior_prep.py` → `infer_all.py` → `concat_results.py` (needs the catalogs + `cphmm` package).
  Or start from the deposited `cphmm_recombination_summary.tsv` and run `fit_divergence_years.py`
  for the clonal-divergence trend (`genome_pair_clonal_divergence.csv`).
- **Moments demography → Fig 4** — `moments_analysis/`. Run `snakemake -j <N>` over the deposited
  SFS files to produce the demographic parameters (`moments_demographic_parameters.csv`), split
  times, and residuals; `plot_private_SNVs.py` reproduces the private-SNV analysis (Fig 4b).

## Figures

The main-text figure panels are assembled in `figure_generation/figure_generation.Rmd` from the
per-analysis outputs above. The exact per-panel values are additionally provided as the paper's
figure **source data**, so each panel can also be re-plotted directly from those tables.
