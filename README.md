# Comigration Metagenomics Analysis Pipeline

This repository contains the analysis pipeline for investigating population structure and co-migration patterns in human gut microbiomes using metagenomics data.

The pipeline is organized into several analysis modules that can be run independently. For detailed instructions, refer to the `Snakefile` and `README.md` within each module directory.

**Module Overview:**

1. **SNV Catalog Pipeline** - Generates single nucleotide variant (SNV) catalogs from assembled genomes (MAGs)
2. **Data Processing** - Computes core genome statistics, genome coverage, and site frequency spectra
3. **Moments Analysis** - Fits demographic models and estimates population parameters
4. **Identical Tract Analysis** - Analyzes identity-by-state genomic segments
5. **CP-HMM Analysis** - Infers recombination events and clonal divergence times using closely related genomes

## System Requirements

### Software Dependencies

The analysis modules use a small, standard scientific-Python stack, provided as a conda
environment file ([`environment.yml`](environment.yml)):

```bash
conda env create -f environment.yml
conda activate comig
```

This installs Python 3.8, Snakemake, and the core packages (`pandas`, `numpy`, `scipy`,
`matplotlib`, `seaborn`, `scikit-learn`, `pyarrow`) plus `moments` and `dadi` for the
demographic inference. Figure generation additionally uses R (≥4.0) with `ggplot2`, `dplyr`,
`tidyr`, and `data.table`.

**Additional tools:**

- CP-HMM: For recombination inference (external dependency from [Liu & Good 2024](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3002472)); see [cphmm repo](https://github.com/zhiru-liu/close_pair_hmm). Install it as a package (clone the repo and run `pip install -e .`); it requires Python ≥3.8 (tested on 3.10) and is numba-accelerated.
- FastSimBac: Bacterial genome simulation (optional, for validation analyses); for installation instructions see [FastSimBac bitbucket](https://bitbucket.org/nicofmay/fastsimbac/src/master/)

> **Note — building SNV catalogs from MAGs.** The `snv_catalog_pipeline/` module additionally
> requires Biopython and the bioinformatics tools MMseqs2 (protein clustering), Prodigal (gene
> prediction), and MUMmer (`nucmer`, `delta-filter`, `show-coords`, `show-snps`, `show-diff`;
> whole-genome alignment). These are not needed if you start from the provided SNV catalogs.

### Apple Silicon Compatibility

- The `moments` package does not natively support ARM architecture. Install it using x86 emulation via the `osx-64` channel and refer to the [moments installation guide](https://momentsld.github.io/moments/). This is the only module that requires x86 emulation; the other analyses (including CP-HMM) run natively on Apple Silicon.


## Reproducing the analyses

Each analysis module has its own `README.md` and `Snakefile` describing how to run it. For an
end-to-end guide — starting from the **Supplementary Data** deposited on Zenodo and mapping each
analysis to the figure it produces — see [`REPRODUCE.md`](REPRODUCE.md).

## Citation
If you use this pipeline, please cite the associated publication:

[Citation information to be added upon publication]

## Contact

For questions or issues, please contact:

- Matt Carter - matthewmcarter2 at gmail.com
- Zhiru Liu - zhiru at stanford.edu