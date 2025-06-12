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

### Software Dependencies and Operating Systems

**Core Requirements:**

- **Python**: 3.8 or higher (tested on 3.8.20)
- **R**: 4.0 or higher
- **Snakemake**: 6.0 or higher for workflow management

**Python Packages:**

- `pandas`
- `numpy`
- `matplotlib`
- `seaborn`
- `scipy`
- `moments` - [Population genetic modeling and demographic inference](https://github.com/MomentsLD/moments)
- `dadi` - [Diffusion approximation for demographic inference](https://dadi.readthedocs.io/en/latest/)
- `Bio` (Biopython)

**R Packages:**

- `ggplot2`
- `dplyr`
- `tidyr`
- `data.table`

**Bioinformatics Tools:**

- MMseqs2: Protein clustering and sequence searching
- Prodigal: Gene prediction in prokaryotic genomes  
- MUMmer (nucmer, delta-filter, show-coords, show-snps, show-diff): Whole genome alignment

**Additional Tools:**

- CP-HMM: For recombination inference (external dependency from [Liu & Good 2024](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3002472)); see [cphmm repo](https://github.com/zhiru-liu/close_pair_hmm). This package requires python 3.8 or *lower* for compatibility with the pipeline. Add to Python path as indicated in `cphmm_analysis/infer_all.py`.
- FastSimBac: Bacterial genome simulation (optional, for validation analyses); for installation instructions see [FastSimBac bitbucket](https://bitbucket.org/nicofmay/fastsimbac/src/master/)

**Typical Install Time:**
None of the dependencies should take more than 30 minutes to install, assuming a good internet connection and no major issues with package repositories.

### Versions Tested On

- **Operating Systems**: 
  - macOS (primary development environment)
  - Linux (HPC cluster environments)
- **Python**: 3.8.20
- **Hardware**: Standard desktop/laptop computers with 8+ GB RAM recommended

### Apple Silicon Compatibility

**Important Note for Apple Silicon (M1/M2/M3) Users:**

- The `moments` package does not natively support ARM architecture. Install it using x86 emulation via the `osx-64` channel and refer to the [moments installation guide](https://momentsld.github.io/moments/).
- Performance is significantly degraded when running CP-HMM under x86 emulation. If you plan to perform CP-HMM analysis, we recommend setting up a separate ARM-native environment for this module while using the x86 environment for other analyses.


## Instructions for Use

#### 1. Prepare Input Data
Your SNV catalog data should be organized as follows:
```
data/
├── [batch_name]_snv_catalog_mag_metadata.tsv
├── [batch_name]_snv_catalog_species_metadata.tsv
└── [batch_name]/
    ├── species1/
    ├── species2/
    └── ...
```

#### 2. Configure Analysis Parameters
Edit `config.py` to set:

- Data paths (`data_base_path`, `snv_catalog_path`, etc.)
- Mutation rates and generation times (for rescaling inferred parameters in moments analysis)

#### 3. Run Analysis Modules
See individual module directories for specific configuration options and instructions.

#### 4. Generate Figures
```bash
cd figure_generation
# Open and run figure_generation.Rmd in RStudio or via command line
Rscript -e "rmarkdown::render('figure_generation.Rmd')"
```

## Citation
If you use this pipeline, please cite the associated publication:

[Citation information to be added upon publication]

## Contact

For questions or issues, please contact:

- Matt Carter - matthewmcarter2 at gmail.com
- Zhiru Liu - zhiru at stanford.edu