# TreeHarmonizer

TreeHarmonizer is a utility for placing called variants onto a pre-existing phylogenetic tree, enabling visualization of variant trajectories and evolutionary progression. It supports single nucleotide variants (SNVs), structural variants (SVs), and copy number alterations (CNAs), allowing for integrated analysis of multiple variant types.

TreeHarmonizer was originally developed for the paper "Long-read sequencing of single cell-derived melanoma subclones reveals divergent and parallel genomic and epigenomic evolutionary trajectories" by Liu & Goretsky, et al.

**Links:**
- Preprint: [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.08.28.672865v1)
- Sequencing data: [NCBI SRA](https://www.ncbi.nlm.nih.gov/bioproject/1307171)
- Additional project files: [Zenodo](https://doi.org/10.5281/zenodo.16883901)

## Installation

### Via Conda (Recommended)

```bash
git clone git@github.com:KolmogorovLab/TreeHarmonizer.git
cd TreeHarmonizer
conda env create -f th_environment.yml -n tree_harmonizer
conda activate tree_harmonizer
```

### Manual Installation

TreeHarmonizer requires Python 3.12+ with the following packages:

| Package | Version | Notes |
|---------|---------|-------|
| pandas | 2.3.3 | Install via conda-forge |
| ete3 | 3.1.3 | Install via conda from etetoolkit channel |
| intervaltree | 3.1.0 | Install via pip |

TreeHarmonizer was tested with Python 3.6.15 with pandas 1.3.3 for both the notebook and standard versions, and Python 3.12.0 with pandas 2.3.3 for the standard version.

## Usage

```bash
python run_th.py --tree-newick <newick_string_or_file> --reference-species <species> [required-path-arguments] [options]
```

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--tree-newick` | Newick-format string OR path to a file containing the tree structure. Accepts `.nwk`, `.txt`, `.tree`, or `.newick` files. Leaf names must match sample folder names in the data directories. |
| `--reference-species` | Reference genome assembly (REQUIRED). Use a predefined species (`mm10` for mouse chr1-19,X,Y or `grch38` for human chr1-22,X,Y) OR provide a path to a UCSC-format chromosome sizes file for custom genomes. See [Custom Genome Support](#custom-genome-support) below. **TreeHarmonizer is species agnostic. This argument is only used for CNA percentage calculations and chromosome validation.** |

At least one variant data path must be provided:

| Argument | Description |
|----------|-------------|
| `--snv-path` | Path to directory containing per-sample SNV VCF files (e.g., from DeepVariant). Each sample should have its own subdirectory with a VCF file. Expected structure: `{sample}/{sample}.vcf`|
| `--cna-path` | Path to directory containing per-sample CNA output. Supports two directory structures: `{sample}/{sample}.bed` (simplified) or `{sample}/bed_output/{sample}_copynumbers_segments.bed` (legacy Wakhan). See [CNA Input Formats](#cna-input-formats) below. |
| `--sv-path` | Path to multi-sample SV VCF file (e.g., from Severus). <br> **Only Severus files are currently supported.** |

### Custom Genome Support

For non-standard genomes or species other than mouse/human, provide a UCSC-format chromosome sizes file:

```bash
# Example: using a custom chromosome sizes file
python run_th.py \
    --tree-newick ./my_tree.nwk \
    --reference-species /path/to/my_genome.chrom.sizes \
    --snv-path ./snv_data/ \
    --output-path ./output/
```

**File format** (tab-separated, two columns):
```
chr1    248956422
chr2    242193529
chrX    156040895
```

- The `chr` prefix is optional and will be stripped internally
- Comments (lines starting with `#`) are ignored
- Both tab-separated and comma-separated formats are supported

### CNA Input Formats

TreeHarmonizer supports two CNA file formats with auto-detection:

#### Directory Structure

CNA files can be organized in either structure (checked in order):
1. **Simplified**: `{cna_path}/{sample}/{sample}.bed`
2. **Legacy (Wakhan)**: `{cna_path}/{sample}/bed_output/{sample}_copynumbers_segments.bed`

#### File Formats

**Wakhan Format** (7 columns, tab-separated, positional):
```
chr	start	end	coverage	copynumber_state	confidence	svs_breakpoints_ids
1	0	3000000	0.0	0	1.0	[]
1	3000001	195471970	34.28	3	0.9634253	[]
```
- Header row (with or without `#` prefix) is optional and will be auto-detected
- Columns are identified by position, not name

**Generic Format** (4 columns, with header):
```
chrom,start,end,copy_number
chr1,0,3000000,0
chr1,3000001,195471970,3
```
- Header row is required
- Supports both CSV (comma) and TSV (tab) separators
- Column names are flexible (case-insensitive):
  - Chromosome: `chrom`, `chr`, `chromosome`
  - Start: `start`, `begin`, `pos_start`
  - End: `end`, `stop`, `pos_end`
  - Copy number: `copy_number`, `cn`, `copynumber`
- The `chr` prefix is optional and will be stripped internally

## Optional Arguments

### Variant Placement Control

By default, placement is performed on all variant types that are provided. In order to limit placement to select variant types, use one of the following arguments.

| Argument | Default | Description |
|----------|---------|-------------|
| `--no-snv-placement` | False | Skip SNV placement even if path is provided |
| `--no-cna-placement` | False | Skip CNA placement even if path is provided |
| `--no-sv-placement` | False | Skip SV placement even if path is provided |

### Sample Selection

| Argument | Description |
|----------|-------------|
| `--sample-list` | Space-separated list of sample names to include. By default, all samples found in data directories are used. <br>**Note:** *this argument is intended for use if you have samples / folders unrelated to the tree samples or project at hand within the input directories. Not including samples that exist on the tree can lead to undefined behavior and placement.* **It is recommended that if the sample-list given does not represent all samples on the newick string, that the tree supplied be modified to reflect this change.** Example: `--sample-list C1 C2 C3` or `--sample-list C1,C2,C3` |

### Chromosome Selection

| Argument | Description |
|----------|-------------|
| `--chromosomes` | Chromosomes to analyze (space or comma-separated). Example: `--chromosomes 1 2 3 X` or `--chromosomes 1,2,3,X`. Default: all chromosomes for the selected species. |
| `--exclude-chromosomes` | Chromosomes to exclude from analysis. Mutually exclusive with `--chromosomes`. |

**Note:** Chromosome names are normalized internally by stripping any `chr` prefix. You can use either `--chromosomes 1,2,X` or `--chromosomes chr1,chr2,chrX` - both will work identically. Input VCF and BED files with `chr` prefixes are handled automatically.

### Regenotyping Options

Regenotyping uses SV and/or CNA data to rescue SNVs that may have been missed due to overlapping deletions.

| Argument | Default | Description |
|----------|---------|-------------|
| `--disable-regenotyping` | False | Disable the regenotyping step |
| `--regenotype-with-sv-only` | False | Use only SV deletions for regenotyping |
| `--regenotype-with-cna-only` | False | Use only CNA data for regenotyping |
| `--fn-rate` | 0.15 | False negative rate for support threshold calculation. Higher values allow more tolerance for missing variant calls. |

### Output Options

| Argument | Default | Description |
|----------|---------|-------------|
| `--output-path` | `./th_output/` | Base output directory for all results |
| `--write-exclusive-vcfs` | True | Write VCFs containing variants exclusive to each tree node |
| `--no-write-exclusive-vcfs` | | Disable exclusive VCF output |
| `--write-cumulative-vcfs` | True | Write VCFs containing variants at each node plus all descendant nodes |
| `--no-write-cumulative-vcfs` | | Disable cumulative VCF output |

### Miscellaneous

| Argument | Description |
|----------|-------------|
| `--verbose` | Enable detailed per-node variant counts during processing |
| `--version` | Display version information |

## Output

TreeHarmonizer generates the following output structure:

```
{output_path}/
├── snv/
│   ├── placed_snv_variants.tsv      # All successfully placed SNV variants with metadata
│   ├── unplaced_snv_variants.tsv    # SNV variants that failed placement thresholds
│   ├── exclusive/                    # VCFs with variants exclusive to each node
│   │   ├── N1.vcf
│   │   ├── N2.vcf
│   │   └── ...
│   └── cumulative/                   # VCFs with variants at node + all descendants
│       ├── N1.vcf
│       ├── N2.vcf
│       └── ...
├── sv/
│   ├── placed_sv_variants.tsv
│   ├── unplaced_sv_variants.tsv
│   ├── exclusive/
│   └── cumulative/
└── cna/
    ├── node_amplification_percentages.tsv    # Amplification % per tree node (by chromosome + total)
    ├── node_loss_percentages.tsv             # Loss % per tree node (by chromosome + total)
    ├── subline_amplification_percentages.tsv # Amplification % per sample (by chromosome + total)
    ├── subline_loss_percentages.tsv          # Loss % per sample (by chromosome + total)
    └── average_percentages.tsv               # Average amplification/loss % across all samples
```

**Exclusive VCFs**: Contain variants that arose specifically at each tree node (branching point).

**Cumulative VCFs**: Contain all variants present in each lineage (node variants plus all descendant node variants).

## Quick Run with Example Data

An example dataset is provided in the `example_data/` directory, containing chromosome 1 data from 23 mouse melanoma subclones.

### Input Files

- `example_data/snv/` - Per-sample SNV VCF files (DeepVariant output)
- `example_data/cna/` - Per-sample CNA segment files (Wakhan output)
- `example_data/sv/severus_chr1.vcf` - Merged SV calls (Severus output)
- `example_data/original_tree.nwk` - Phylogenetic tree in Newick format

### Running the Example

```bash
# Activate the environment
conda activate th_env

# Run with all variant types (using tree file path), regenotype only using SVs.
python run_th.py \
    --tree-newick ./example_data/original_tree.nwk \
    --reference-species mm10 \
    --chromosomes 1 \
    --snv-path ./example_data/snv/ \
    --cna-path ./example_data/cna/ \
    --sv-path ./example_data/sv/severus_chr1.vcf \
    --regenotype-with-sv-only \
    --output-path ./th_output_example/

# Run SNV placement only without regenotyping
python run_th.py \
    --tree-newick ./example_data/original_tree.nwk \
    --reference-species mm10 \
    --chromosomes 1 \
    --snv-path ./example_data/snv/ \
    --disable-regenotyping \
    --output-path ./th_output_example/
```

### Expected Output

The example run will produce:
- `th_output_example/snv/placed_snv_variants.tsv` - Placed SNV variants with tree node assignments
- `th_output_example/snv/unplaced_snv_variants.tsv` - SNV variants that did not meet placement criteria
- `th_output_example/snv/exclusive/` - Per-node exclusive VCF files
- `th_output_example/snv/cumulative/` - Per-node cumulative VCF files
- **Similar output for SV per the general output structure**
- **Similar output for CNAs per the general output structure**

### Expected Runtime

- Example dataset (chr1 only): ~2 minutes
- Full mouse melanoma dataset (all chromosomes): ~10 minutes

Runtimes may vary depending on system specifications.

### Legacy notebook version (deprecated)

#### Installation Via Conda (Recommended)

```bash
git clone git@github.com:KolmogorovLab/TreeHarmonizer.git
cd TreeHarmonizer
conda env create -f notebook_environment.yml -n tree_harmonizer_nb
conda activate tree_harmonizer_nb
```

#### Running

* All notebooks may be run as is, input paths are pre-populated. Input and output for this version is significantly more limited than the standard version and is considered deprecated.

## In Development

- Additional input support for other SV callers.
- CNA visualization improvements.
- Tree visualizations.

## Contact

For bug reports, feature requests, or questions:
- Submit an issue: [GitHub Issues](https://github.com/KolmogorovLab/TreeHarmonizer/issues)
- Primary developer: [anton.goretsky@nih.gov](mailto:anton.goretsky@nih.gov)
