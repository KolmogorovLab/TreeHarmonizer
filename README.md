# TreeHarmonizer

TreeHarmonizer is a utility that is used to place called variants onto a pre-existing phylogenetic tree, allowing for visualization of variant trajectories and evolutionary progression. TreeHarmonizer was developed with single nucleotide (SNV), structural (SV), and copy number (CNA) variants in mind, allowing for placement of each variant type.

As of this verion 1.0, TreeHarmonizer works as a jupyter notebook, designed for the paper "Long-read sequencing of single cell-derived melanoma subclones reveals divergent and parallel genomic and epigenomic evolutionary trajectories", by Liu & Goretsky, et al. 

It is being actively developed to serve as a standalone utility for multiple variant calling inputs and trees alike.

## Installation (via Conda environment creation)

It is easiest to install TreeHarmonizer via the provided conda environment yaml files.
For the notebook version, use `notebook_environment.yml`

```
git clone git@github.com:KolmogorovLab/TreeHarmonizer.git
cd TreeHarmonizer
conda env create -f notebook_environment.yml -n tree_harmonizer_nb
conda activate tree_harmonizer_nb
```

#### Dependencies (for manual installation)

If you would like to create the environment manually, TreeHarmonizer requires a Python 3.6.15 environment with the following packages -
* pandas (v1.1.5)
* intervaltree (v3.1.0) (https://pypi.org/project/intervaltree/)
* ete3 (v3.1.3) (Can be installed with conda - https://etetoolkit.org/download/)
* ete_toolchain (v3.0.0)
* jupyter

Versions listed are those that were tested to be stable and work together.

## In development
* Package version of TreeHarmonizer
* TreeHarmonizer is being updated to work with ete4, which will allow for Python versions > 3.6 when used with jupyter.

##### Contact

For bug reporting, help, and advising, please submit an [[issue](https://github.com/KolmogorovLab/TreeHarmonizer/issues)]. You may also contact the primary developer at [[anton.goretsky@nih.gov](mailto:anton.goretsky@nih.gov)].