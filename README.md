# TreeHarmonizer

TreeHarmonizer is a utility that is used to place called variants onto a pre-existing phylogenetic tree, allowing for visualization of variant trajectories and evolutionary progression. TreeHarmonizer was developed with single nucleotide (SNV), structural (SV), and copy number (CNA) variants in mind, allowing for placement of each variant type.

As of this verion 1.0, TreeHarmonizer works as a jupyter notebook, designed for the paper "Long-read sequencing of melanoma subclones reveals multifaceted and parallel tumor progression", by Liu & Goretsky, et al. 

It is being actively developed to serve as a standalone utility for multiple variant calling inputs and trees alike.

### Dependencies

TreeHarmonizer requires a Python 3.6 environment with the following packages -
* pandas
* io
* functools 
* os
* intervaltree (https://pypi.org/project/intervaltree/)
* ete3 (Can be installed with conda - https://etetoolkit.org/download/)
* jupyter

TreeHarmonizer is being updated to work with ete4, which will allow for Python versions > 3.6 when used with jupyter.
