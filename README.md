# CRISPRi_data_analysis
Tool for working with data from CRISPRi experiments. Goal is to improve gold standards validation, reproducibility, and accessiblity of this type of analysis.

> # Please see ***tutorials/tutorial.ipynb*** for a walkthrough of installation, implementation, and full functionality.

#### Inputs for typical workflow
- folder storing raw, paired-end fastq read files from a CRISPRi experiment, with naming convention for replicates and controls
- guide RNA database file for alignment reference
- wash control name
- optional parameters for merging, alignement, and relative abundance analysis
- desired output file

#### Features & Uses
- possible to run everything from a Jupyter notebook (python)
- plotting features: polar, volcano, violin, pathways, individual guide enrichment, kmeans clustering
- generalizable user inputs
- different branches for gold standards comparisons: normalization, fitness scores, differential expression, false discovery rate
- non-targeting (negative control guide) and wash control diagnostic functions
