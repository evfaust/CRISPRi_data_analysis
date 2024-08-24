# CRISPRi_data_analysis
(Work in progress / jumping off point) tool for working with data from CRISPRi experiments. Goal is to improve gold standards validation, reproducibility, and accessiblity of this type of analysis. Please see TODO's in the tutorial (and create issues from these on GitHub)

> # Please see ***tutorial.ipynb*** for a walkthrough of installation, implementation, and full functionality. IGNORE the warnings - that was just a result of a highly similar test dataset - should be updated with a better test dataset upon further development. See also the plotting .ipynb files for help. These tutorials implement functions defined in .py files (see definitions for full documentation of parameters and outputs - code is modular).

#### Inputs for typical workflow
- folder storing raw, paired-end fastq read files from a CRISPRi experiment, with naming convention for replicates and controls
- guide RNA database file for alignment reference
- wash control name
- optional parameters for merging, alignement, and relative abundance analysis
- desired output file

#### Features & Uses
- possible to run everything from a Jupyter notebook (python)
- plotting features: polar, bar, volcano, top enrichment/depletion, individual guide enrichment along a gene
- generalizable user inputs
- different normalization options for gold standards comparison, and directions for fitness scores, differential abundance, false discovery rate
