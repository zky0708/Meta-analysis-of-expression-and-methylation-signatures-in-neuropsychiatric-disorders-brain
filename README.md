# Meta-analysis-of-expression-and-methylation-signatures-in-neuropsychiatric-disorders-brains

This depository contains the consensus probe-selection attractor finding algorithm we used for identifying the signatures 
discussed in Zhu et al., "Meta-analysis of expression and methylation signatures indicates a stress-related 
epigenetic mechanism in multiple neuropsychiatric disorders".

## 1. Required packages
To run the algorithm, the following two packages are required.

- "cafr" -- [GitHub](https://github.com/weiyi-bitw/cafr)
- "limma" -- [Bioconductor](https://bioconductor.org/packages/release/bioc/html/limma.html)

## 2. Data files
- ```findMultipleConsPSAttractors.R``` is the main function we need to call for running the algorithm.
- ```dataList_exprs_PFC.rda``` is the discovery expression data cohort we used in the meta-analysis,
and its form is a R "list" of which each item is a expression data matrix from an individual study.
- ```map_exprs.rda``` contains a vector mapping each probe of Affymetrix Human Genome U133 array to a gene symbol name
(note that multiple probes might be mapped the same gene).
- ```sample_seeds.csv``` contains probes which can be used as seeds to generate the three co-expression signatures in the paper, 
though you can try any other probe as seed to run the algorithm as long as the probe is among the row names 
of the expression matrices.
- ```sample_script.R``` contains a sample code for running the algorithm.
