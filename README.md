Meltos: Multi-Sample Tumor Phylogeny Reconstruction for Structural Variants
-----------------------

### About
Meltos is a novel approach to estimate the variant allele frequency of somatic SVs from whole genome sequencing (WGS) signals and using these VAFs to identify the phylogenetic relationship between SSVs in a multi-sample cell lineage tree using an SNV-only tree as a basis. Our probabilistic framework allows us to assess multiple types of signals taken from the data simultaneously and more accurately calculate the VAF of SV events. Following the maximum parsimony principle, Meltos uses a novel combinatorial algorithm to assign SV clusters on the branches of the given lineage tree, while modestly augmenting the tree topology if needed.

###Required Files

- SSNV Tree file produced by LICHEE (https://github.com/viq854/lichee)
- Tab delimited file of SSNVs used to build tree
- Configuration file produced by Meltos input preparation script


