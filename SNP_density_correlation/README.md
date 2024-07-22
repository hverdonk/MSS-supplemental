# Drosophila-polymorphism

This repository archives the *Drosophila melanogaster* polymorphism data and code used in Verdonk *et al.* 2024. The analyses presented in the paper use only SNPs annotated as causing synonymous amino acid changes. 

## mean_SNP_densities.csv
Mean density of single nucleotide polymorphisms (SNPs) calculated for each pair of codons reachable by a single nucleotide substitution in *D. melanogaster*.

## calculte_SNP_densities
Contains the pipeline and intermediate data used to calculate mean SNP densities.

## src
Contains the code used to:
- Obtain the frequency (or density) of each codon pair;
- Correlate with codon change rates;
- Produce the plots.


