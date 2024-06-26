# MANIFEST

## tRNA analysis code
Code to analyze the correlation between *Drosophila melanogaster* tRNA abundance and SynREVCodon estimated rates as described in Verdonk *et. al*. The script `correlate_rates.R` takes the tRNA transcript counts in `GSE152621_Dmel_anticodon_counts.csv` (GEO accession number GSE152621) and a csv file of SynREVCodon rates estimated for 12,062 *Drosophila* genes (`Drosophila_SynREVCodon_parsed.csv`). The script can be run on the command line as follows:

```
Rscript correlate_rates.R
```

`correlate_rates.R` should be run in the same directory as  `GSE152621_Dmel_anticodon_counts.csv`, `Drosophila_SynREVCodon_parsed.csv`, and `tRNA-predicted-codons.csv`. You should see the following output:

```
	Spearman's rank correlation rho

data:  median_MSS_rates and tRNA_abundance_ratios
S = 71640, p-value = 0.0002869
alternative hypothesis: true rho is not equal to 0
sample estimates:
       rho 
-0.4294888 

Warning message:
In cor.test.default(median_MSS_rates, tRNA_abundance_ratios, use = "pairwise.complete.obs",  :
  Cannot compute exact p-value with ties
```

## tRNA abundance data
- stored in `GSE152621_Dmel_anticodon_counts.csv`
- quantified by tRNA sequencing in *Drosophila* BG3-c2 cells and represented as transcript counts
- no modifications, identical to data stored in GEO
- GEO accession number GSE152621

## SynREVCodon rate data
- stored in `Drosophila_SynREVCodon_parsed.csv`
- each column is a single estimated substitution rate between a pair of synonymous codons, formatted as [amino acid single letter code]_CODON1:CODON2
- each row contains the rate estimates for a single gene

## tRNA predictions
- Some tRNAs serve more than one codon. We must adjust the pool of tRNAs available to each codon accordingly.
- Predictions of which tRNAs serve which other codons are based on common tRNA modifications, and are more fully described in Verdonk *et. al*
- predictions are stored in `tRNA-predicted-codons.csv`
