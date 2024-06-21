# MANIFEST

## tRNA analysis code
Code to analyze the correlation between *Drosophila melanogaster* tRNA abundance and SynREVCodon estimated rates as described in Verdonk *et. al*. The script `correlate_rates.R` takes the tRNA transcript counts in `GSE152621_Dmel_anticodon_counts.csv` (GEO accession number GSE152621) and a user-supplied `csv` file of SynREVCodon rates (`Drosophila-SynREVCodon.csv` is provided as an example). The script can be run on the command line as follows:

```

```

Both `GSE152621_Dmel_anticodon_counts.csv` and your csv of SynREVCodon rates should be in the same directory as `correlate_rates.R`

## tRNA abundance data
- stored in `GSE152621_Dmel_anticodon_counts.csv`
- quantified by tRNA sequencing in *Drosophila* BG3-c2 cells and represented as transcript counts (unmodified from source)
- GEO accession number GSE152621

## SynREVCodon rate data
- stored in `Drosophila-SynREVCodon.csv`
