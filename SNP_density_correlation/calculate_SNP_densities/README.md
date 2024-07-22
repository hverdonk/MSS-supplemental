If you wish to recalculate SNP densities yourself, you can do the following:
1. Use the `remake_vcf` and `anotate_vcf` pipelines in the [Dmel_data](https://github.com/vitorpavinato/dmel_data) repository to generate a tsv of *D. melanogaster* SNPs from the original *D. melanogaster* data from the Zambia (DPGP3) population
2. Use `subset_drosophila_data.ipynb` to filter for synonymous and non-synonymous SNPs and associated gene information for the four major *D. melanogaster* chromosomes. This will generate the same files contained in `data`.
3. Calculate the SNP density for all exonic SNPs for the four major *D. melanogaster* chromosomes.

## data
Contains a .tsv file with only exonic SNPs for the four major D. melanogaster chromosomes.

## `subset_drosophila_data.ipynb`
A Jupyter notebook used to filter the tsv files in `data` for only synonymous and non-synonymous SNPs and associated gene information.
