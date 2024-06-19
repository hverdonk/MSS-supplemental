# MANIFEST

## tRNA gene count numbers
- Extended_Data_Table_1.xlsx
- from doi: 10.17912/micropub.biology.000560    

## tRNA abundance
- GSE152621_Dmel_anticodon_counts.csv
- quantified by tRNA sequencing in *Drosophila* BG3-c2 cells
- from https://doi.org/10.1016/j.molcel.2021.01.028
    - per this source, "The correlation between gene copy number and tRNA abundance was also lower in Drosophila BG3-c2 cells (adjusted R2 = 0.79)...consistent with differential tRNA gene use in distinct cell types (Dittmar et al., 2006; Ishimura et al., 2014; Kutter et al., 2011; Schmitt et al., 2014) and highlight that mechanisms beyond gene copy number shape metazoan tRNA pools."
- GEO accession number GSE152621

## isodecoder abundance
- GSE152621_Dmel_isodecoder_counts.csv
- quantified by tRNA sequencing in *Drosophila* BG3-c2 cells
- from https://doi.org/10.1016/j.molcel.2021.01.028
- GEO accession number GSE152621
- count transcripts per tRNA gene, rather than per anticodon
- "Single_isodecoder" is because authors cluster isodecoders, align sequences, and then do deconvolution to assess which specific reference a sequence aligns with. "FALSE" means that the given gene was the dumping ground for all sequences they couldn't conclusively assign (and also for its own sequences)
