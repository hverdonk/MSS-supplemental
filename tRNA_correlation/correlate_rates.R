suppressMessages(library(tidyverse))


# original filename 20240321_drosophila_SynREVCodon_parsed.csv
SRC <- read_csv("./Drosophila_SynREVCodon_parsed.csv", col_types = cols())

# Drosophila tRNA isodecoder count numbers from https://doi.org/10.1016/j.molcel.2021.01.028
# GEO Accession number GSE152621
tRNA_abundance <- read_csv("./GSE152621_Dmel_anticodon_counts.csv", col_types = cols())

# my hypotheses for which codons not present in tRNA are served by which isoacceptors
# assume codons share an isoacceptor if they only differ at the wobble base and 
tRNA_predictions <- read_csv("./tRNA-predicted-codons.csv", col_types = cols())

# translation_table <- read_csv("./universal-code.csv")


###############
## Utilities ##
###############
reverse_complement <- function(seqs){
  # input: vector of DNA sequences
  # output: vector of reverse complements to input DNA sequences
  # split each string into its individual characters, reverse them using the rev function, and then collapse them back into strings using sapply and paste
  seqs <- sapply(seqs, function(x) paste(rev(str_split(x, pattern = "")[[1]]), collapse = ""))
  # get complimentary DNA sequence
  seqs <- seqs %>%
    str_replace_all(pattern="A", "X") %>%
    str_replace_all(pattern="T", "A") %>%
    str_replace_all(pattern="X", "T") %>%
    str_replace_all(pattern="G", "Z") %>%
    str_replace_all(pattern="C", "G") %>%
    str_replace_all(pattern="Z", "C")
  seqs
}



###############
## tRNA DATA ##
###############

# normalize read counts by total number of reads sequenced
# tRNA genes are all roughly the same length (76-90 bp) so I chose NOT to normalize
# by transcript length
tRNA_abundance_clean <- tRNA_abundance %>%
  rename("abundance_rep1" = `Dmel_BG3-c2_rep1`, "abundance_rep2" = `Dmel_BG3-c2_rep2`) %>%
  mutate(cpm_rep1 = (abundance_rep1/sum(abundance_rep1))*1000000,
         cpm_rep2 = (abundance_rep2/sum(abundance_rep2))*1000000) %>%
  filter(!str_detect(Anticodon, "mito")) %>%
  filter(!str_detect(Anticodon, "eColiLys")) %>%
  separate_wider_delim(Anticodon, delim='-', names=c("aa", "anticodon"))


# drop raw count data (optional, but makes the data cleaner and easier to work with)
# filter out nonredundant and nonstandard aa
tRNA <- tRNA_abundance_clean %>%
  select(-abundance_rep1, -abundance_rep2)  %>%
  mutate(codon = reverse_complement(anticodon)) %>%
  filter(!str_detect(aa, "iMet")) %>%
  filter(!str_detect(aa, "Met")) %>%
  filter(!str_detect(aa, "SeC")) %>%
  filter(!str_detect(aa, "Trp"))


###############################
## SynREVCodon per-gene fits ##
###############################
# format SynREVCodon per-gene data
medians <- SRC %>%
  summarise(across(everything(), median)) %>%
  pivot_longer(cols = everything(), names_to = "codon_pair", values_to = "median_rate") %>%
  separate_wider_delim(codon_pair, delim="_", names=c("aa", "codons")) %>%
  separate_wider_delim(codons, delim=":", names=c("codon1", "codon2"))



###################
## Add tRNA data ##
###################

# merge tRNA isoacceptor count data for each codon we have info on
# NA codon entries are because not all codons were listed in the abundance data

## CODON 1 ##
SRC_tRNA <- medians %>% # joint_fit_clean %>%
  left_join(tRNA, by=c("codon1"="codon")) %>% 
  rename("anticodon1"=anticodon, 
         "codon1_cpm_rep1" = cpm_rep1, 
         "codon1_cpm_rep2" = cpm_rep2
         ) %>%
  mutate(anticodon1 = reverse_complement(codon1)) %>%
  select(-aa.x, -aa.y)


## CODON 2 ##
SRC_tRNA <- SRC_tRNA %>% 
  left_join(tRNA, by=c("codon2"="codon")) %>%
  rename("anticodon2"=anticodon, 
         "codon2_cpm_rep1" = cpm_rep1, 
         "codon2_cpm_rep2" = cpm_rep2
         ) %>%
  mutate(anticodon2 = reverse_complement(codon2)) %>%
  select(-aa)



###################################
## Add predicted tRNA abundances ##
###################################

# some tRNAs serve more than one codon. I've predicted which tRNAs serve which other codons
# copy tRNA abundance data from "main codons" served by a tRNA
# associate that abundance data with other codons served by the same tRNA
# based on predictions of which codons share an isoacceptor
predictions_subset <- filter(SRC_tRNA, codon1 %in% tRNA_predictions$main_codon) %>%  # find main_codons in codon1 column to extract their abundance data
  distinct(codon1, .keep_all=TRUE) %>%
  select(codon1, codon1_cpm_rep1, codon1_cpm_rep2) %>%
  rename(main_codon="codon1", cpm_rep1 = "codon1_cpm_rep1", cpm_rep2 = "codon1_cpm_rep2")

predictions <- filter(SRC_tRNA, codon2 %in% tRNA_predictions$main_codon) %>%
  distinct(codon2, .keep_all=TRUE) %>%
  select(codon2, codon2_cpm_rep1, codon2_cpm_rep2) %>%
  rename(main_codon="codon2", cpm_rep1 = "codon2_cpm_rep1", cpm_rep2 = "codon2_cpm_rep2") %>%
  bind_rows(predictions_subset) %>%
  distinct(main_codon, .keep_all = TRUE)

tRNA_predictions_data <- left_join(tRNA_predictions, predictions, by="main_codon")


# Merge and update the data with tRNA abundance predictions

## CODON 1 ##
updated_SRC_tRNA <- SRC_tRNA %>%
  replace(is.na(.), 0) %>%
  left_join(tRNA_predictions_data, by = c("codon1" = "other_codon")) %>%
  mutate(codon1_cpm_rep1 = if_else(is.na(cpm_rep1), codon1_cpm_rep1, codon1_cpm_rep1 + cpm_rep1),
         codon1_cpm_rep2 = if_else(is.na(cpm_rep2), codon1_cpm_rep2, codon1_cpm_rep2 + cpm_rep2)
         ) %>%
  select(-main_codon, -cpm_rep1, -cpm_rep2, -mod)


## CODON 2 ##
updated_SRC_tRNA <- updated_SRC_tRNA %>%
  left_join(tRNA_predictions_data, by = c("codon2" = "other_codon")) %>%
  mutate(codon2_cpm_rep1 = if_else(is.na(cpm_rep1), codon2_cpm_rep1, codon2_cpm_rep1 + cpm_rep1),
         codon2_cpm_rep2 = if_else(is.na(cpm_rep2), codon2_cpm_rep2, codon2_cpm_rep2 + cpm_rep2)
  ) %>%
  select(-main_codon, -cpm_rep1, -cpm_rep2, -mod, -aa.x, -aa.y)


# calculate difference in isoacceptor abundance for each susbtitution
updated_SRC_tRNA <- updated_SRC_tRNA %>%
  mutate(cpm_ratio_rep1 = if_else((codon1_cpm_rep1/codon2_cpm_rep1) > 1, true=codon1_cpm_rep1/codon2_cpm_rep1, false=codon2_cpm_rep1/codon1_cpm_rep1)) %>%
  mutate(cpm_ratio_rep2 = if_else((codon1_cpm_rep2/codon2_cpm_rep2) > 1, true=codon1_cpm_rep2/codon2_cpm_rep2, false=codon2_cpm_rep2/codon1_cpm_rep2))


# For substitutions where each codon in a pair is served by a different tRNA with a different abundance
# I would expect to see slower rates than for a pair of codons served by the same tRNA



## tRNA ABUNDANCE CORRELATION WITH JOINT RATE ESTIMATES ##
median_MSS_rates <- updated_SRC_tRNA$median_rate
tRNA_abundance_ratios <- updated_SRC_tRNA$cpm_ratio_rep1

cor.test(median_MSS_rates, tRNA_abundance_ratios, use="pairwise.complete.obs", method="spearman") 






