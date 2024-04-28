library(tidyverse)
library(readxl)

# Drosophila tRNA gene count numbers from doi: 10.17912/micropub.biology.000560
tRNA_genecounts <- read_excel("~/Dropbox/Temple/pond-lab-research/Drosophila-tRNA/Extended_Data_Table_1.xlsx")
# Drosophila tRNA isodecoder count numbers from https://doi.org/10.1016/j.molcel.2021.01.028
# GEO Accession number GSE152621
tRNA_abundance <- read_csv("~/Dropbox/Temple/pond-lab-research/Drosophila-tRNA/GSE152621_Dmel_anticodon_counts.csv")
SRC <- read_csv("~/Dropbox/Temple/pond-lab-research/Drosophila-tRNA/20240321_drosophila_SynREVCodon_parsed.csv")
joint_fit <- read_csv("https://raw.githubusercontent.com/veg/MSS-results/master/empirical-validation/drosophila/20240322-droso-500RandomAlignments-SynREVCodonjointfit-0.csv")
# my hypotheses for which codons not present in tRNA are served by which isoacceptors
# assume codons share an isoacceptor if they only differ at the wobble base and 
#     - there is no other isoacceptor that only differs at the wobble base.
#         - Tyr, His, Asn, Asp get their wobble bases modified to Queosine (https://www-ncbi-nlm-nih-gov.libproxy.temple.edu/pmc/articles/PMC8897278/)
#     - one anticodon begins with A (which is modified to Inosine, a Guanine analog, for nearly all mature tRNAs) and the other anticodon begins with G
#     - Crick's wobble hypothesis predicts that...
#         - Gly: tRNA recognizing GGC also serves GGU, and the tRNA recognizing GGA also serves GGG
#         - Arg: the tRNA recognizing CGU also serves CGC (confirmed by inosine modification) and the tRNA recognizing CGA also serves CGG
# If a twofold redundant aa only has one isoacceptor listed, assume it translates both codons
tRNA_predictions <- tibble(
  aa = c(
    "Ala", "Ala", "Arg", "Arg", 
    "Asn", "Gly", "Gly", "Asp", 
    "Cys", "His", "Ile", "Ile",
    "Phe", "Pro", "Pro", "Ser", 
    "Ser", "Ser", "Tyr", "Val", 
    "Val", "Leu", "Leu", "Thr", 
    "Thr"
  ),
  main_codon = c(  # the main codon recognized by a tRNA
    "GCT", "GCT", "CGT", "CGA", 
    "AAC", "GGC", "GGA",
    "GAC", "TGC", "CAC", 
    "ATT", "ATT",
    "TTC", 
    "CCT", "CCT", "TCT", "TCT", "AGC", 
    "TAC", "GTT", "GTT", "CTT", "CTT", 
    "ACT", "ACT"
  ),
  other_codon = c(  # other codons served by the same tRNA
    "GCC", "GCA",
    "CGC", "CGG",
    "AAT", "GGT", "GGG",
    "GAT", "TGT", "CAT", 
    "ATC", "ATA", 
    "TTT", 
    "CCC", "CCA", "TCC", "TCA", "AGT", 
    'TAT', "GTC", "GTA", "CTC", "CTA",
    "ACC", "ACA"
  ),
  mod = c(  # any modification (real or predicted) that allows a tRNA to recognize the other codon
    "Inosine", "Inosine", 
    "Inosine", "Crick wobble",
    "Queuosine", "Crick wobble", "Crick wobble",
    "Queuosine", "twofold", "Queuosine", 
    "Inosine", "Inosine",
    "twofold", 
    "Inosine", "Inosine", "Inosine", "Inosine", "twofold",
    "Queuosine", "Inosine", "Inosine", "Inosine", "Inosine", 
    "Inosine", "Inosine"
  )

)

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


# Create a tibble with amino acid codes
amino_acids <- tibble(
  Three_Letter_Code = c("Ala", "Arg", "Asn", "Asp", "Cys",
                        "Gln", "Glu", "Gly", "His", "Ile",
                        "Leu", "Lys", "Met", "Phe", "Pro",
                        "Ser", "Thr", "Trp", "Tyr", "Val"),
  One_Letter_Code = c("A", "R", "N", "D", "C",
                      "Q", "E", "G", "H", "I",
                      "L", "K", "M", "F", "P",
                      "S", "T", "W", "Y", "V")
)

###############
## tRNA DATA ##
###############

# format tRNA gene count data
tRNA_genecounts_clean <- tRNA_genecounts %>% 
  select(`FlyBase gene symbol`) %>%
  separate_wider_delim(`FlyBase gene symbol`, delim="-", names=c("aa", "anticodon", NA, "pseudo")) %>%
  separate_wider_delim(aa, delim=":", names=c(NA, "aa")) %>%
  filter(str_detect(pseudo, pattern="Psi", negate=TRUE)) %>%
  select(-pseudo) %>%
  mutate(codon = reverse_complement(anticodon)) %>%
  # count number of tRNA genes per anticodon
  group_by(aa, anticodon, codon) %>%
  summarise(gene_count = n(), .groups = "drop")


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
tRNA_abundance_clean <- tRNA_abundance_clean %>%
  select(-abundance_rep1, -abundance_rep2)

tRNA <- full_join(tRNA_genecounts_clean, tRNA_abundance_clean, by=c("aa", "anticodon"))

###########################
## SynREVCodon joint fit ##
###########################
# format SynREVCodon joint fit data
joint_fit_clean <- joint_fit %>% 
  mutate(rate=str_remove(rate, "SynREVCodon_")) %>%
  filter(!(rate %in% c("AIC", "omega", "non-synonymous/synonymous rate ratio"))) %>%
  mutate(rate=str_remove(rate, "synonymous rate between codon classes ")) %>%
  separate_wider_delim(rate, delim=" and ", names=c("codon1", "codon2")) %>%
  rename(jointfit = "MLE")


###############################
## SynREVCodon per-gene fits ##
###############################
# format SynREVCodon per-gene data
means <- SRC %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols = everything(), names_to = "codon_pair", values_to = "mean_rate") %>%
  separate_wider_delim(codon_pair, delim="_", names=c("aa", "codons")) %>%
  separate_wider_delim(codons, delim=":", names=c("codon1", "codon2")) %>%
  full_join(amino_acids, by=c("aa"="One_Letter_Code"))
medians <- SRC %>%
  summarise(across(everything(), median)) %>%
  pivot_longer(cols = everything(), names_to = "codon_pair", values_to = "median_rate") %>%
  separate_wider_delim(codon_pair, delim="_", names=c("aa", "codons")) %>%
  separate_wider_delim(codons, delim=":", names=c("codon1", "codon2")) %>%
  full_join(amino_acids, by=c("aa"="One_Letter_Code"))


#########################################
## Combine SRC and Joint Fit estimates ##
#########################################
SRC_rollup <- full_join(means, medians) %>% 
  filter(aa != "M", aa != "W") %>%
  full_join(joint_fit_clean, by=c("codon1", "codon2"))


###################
## Add tRNA data ##
###################

# merge tRNA isoacceptor count data for each codon we have info on
# NA codon entries are because not all codons were listed in the abundance data

## CODON 1 ##
SRC_tRNA <- SRC_rollup %>%
  left_join(tRNA, by=c("codon1"="codon", "Three_Letter_Code"="aa")) %>% 
  select(-aa) %>%
  rename("codon1_tRNA_count"=gene_count, 
         "anticodon1"=anticodon, 
         "aa"=Three_Letter_Code, 
         "codon1_cpm_rep1" = cpm_rep1, 
         "codon1_cpm_rep2" = cpm_rep2
         ) %>%
  mutate(anticodon1 = reverse_complement(codon1))

## CODON 2 ##
SRC_tRNA <- SRC_tRNA %>% 
  left_join(tRNA, by=c("codon2"="codon", "aa"="aa")) %>%
  rename("codon2_tRNA_count"=gene_count, 
         "anticodon2"=anticodon, 
         "codon2_cpm_rep1" = cpm_rep1, 
         "codon2_cpm_rep2" = cpm_rep2
         ) %>%
  mutate(anticodon2 = reverse_complement(codon2))


# calculate difference in isoacceptor abundance for each susbtitution
SRC_tRNA <- SRC_tRNA %>%
  mutate(count_diff_ratio = if_else((codon1_tRNA_count/codon2_tRNA_count) > 1, true=codon1_tRNA_count/codon2_tRNA_count, false=codon2_tRNA_count/codon1_tRNA_count)) %>%
  mutate(cpm_ratio_rep1 = if_else((codon1_cpm_rep1/codon2_cpm_rep1) > 1, true=codon1_cpm_rep1/codon2_cpm_rep1, false=codon2_cpm_rep1/codon1_cpm_rep1)) %>%
  mutate(cpm_ratio_rep2 = if_else((codon1_cpm_rep2/codon2_cpm_rep2) > 1, true=codon1_cpm_rep2/codon2_cpm_rep2, false=codon2_cpm_rep2/codon1_cpm_rep2))
# For substitutions where each codon in a pair is served by a different tRNA with a different abundance
# I would expect to see slower rates than for a pair of codons served by the same tRNA



## CORRELATION WITH MEAN RATES ##
# ignore NAs as missing data, for now
cor(SRC_tRNA$mean_rate, SRC_tRNA$count_diff_ratio, use="pairwise.complete.obs", method="spearman") 
cor(SRC_tRNA$mean_rate, SRC_tRNA$cpm_ratio_rep1, use="pairwise.complete.obs", method="spearman") 
cor(SRC_tRNA$mean_rate, SRC_tRNA$cpm_ratio_rep2, use="pairwise.complete.obs", method="spearman")

## CORRELATION WITH MEDIAN RATES ##
# ignore NAs as missing data, for now
cor(SRC_tRNA$median_rate, SRC_tRNA$count_diff_ratio, use="pairwise.complete.obs", method="spearman")
cor(SRC_tRNA$median_rate, SRC_tRNA$cpm_ratio_rep1, use="pairwise.complete.obs", method="spearman") 
cor(SRC_tRNA$median_rate, SRC_tRNA$cpm_ratio_rep2, use="pairwise.complete.obs", method="spearman")

## CORRELATION WITH JOINT RATES ##
# ignore NAs as missing data, for now
cor(SRC_tRNA$jointfit, SRC_tRNA$count_diff_ratio, use="pairwise.complete.obs", method="spearman")
cor(SRC_tRNA$jointfit, SRC_tRNA$cpm_ratio_rep1, use="pairwise.complete.obs", method="spearman") 
cor(SRC_tRNA$jointfit, SRC_tRNA$cpm_ratio_rep2, use="pairwise.complete.obs", method="spearman")


# cor(SRC_tRNA$codon1_tRNA_count, SRC_tRNA$codon1_cpm_rep1, use="pairwise.complete.obs", method="spearman")
# cor(SRC_tRNA$mean_rate, SRC_tRNA$median_rate, use="pairwise.complete.obs", method="spearman")


ggplot(SRC_tRNA, aes(x=mean_rate, 
                     y=cpm_ratio_rep1, 
                     )
) +
  geom_point() + 
  geom_text(aes(label=paste(aa, paste(codon1, codon2, sep=":"), sep="_")), hjust=0, vjust=0) +
  xlab("Mean SynREVCodon rate") +
  ylab("tRNA abundance ratio") +
  theme_bw()

ggplot(SRC_tRNA, aes(x=median_rate, 
                     y=cpm_ratio_rep1, 
                     label=paste(aa, paste(codon1, codon2, sep=":"), sep="_"))
) +
  geom_point() + 
  geom_text(hjust=1, vjust=0) +
  xlab("Median SynREVCodon rate") +
  ylab("tRNA abundance ratio") +
  theme_bw()


ggplot(SRC_tRNA, aes(x=jointfit, 
                     y=cpm_ratio_rep1
                     )
) +
  geom_point() + 
  geom_text(aes(label=paste(aa, paste(codon1, codon2, sep=":"), sep="_")), hjust=1, vjust=0) +
  xlab("Median SynREVCodon rate") +
  ylab("tRNA abundance ratio") +
  theme_bw()

# filter(SRC_tRNA, !(codon1 == "CTA" & codon2 == "CTT"))
# test <- filter(SRC_tRNA, !(codon1 == "CTA" & codon2 == "CTT"))
# cor(test$jointfit, test$cpm_ratio_rep1, use="pairwise.complete.obs", method="spearman")


###################################
## Add predicted tRNA abundances ##
###################################

# some tRNAs serve more than one codon. I've predicted which tRNAs serve which other codons
# copy tRNA abundance data and gene count data from "main codons" served by a tRNA
# associate that abundance data with other codons served by the same tRNA
# based on predictions of which codons share an isoacceptor
predictions_subset <- filter(SRC_tRNA, codon1 %in% tRNA_predictions$main_codon) %>%  # find present_codons in codon1 column to extract their data
  distinct(codon1, .keep_all=TRUE) %>%
  select(codon1, codon1_tRNA_count, codon1_cpm_rep1, codon1_cpm_rep2) %>%
  rename(main_codon="codon1", tRNA_count = "codon1_tRNA_count", cpm_rep1 = "codon1_cpm_rep1", cpm_rep2 = "codon1_cpm_rep2")

predictions <- filter(SRC_tRNA, codon2 %in% tRNA_predictions$main_codon) %>%
  distinct(codon2, .keep_all=TRUE) %>%
  select(codon2, codon2_tRNA_count, codon2_cpm_rep1, codon2_cpm_rep2) %>%
  rename(main_codon="codon2", tRNA_count = "codon2_tRNA_count", cpm_rep1 = "codon2_cpm_rep1", cpm_rep2 = "codon2_cpm_rep2") %>%
  bind_rows(predictions_subset) %>%
  distinct(main_codon, .keep_all = TRUE)

tRNA_predictions_data <- left_join(tRNA_predictions, predictions, by="main_codon")


# Merge and update the data with tRNA abundance predictions

## CODON 1 ##
updated_SRC_tRNA <- SRC_tRNA %>%
  replace(is.na(.), 0) %>%
  left_join(tRNA_predictions_data, by = c("codon1" = "other_codon", "aa"="aa")) %>%
  mutate(codon1_tRNA_count = if_else(is.na(tRNA_count), codon1_tRNA_count, codon1_tRNA_count + tRNA_count),
         codon1_cpm_rep1 = if_else(is.na(cpm_rep1), codon1_cpm_rep1, codon1_cpm_rep1 + cpm_rep1),
         codon1_cpm_rep2 = if_else(is.na(cpm_rep2), codon1_cpm_rep2, codon1_cpm_rep2 + cpm_rep2)
         ) %>%
  # rename(codon1_tRNA_mod="mod") %>%
  # mutate(codon1_tRNA_mod = replace_na(codon1_tRNA_mod, "none known")) %>%
  select(-tRNA_count, -main_codon, -cpm_rep1, -cpm_rep2, -mod)

## CODON 2 ##
updated_SRC_tRNA <- updated_SRC_tRNA %>%
  left_join(tRNA_predictions_data, by = c("codon2" = "other_codon", "aa"="aa")) %>%
  mutate(codon2_tRNA_count = if_else(is.na(tRNA_count), codon2_tRNA_count, codon2_tRNA_count + tRNA_count),
         codon2_cpm_rep1 = if_else(is.na(cpm_rep1), codon2_cpm_rep1, codon2_cpm_rep1 + cpm_rep1),
         codon2_cpm_rep2 = if_else(is.na(cpm_rep2), codon2_cpm_rep2, codon2_cpm_rep2 + cpm_rep2)
  ) %>%
  select(-tRNA_count, -main_codon, -cpm_rep1, -cpm_rep2, -mod)



# calculate difference in isoacceptor abundance for each susbtitution
updated_SRC_tRNA <- updated_SRC_tRNA %>%
  mutate(count_diff_ratio = if_else((codon1_tRNA_count/codon2_tRNA_count) > 1, true=codon1_tRNA_count/codon2_tRNA_count, false=codon2_tRNA_count/codon1_tRNA_count)) %>%
  mutate(cpm_ratio_rep1 = if_else((codon1_cpm_rep1/codon2_cpm_rep1) > 1, true=codon1_cpm_rep1/codon2_cpm_rep1, false=codon2_cpm_rep1/codon1_cpm_rep1)) %>%
  mutate(cpm_ratio_rep2 = if_else((codon1_cpm_rep2/codon2_cpm_rep2) > 1, true=codon1_cpm_rep2/codon2_cpm_rep2, false=codon2_cpm_rep2/codon1_cpm_rep2))
# For substitutions where each codon in a pair is served by a different tRNA with a different abundance
# I would expect to see slower rates than for a pair of codons served by the same tRNA



# Create an identifier for easy matching in tRNA_predictions_data
tRNA_predictions_data <- tRNA_predictions_data %>%
  mutate(join_key = pmin(main_codon, other_codon),
         join_key2 = pmax(main_codon, other_codon))

# Create corresponding identifiers in SRC_tRNA
updated_SRC_tRNA <- updated_SRC_tRNA %>%
  mutate(join_key = pmin(codon1, codon2),
         join_key2 = pmax(codon1, codon2))

# Perform the join
updated_SRC_tRNA <- updated_SRC_tRNA %>%
  left_join(tRNA_predictions_data, by = c("join_key", "join_key2", "aa")) %>%
  select(-join_key, -join_key2, -main_codon, -other_codon, -tRNA_count, -cpm_rep1, -cpm_rep2) %>%
  mutate(mod = replace_na(mod, "none"))
  # -everything from tRNA_predictions_data except mod and identifiers

# # Add the 'mod' column
# updated_SRC_tRNA_2 <- updated_SRC_tRNA_2 %>%
#   mutate(mod = coalesce(mod.x, mod.y)) %>%
#   select(-mod.x, -mod.y)  # Remove temporary columns if they were created during the join




## CORRELATION WITH MEAN RATES ##
#cor(SRC_tRNA$mean_rate, SRC_tRNA$count_diff_sq, use="pairwise.complete.obs", method="spearman") # ignore NAs as missing data, for now
cor(updated_SRC_tRNA$mean_rate, updated_SRC_tRNA$count_diff_ratio, use="pairwise.complete.obs", method="spearman") # ignore NAs as missing data, for now
cor(updated_SRC_tRNA$mean_rate, updated_SRC_tRNA$cpm_ratio_rep1, use="pairwise.complete.obs", method="spearman") # ignore NAs as missing data, for now
cor(updated_SRC_tRNA$mean_rate, updated_SRC_tRNA$cpm_ratio_rep2, use="pairwise.complete.obs", method="spearman") # ignore NAs as missing data, for now

## CORRELATION WITH MEDIAN RATES ##
#cor(SRC_tRNA$median_rate, SRC_tRNA$count_diff_sq, use="pairwise.complete.obs", method="spearman") # ignore NAs as missing data, for now
cor(updated_SRC_tRNA$median_rate, updated_SRC_tRNA$count_diff_ratio, use="pairwise.complete.obs", method="spearman")
cor(updated_SRC_tRNA$median_rate, updated_SRC_tRNA$cpm_ratio_rep1, use="pairwise.complete.obs", method="spearman") 
cor(updated_SRC_tRNA$median_rate, updated_SRC_tRNA$cpm_ratio_rep2, use="pairwise.complete.obs", method="spearman")

## CORRELATION WITH JOINT RATES ##
#cor(SRC_tRNA$median_rate, SRC_tRNA$count_diff_sq, use="pairwise.complete.obs", method="spearman") # ignore NAs as missing data, for now
cor(updated_SRC_tRNA$jointfit, updated_SRC_tRNA$count_diff_ratio, use="pairwise.complete.obs", method="spearman")
cor(updated_SRC_tRNA$jointfit, updated_SRC_tRNA$cpm_ratio_rep1, use="pairwise.complete.obs", method="spearman") 
cor(updated_SRC_tRNA$jointfit, updated_SRC_tRNA$cpm_ratio_rep2, use="pairwise.complete.obs", method="spearman")


# cor(SRC_tRNA$codon1_tRNA_count, SRC_tRNA$codon1_cpm_rep1, use="pairwise.complete.obs", method="spearman")
# cor(SRC_tRNA$mean_rate, SRC_tRNA$median_rate, use="pairwise.complete.obs", method="spearman")


ggplot(updated_SRC_tRNA, aes(x=mean_rate, 
                             y=cpm_ratio_rep1, 
                             #color=mod
                             )
       ) +
  geom_point() + 
  geom_text(aes(label=paste(aa, paste(codon1, codon2, sep=":"), sep="_")), hjust=0, vjust=0) +
  xlab("Mean SynREVCodon rate") +
  ylab("tRNA abundance ratio") +
  theme_bw()

ggplot(updated_SRC_tRNA, aes(x=median_rate, 
                             y=cpm_ratio_rep1, 
                             color=mod, 
                             label=paste(aa, paste(codon1, codon2, sep=":"), sep="_"))
       ) +
  geom_point() + 
  geom_text(hjust=1, vjust=0) +
  xlab("Median SynREVCodon rate") +
  ylab("tRNA abundance ratio") +
  theme_bw()



ggplot(updated_SRC_tRNA, aes(x=jointfit, 
                             y=cpm_ratio_rep1, 
                             color=mod
                             )
) +
  geom_point() + 
  #geom_text(aes(label=paste(aa, paste(codon1, codon2, sep=":"), sep="_")), hjust=1, vjust=0) +
  xlab("Median SynREVCodon rate") +
  ylab("tRNA abundance ratio") +
  theme_bw()





