#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                                                                ##
# Antarctica Algae FA Bsides SI Data Pipeline                                    ##
# Script Created 2023-03-17                                                      ##
# Data source: Antarctic Gradients 2019                                          ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2023-03-21                                                        ##
#                                                                                ##
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# SUMMARY:

# Script to combine all current fatty acid biomarker data from the 'bside' list
# into a single spreadsheet with SI data and QAQC standardization of data. 


# Required Files:
# bside_algae_spp_final_concs.csv (no metadata, FA values)
# bside_full_list_KI.csv (updated metadata, SI values)


# Associated Scripts:
# none

# TO DO 

# 

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# TABLE OF CONTENTS                                                            ####
#                                                                                 +
# RECENT CHANGES TO SCRIPT                                                        +
# LOAD PACKAGES                                                                   +
# READ IN AND PREPARE DATA                                                        +
# MANIPULATE DATA                                                                 +
#                                                                                 +
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# RECENT CHANGES TO SCRIPT                                                     ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 2023-03-17 Script created

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LOAD PACKAGES                                                                ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(tidyverse)
library(stringr)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# READ IN AND PREPARE DATA                                                     ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


bsides_FA <- read_csv("Data/Biomarkers/FattyAcids/bside_algae_spp_final_concs.csv")
bsides_SI <- read_csv("Data/Biomarkers/SI/bside_full_list_KI.csv") %>%
  drop_na(ProjID)
core_FA <- read_csv("Data/Biomarkers/FattyAcids/gradients2019_corespecies_FA_QAQC.csv")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# MANIPULATE DATA                                                              ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# PREPARE FA DATASET FOR JOINING

# extract ProjID value into new column and pivot data wider
FA_step1 <- bsides_FA %>% 
  mutate(ProjID = str_extract(`sample`, '(?<=_).*?(?=_)')) %>%
  pivot_wider(names_from = FA, values_from = proportion, values_fill = 0)


# PREPARE SI DATASET WITH METADATA FOR JOINING

# extract ProjID value into new column
SI_step1 <- bsides_SI
  


# JOIN DATASETS

# find differences in column names between FA and SI (in first vector, not in second)
diffSI <- SI_step1 %>%
  select(ProjID, `CN ratio`) %>%
  filter(!is.na(`CN ratio`))
FA_noSI <- as_tibble(setdiff(FA_step1$ProjID, diffSI$ProjID))
names(FA_noSI) <- "ProjID"
SI_noFA <- as_tibble(setdiff(diffSI$ProjID, FA_step1$ProjID))
names(SI_noFA) <- "ProjID"


# join FA and SI data with metadata by ProjID column
full_step1 <- SI_step1 %>%
  left_join(FA_step1, by = "ProjID")



# make sure all columns have been sorted and selected properly
colnames(full_step1)

# separate out my FA sample names w/ Katrin's sample IDs
species <- full_step1 %>%
  select(`ProjID`, `genusSpecies...9`, `genusSpecies...35`, `sample`, `CN ratio`, `8:0`)

# which samples do I have that Katrin doesn't?
FA_noSI <- full_step1 %>%
  select(`genusSpecies...35`, `sample`) %>%
  filter(`genusSpecies...35` =="n/a") %>%
  select(`sample`)



# which samples does Katrin have that I don't?
SI_noFA <- full_step1 %>%
  select(`ProjID`, `genusSpecies...35`, `sample`) %>%
  filter(is.na(`sample`)) %>%
  select(`ProjID`, `genusSpecies...35`)

# use this list of SI that is not in my bsides to search my core for those samples

FAbsides_additions <- SI_noFA %>%
  select(ProjID) %>%
  left_join(core_FA, by = "ProjID") %>%
  select(ProjID, FAsampleName, batch, `8:0`:`22:4w3`)









full_step2 <- full_step1 %>%
  full_join(FAbsides_additions, by = intersect("ProjID")) %>%
  group_by(ProjID) %>%
  summarize_all(na.omit)
















# save final joined FA dataset with different name
gradients2019_bsides_FA-SI_QAQC <- full_step1

# write .csv with current joined data
write_csv(gradients2019_bsides_FA-SI_QAQC, "Data/Biomarkers/FattyAcids/gradients2019_bsides_FA-SI_QAQC.csv")

####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#

# SCRATCH PAD ####

# temp script to extract list of every species and rep count per site
species_list <- gradients2019_corespecies_FA_QAQC %>%
  select(LetterID, genusSpecies, phylum, Replicate) %>%
  distinct() %>%
  group_by(genusSpecies, LetterID, phylum) %>%
  top_n(1, abs(Replicate)) %>%
  arrange(LetterID, phylum)

# temp script to copy to bsides script extracting full list of algae used in bsides MS
bside_algae_list <- read_csv("Data/Biomarkers/FattyAcids/bside_algae_list.csv")

# join code list with species and collection attributes
bside_algae_list <- bside_algae_list %>%
  rename("ProjID" = "code")
bside_full <- left_join(bside_algae_list, sample_metadata_step2, by = "ProjID") %>%
  ungroup() %>%
  arrange(genusSpecies)

write_csv(bside_full, "bside_full_list.csv")





