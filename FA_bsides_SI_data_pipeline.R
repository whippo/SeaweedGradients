#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                                                                ##
# Antarctica Algae FA Bsides SI Data Pipeline                                    ##
# Script Created 2023-03-17                                                      ##
# Data source: Antarctic Gradients 2019                                          ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2023-12-29                                                        ##
#                                                                                ##
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# SUMMARY:

# Script to combine all current fatty acid biomarker data from the 'bside' list
# into a single spreadsheet with SI data and QAQC standardization of data. 


# Required Files:
# bside_algae_spp_final_prop.csv (no metadata, FA values)
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
# 2023-03-21 First complete working pipeline completed, may need to finesse

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LOAD PACKAGES                                                                ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(tidyverse)
library(stringr)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# READ IN AND PREPARE DATA                                                     ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


bsides_FA <- read_csv("Data/Biomarkers/FattyAcids/bside_algae_spp_final_prop.csv")
bsides_SI <- read_csv("Data/Biomarkers/SI/bside_full_list_KI.csv") 

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

# pull up all SI projectID to first column
SI_step2 <- SI_step1 %>%
  mutate(ProjID = coalesce(ProjID, `SI-ID`)) %>%
  drop_na(ProjID) %>%
  filter(!row_number() %in% c(127:131))
# fix incorrect site and GPS for sample 154 (incorrectly listed as site "I", 
# actual site is "M" crossreferenced by date
SI_step2[47,38] <- -68.17577
SI_step2[47,39] <- -67.26823
  

  
  


# JOIN DATASETS

# find differences in column names between FA and SI (in first vector, not in second)
diffSI <- SI_step2 %>%
  select(ProjID, `CN ratio`) %>%
  filter(!is.na(`CN ratio`))
FA_noSI <- as_tibble(setdiff(FA_step1$ProjID, diffSI$ProjID))
names(FA_noSI) <- "ProjID"
SI_noFA <- as_tibble(setdiff(diffSI$ProjID, FA_step1$ProjID))
names(SI_noFA) <- "ProjID"


# join FA and SI data with metadata by ProjID column
full_step1 <- SI_step2 %>%
  full_join(FA_step1, by = "ProjID")
# update species names to most current values and all site names
full_step2 <- full_step1 %>%
  mutate(siteName = case_when(SiteID == "12" ~ "A",
                              SiteID == "14" ~ "B",
                              SiteID == "11" ~ "C",
                              SiteID == "10" ~ "D",
                              SiteID == "9" ~ "E",
                              SiteID == "15" ~ "F",
                              SiteID == "8" ~ "G",
                              SiteID == "1" ~ "H",
                              SiteID == "13" ~ "I",
                              SiteID == "7" ~ "J",
                              SiteID == "6" ~ "K",
                              SiteID == "5" ~ "L",
                              SiteID == "3" ~ "M",
                              SiteID == "4" ~ "N",
                              SiteID == "2" ~ "X",
                              SiteID == "XX" ~ "XX",
                              sample == "IRCO_09F0739_0101" ~ "E",
                              sample == "IRCO_11F1030_0076" ~ "C",
                              sample == "IRCO_13F1251_0094" ~ "I",
                              sample == "PHAN_13F1231_0150" ~ "I",
                              sample == "PHAN_13F1230_0143" ~ "I")) %>%
  mutate(coreadditions = case_when(sample == "IRCO_09F0739_0101" ~ "Iridaea cordata",
                                   sample == "IRCO_11F1030_0076" ~ "Iridaea cordata",
                                   sample == "IRCO_13F1251_0094" ~ "Iridaea cordata",
                                   sample == "PHAN_13F1231_0150" ~ "Callophyllis atrosanguinea",
                                   sample == "PHAN_13F1230_0143" ~ "Callophyllis atrosanguinea")) %>%
  mutate(revisedSpecies = coalesce(coreadditions, `REVISED NAME`, genusSpecies...9, genusSpecies...36))
# reduce to required columns for analysis and delete comment lines
full_step3 <- full_step2 %>%
  select(ProjID, sample, `SI-ID`, LinkedTo, Comments, siteName, revisedSpecies, phylum,
         Date, Transect, Depth.m, Replicate, Tissue, VialAmount, batch:type,
         `Latitude (dec)`:`Ice cover (NIC-Midpoint-Annual)`, 
         `CN ratio`:d13C, `8:0`:`24:1w9`) %>%
  filter(!is.na(revisedSpecies))
# remove all FA columns that are zero
full_step4 <- full_step3 %>%
  select_if(~ !is.numeric(.) || sum(., na.rm = TRUE) != 0)
# fix species typo and uncertain ID
full_step5 <- full_step4 %>% 
  mutate(revisedSpecies = case_when(
    revisedSpecies == "Myriogramme mangini" ~ "Myriogramme manginii", TRUE ~ revisedSpecies)) %>%
  filter(revisedSpecies != "eliminate b/c uncertain ID")
# fill in missing phylum values (currently only Rhodophyta, check before changing this code)
full_step6 <- full_step5 %>%
  mutate(phylum1 = case_when(revisedSpecies == "Benthic diatoms" ~ "Ochrophyta",
                             revisedSpecies == "Meridionella antarctica" ~ "Rhodophyta",
                             revisedSpecies == "Georgiella confluens" ~ "Rhodophyta",
                             revisedSpecies == "Porphyra plocamiestris" ~ "Rhodophyta",
                             revisedSpecies == "Myriogramme smithii" ~ "Rhodophyta",
                             revisedSpecies == "Pachymenia orbicularis" ~ "Rhodophyta",
                             revisedSpecies == "Palmaria decipiens" ~ "Rhodophyta",
                             revisedSpecies == "Pantoneura plocamioides" ~ "Rhodophyta",
                             revisedSpecies == "Sarcopeltis antarctica" ~ "Rhodophyta",
                             revisedSpecies == "Iridaea cordata" ~ "Rhodophyta",
                             revisedSpecies == "Callophyllis atrosanguinea" ~ "Rhodophyta")) %>%
  mutate(phylum = coalesce(phylum, phylum1)) %>%
  select(!phylum1)

# save final joined FA dataset
gradients2019_bsides_FASI_QAQC <- full_step6
  
rm(full_step1, full_step2, full_step3, full_step4, full_step5, full_step6, diffSI, 
   FA_noSI, FA_step1, SI_noFA, SI_step1, SI_step2)
# write .csv with current joined data
write_csv(gradients2019_bsides_FASI_QAQC, "Data/Biomarkers/FattyAcids/gradients2019_bsides_FASI_QAQC_new.csv")


####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#

# SCRATCH PAD ####






