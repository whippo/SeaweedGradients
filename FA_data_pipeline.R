#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                                                                ##
# Antarctica Algae FA Data Pipeline                                              ##
# Data are current as of 2022-03-24                                              ##
# Data source: Antarctic Gradients 2019                                          ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2022-03-24                                                        ##
#                                                                                ##
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# SUMMARY:

# Script to combine all current fatty acid biomarker data into single spreadsheet
# with QAQC and standardization of data. 


# Required Files:
# core_algae_spp_final_concs.csv
# core_invert_proportions - Sheet1.csv
# B-236_Fatty-Acid_Collections.csv


# Associated Scripts:
# none

# TO DO 

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

# 2022-03-24  Script created

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LOAD PACKAGES                                                                ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(tidyverse)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# READ IN AND PREPARE DATA                                                     ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


algae_raw <- read_csv("Data/Biomarkers/FattyAcids/core_algae_spp_final_concs.csv")
invert_raw <- read_csv("Data/Biomarkers/FattyAcids/core_invert_proportions - Sheet1.csv")
sample_metadata <- read_csv("Data/Biomarkers/FattyAcids/B-236_Fatty-Acid_Collections.csv")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# MANIPULATE DATA                                                              ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# PREPARE ALGAL DATASET FOR JOINING

# extract ProjID value into new column and pivot data wider
algae_step1 <- algae_raw %>% 
  mutate(ProjID = str_extract(`sample`, '(?<=_).*?(?=_)')) %>%
  pivot_wider(names_from = FA, values_from = proportion, values_fill = 0)
# join algae data with metadata by ProjID column
algae_step2 <- algae_step1 %>%
  left_join(sample_metadata, by = "ProjID")
# create list of colnames for comparison to inverts
algae_cols <- colnames(algae_step2)


# PREPARE INVERT DATASET FOR JOINING

# extract ProjID value into new column
invert_step1 <- invert_raw %>%
  mutate(ProjID = str_extract(`ID`, '.*?(?=_)'))
# filter out non-CORE project data
invert_step2 <- invert_step1 %>%
  filter(str_length(ProjID) > 6)
# create list of colnames for comparison to algae
invert_cols <- colnames(invert_step2)


# find differences in column names between algae and inverts
setdiff(algae_cols, invert_cols)
setdiff(invert_cols, algae_cols)


# remove "C" before colnames in inverts dataframe
colnames(invert_step2) <- gsub(colnames(invert_step2), pattern = "C", replacement = "") 
# create new list of colnames for comparison to algae
invert_cols <- colnames(invert_step2)
# check diff after change
setdiff(invert_cols, algae_cols)

shared_cols <- intersect(invert_cols, algae_cols)


# join inverts to algae and metadata across ProjID column and coalesce columns
joined_FA_step1 <- algae_step2 %>%
  full_join(invert_step2, by = shared_cols) %>%
  mutate(sample = coalesce(sample, ID)) %>%
  mutate(IceCoverCat = coalesce(IceCoverCat, Ice.cover)) %>%
  mutate(Depth.m = coalesce(Depth.m, depth)) %>%
  mutate(Genus = coalesce(Genus, genus)) %>%
  mutate(Tissue = coalesce(Tissue, tissue.type))
# drop redundant columns and move metadata to first columns
joined_FA_step2 <- joined_FA_step1 %>%
  select(-c("ID", "Ice.cover", "depth", "genus", "tissue.type")) %>%
  relocate(SiteID:sp.abrev, .after = ProjID) 
# fill NA in FA values with 0
joined_FA_step3 <- joined_FA_step2 %>%
  mutate(across("8:0":"22:4w3", ~replace_na(., 0)))

# make sure all columns have been sorted and selected properly
colnames(joined_FA_step3)

####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#

# SCRATCH PAD ####