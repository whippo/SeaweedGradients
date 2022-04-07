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

# 03F0146 in metadata is duplicated, second one should be 03F0147?
# 12F1103-4 and 12F1111-12 mislabeled with 'S'

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
# 2022-03-31  Final first version of output file produced

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LOAD PACKAGES                                                                ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(tidyverse)
library(stringr)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# READ IN AND PREPARE DATA                                                     ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


algae_raw <- read_csv("Data/Biomarkers/FattyAcids/core_algae_spp_final_concs.csv")
invert_raw <- read_csv("Data/Biomarkers/FattyAcids/core_invert_proportions - Sheet1.csv")
sample_metadata <- read_csv("Data/Biomarkers/FattyAcids/B-236_Fatty-Acid_Collections_QAQC.csv")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# MANIPULATE DATA                                                              ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# PREPARE METADATASET FOR JOINING

# remove duplicated columns, make ice cover cat character
sample_metadata_step1 <- sample_metadata %>%
  select(!(X:LON_DD.1)) %>%
  mutate(IceCoverCat = as.character(IceCoverCat))
# fix ProjID names with incorrect values
sample_metadata_step1$ProjID <- sample_metadata_step1$ProjID %>%
  recode("12S1103" = "12F1103",
         "12S1104" = "12F1104",
         "12S1111" = "12F1111",
         "12S1112" = "12F1112")
sample_metadata_step1[120,1] <- "03F0147"
# Update phylogenetic labels
sample_metadata_step1$genusSpecies <- sample_metadata_step1$genusSpecies %>%
  recode("Gigartina skottsbergii" = "Sarcopeltis antarctica",
         "Rhodokrambe lanigioides" = "Myriogramme manginii")
sample_metadata_step1$Genus <- sample_metadata_step1$Genus %>%
  recode("Gigartina" = "Sarcopeltis",
         "Rhodokrambe" = "Myriogramme")
sample_metadata_step1$species <- sample_metadata_step1$species %>%
  recode("skottsbergii" = "antarctica",
         "lanigioides" = "manginii")
# create new column in metadata to join by sample number
sample_metadata_step2 <- sample_metadata_step1 %>%
  mutate(sampleNum = substr(ProjID, 4, 7)) %>%
  mutate(sampleNum = str_remove(sampleNum, "^0+"))
# reduce to name only columns
sample_metadata_step3 <- sample_metadata_step2 %>%
  select(ProjID, sampleNum)




# PREPARE ALGAL DATASET FOR JOINING


# extract ProjID value into new column and pivot data wider
algae_step1 <- algae_raw %>% 
  mutate(ProjID = str_extract(`sample`, '(?<=_).*?(?=_)')) %>%
  pivot_wider(names_from = FA, values_from = proportion, values_fill = 0)
# join algae data with metadata step 1 by ProjID column
algae_step2 <- algae_step1 %>%
  left_join(sample_metadata_step1, by = "ProjID")
# create sp.abrev and phylum columns 
algae_step3 <- algae_step2 %>%
  mutate(sp.abrev = substr(sample, 1, 4)) %>%
  mutate(phylum = case_when(sp.abrev == "DEME" ~ "ochrophyta",
                           sp.abrev == "HIGR" ~ "ochrophyta",
                           sp.abrev == "IRCO" ~ "rhodophyta",
                           sp.abrev == "PHAN" ~ "rhodophyta",
                           sp.abrev == "PLCA" ~ "rhodophtya",
                           sp.abrev == "MYMA" ~ "rhodophyta"))

# create list of colnames for comparison to inverts
algae_cols <- colnames(algae_step3)




# PREPARE INVERT DATASET FOR JOINING


# extract ProjID value into new column
invert_step1 <- invert_raw %>%
  mutate(ProjID = str_extract(`ID`, '.*?(?=_)'))
# filter out non-CORE project data
invert_step2 <- invert_step1 %>%
  filter(substr(ProjID, 1, 6) != "ANTSEA")
# fix ProjID names with incorrect leading zeros and bad values
invert_step2$ProjID <- invert_step2$ProjID %>%
  recode("0962" = "962",
         "0963" = "963",
         "0964" = "964",
         "0965" = "965",
         "0974" = "974",
         "0975" = "975",
         "0976" = "976",
         "0977" = "977",
         "0978" = "978",
         "0979" = "979",
         "0980" = "980",
         "0981" = "981",
         "0987" = "987",
         "0988" = "988",
         "0989" = "989",
         "0990" = "990",
         "0991" = "991",
         "0996" = "996",
         "0997" = "997",
         "0998" = "998",
         "0999" = "999",
         "01F0421" = "421",
         "01F0422" = "422",
         "01F0067" = "02F0067",
         "01F0068" = "02F0068",
         "01F0069" = "02F0069",
         "01F0070" = "02F0070",
         "01F0071" = "02F0071",
         "01F0072" = "02F0072",
         "01F0073" = "02F0073",
         "01F0074" = "02F0074",
         "01F0075" = "02F0075",
         "01F0076" = "02F0076",
         "01F0086" = "02F0086",
         "01F0087" = "02F0087",
         "01F0092" = "02F0092",
         "01F0093" = "02F0093",
         "01F0094" = "02F0094",
         "01F0100" = "02F0100",
         "01F0101" = "02F0101",
         "01F0112" = "02F0112",
         "01F0113" = "02F0113",
         "01F0114" = "02F0114",
         "01F0115" = "02F0115",
         "01F0116" = "02F0116",
         "01F0122" = "02F0122",
         "01F0123" = "02F0123",
         "01F0124" = "02F0124",
         "01F0130" = "02F0130",
         "01F0131" = "02F0131",
         "01F0216" = "03F0216",
         "01F0217" = "03F0217",
         "01F0218" = "03F0218")
# duplicate and rename ProjID column for joining
invert_step3 <- invert_step2 %>%
  mutate(sampleNum = ProjID)
# add ProjID from metadata column
invert_step4 <- invert_step3 %>%
  left_join(sample_metadata_step3, by = "sampleNum")
# create new column joining correct ProjID names
invert_step5 <- invert_step4 %>%
  mutate(ProjID.z = coalesce(ProjID.y, sampleNum))
# remove old ProjID names and rename column
invert_step6 <- invert_step5 %>%
  mutate(ProjID = ProjID.z) %>%
  select(-ProjID.x, -ProjID.y, -ProjID.z, -sampleNum)
# join invert data with metadata by ProjID column, drop Ice.cover column
invert_step7 <- invert_step6 %>%
  left_join(sample_metadata_step1, by = "ProjID") %>%
  select(!Ice.cover)
# create list of colnames for comparison to algae
invert_cols <- colnames(invert_step7)


# find differences in column names between algae and inverts
setdiff(algae_cols, invert_cols)
setdiff(invert_cols, algae_cols)


# remove "C" before colnames in inverts dataframe
colnames(invert_step7) <- str_replace(colnames(invert_step7), pattern = "[C](?=[12])", replacement = "") 
# create new list of colnames for comparison to algae
invert_cols <- colnames(invert_step7)
# check diff after change
setdiff(invert_cols, algae_cols)

shared_cols <- intersect(invert_cols, algae_cols)


# join inverts to algae and metadata across ProjID column and coalesce columns
joined_FA_step1 <- algae_step3 %>%
  full_join(invert_step7, by = shared_cols) %>%
  mutate(sample = coalesce(sample, ID)) %>%
  mutate(Depth.m = coalesce(Depth.m, depth)) %>%
  mutate(Genus = coalesce(Genus, genus)) %>%
  mutate(species = coalesce(species, species.x)) %>%
  mutate(Tissue = coalesce(Tissue, tissue.type)) %>%
  mutate(phylum = coalesce(phylum, phyla))
# drop redundant columns and move metadata to first columns
joined_FA_step2 <- joined_FA_step1 %>%
  select(-c("ID", 
            "Type",
            "EnteredBy",
            "DATE",
            "depth", 
            "genus", 
            "tissue.type", 
            "species.x", 
            "FAvialNum.y",
            "IceCoverCat",
            "species.y",
            "Comments",
            "phyla")) %>%
  relocate(SiteID:season, .after = ProjID) 
# join FAvialNum columns together
joined_FA_step3 <- joined_FA_step2 %>%
  mutate(FAvialNum = coalesce(FAvialNum, FAvialNum.x))
# fill NA in FA values with 0
joined_FA_step4 <- joined_FA_step3 %>%
  mutate(across("8:0":"22:4w3", ~replace_na(., 0)))
# drop reduntant vial number column
joined_FA_step5 <- joined_FA_step4 %>%
  select(-FAvialNum.x)

# make sure all columns have been sorted and selected properly
colnames(joined_FA_step5)


# save final joined FA dataset with different name
gradients2019_corespecies_FA_QAQC <- joined_FA_step5

# write .csv with current joined data
write_csv(gradients2019_corespecies_FA_QAQC, "Data/Biomarkers/FattyAcids/gradients2019_corespecies_FA_QAQC.csv")

####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#

# SCRATCH PAD ####