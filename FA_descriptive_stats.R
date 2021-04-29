#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                                                                ##
# Antarctica Algae Descriptive Stats                                             ##
# Data are current as of 2021-04-26                                              ##
# Data source: Antarctic Gradients 2019                                          ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2021-04-26                                                        ##
#                                                                                ##
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# SUMMARY:

# Script to generate descriptive statistical summaries of fatty acid data collected
# from Antarctic Gradients 2019 project. 


# Required Files (check that script is loading latest version):
# dummy_sample_values.csv


# Associated Scripts:
# FILE.R

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

# 2021-04-26 Script created

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LOAD PACKAGES                                                                ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(stringr)
library(vegan)
library(viridis)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# READ IN AND PREPARE DATA                                                     ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

dummy_data <- read_csv("Data/Biomarkers/FattyAcids/dummy_sample_values.csv")

Gradients19_FA_Concs <- read_csv("Data/Biomarkers/FattyAcids/Gradients19_FA_Concs.csv", 
                                                 col_types = cols(Conc = col_double(), 
                                                 Date.anal = col_character(), Notes = col_character()))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# MANIPULATE DATA                                                              ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

############### Dummy Data

# convert data to long form
dummy_trans <- t(dummy_data)
long_dummy <- dummy_trans %>%
  pivot_longer(!X1, names_to = "")

############### Actual data

# create species-genus ID column
grad_conc <- Gradients19_FA_Concs %>%
  mutate(species = str_replace_all(Sample, c("HIGR_10F0836_0130_5252020_22" = "H. grandifolius",
                                            "MYMA_09F0778_0166_5252020_7" = "M. manginii",
                                            "LAAN_14F1207_0194_5252020_14" = "L. antarctica")))

grad_conc <- grad_conc %>% 
  drop_na(Conc)

# replace NAs with 0 in Conc
grad_conc$Conc[is.na(grad_conc$Conc)] <- 0

### MDS

# pivot data wide for mds
grad_conc_wide <- grad_conc %>%
  select(Name, Conc, species) %>%
  pivot_wider(names_from = Name, values_from = Conc, values_fill = 0)

metaMDS(grad_conc_wide[2:14])

### Figure 

ggplot(filter(grad_conc, Conc > 5), aes(x = Name, y = Conc, colour = species)) +
  geom_point(size = 4) +
  theme_classic() +
  scale_colour_viridis(discrete = TRUE, end = 0.9) +
  labs(x = "Fatty Acid", y = "Concentration (ng/ul)")
  



####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#

# SCRATCH PAD ####