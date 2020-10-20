###################################################################################
#                                                                                ##
# B-236 Algae run for FA                                                         ##
# Data are current as of 2019-05-19                                              ##
# Data source: B-236 Seaweed Gradients Cruise                                    ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2020-04-09                                                        ##
#                                                                                ##
###################################################################################

# SUMMARY:

# Script to extract list of algae currently in the queue (or already extracted) for
# fatty acid analysis. Requires access to Ross Whippo's FA extraction log file in OSF.


# Required Files (check that script is loading latest version):
# Whippo_FA_extraction_log.csv

# Associated Scripts:
# none

# TO DO

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# RECENT CHANGES TO SCRIPT                                                        #
# LOAD PACKAGES                                                                   #
# READ IN AND PREPARE DATA                                                        #
# SPECIES NAMES & LIST                                                            #
# COLLECTION VISUALIZATIONS                                                       #
# IDENTIFY TARGET SPECIES                                                         #
# EXTRACT SAMPLE COLLECTION DATES                                                 #
#                                                                                 #
###################################################################################

###################################################################################
# RECENT CHANGES TO SCRIPT                                                        #
###################################################################################

# 2020-04-09 - script created

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

library(tidyverse)
library(worrms)
library(reshape2)
library(ggplot2)
library(viridis)



###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################

setwd("~/Git/SeaweedGradients/Data/Biomarkers/FattyAcids")

all_extractions <- read.csv("/home/ross/Dropbox/OSF/Fatty Acid Extractions/Whippo_FA_extraction_log.csv")

str(all_extractions)

# constrain to Gradients 2019 extractions
gradients_extractions <- all_extractions %>%
  filter(projectID == "Gradients2019")

# remove procedural blanks
gradients_noblanks <- gradients_extractions %>%
  filter(!is.na(genus))

# make genusSpecies column
gradients_noblanks <- gradients_noblanks %>%
  unite(genusSpecies, genus, species, sep = " ", remove = FALSE)

# list all species
spp_list <- unique(gradients_noblanks$genusSpecies)

# export list
write_csv(as.data.frame(spp_list), "all_algae")

#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#


######### SCRATCH PAD

