###################################################################################
#                                                                                ##
# B-236 Invertebrate Collections Data                                            ##
# Data are current as of 2019-05-12                                              ##
# Data source: B-236 Seaweed Gradients Cruise                                    ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2019-05-12                                                        ##
#                                                                                ##
###################################################################################

# SUMMARY:

# Script for QAQC of invertebrate collection data made from the B-236 Seaweed Gradients
# Antarctic cruise from 2019-04-12:2019-05-18 including expedition members from the 
# University of Oregon (OIMB), University of Alabama Birmingham, and Texas A&M.


# Required Files (check that script is loading latest version):
# B-236_Invert_Collections.csv

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
#                                                                                 #
###################################################################################

###################################################################################
# RECENT CHANGES TO SCRIPT                                                        #
###################################################################################

# 2019-05-12 Script created and added to git repository:
# https://github.com/whippo/SeaweedGradients.git

###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

library(tidyverse)
library(worrms)
library(reshape2)

###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################

setwd("~/Git/SeaweedGradients/Data/Inverts")

raw_inverts_0 <- read.csv("B-236_Invert_Collections.csv")

str(raw_inverts_0)

# remove empty rows
raw_inverts_1 <- raw_inverts_0 %>%
  filter(SiteID != "NA")


###################################################################################
# SPECIES NAMES & LIST                                                            #
###################################################################################

################ identify incorrect species names

# create genusSpecies column
raw_inverts_1 <- raw_inverts_1 %>%
  unite(genusSpecies, Genus, species, sep = " ", remove = FALSE)

# check names against WoRMS database
taxon_names <- raw_inverts_1$genusSpecies
taxon_names <- unique(taxon_names)
w <-wm_name2id_(taxon_names)
# print unrecognized returns
w1 <- melt(w)
# join nunrecognized species
taxon_table <- melt(taxon_names)
names(taxon_table)[names(taxon_table)=="value"] <- "L1"
species_notfound <- setdiff(taxon_table$L1, w1$L1)
species_notfound <- melt(species_notfound)
names(species_notfound)[names(species_notfound)=="value"] <- "L1"
failed_species <- w1 %>%
  filter(value == 1)

incorrect_spnames <- full_join(species_notfound, failed_species)

####################### correct species names

# fix incorrect species names
raw_inverts_2 <- raw_inverts_1
raw_inverts_2$genusSpecies %>%
  recode("Neosmilaster georginus" = "Neosmilaster georgianus", 
         "Sterechinus neumeyeri" = "Sterechinus neumayeri",
         "Isotaelia antarctica" = "Isotealia antarctica",
         "Heterocucumis cucumaris" = "Heterocucumis steineni",
         "Porania antarctica glabra" = "Glabraster antarctica")



####################### make species list

species_list <-raw_inverts_1 %>% 
  group_by(genusSpecies) %>% 
  tally()

write_csv(species_list, "B236_invert-species-list.csv")



############### SUBSECTION HERE

#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#