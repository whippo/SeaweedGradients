###################################################################################
#                                                                                ##
# B-236 Fatty Acid Collections Data                                              ##
# Data are current as of 2019-05-12                                              ##
# Data source: B-236 Seaweed Gradients Cruise                                    ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2019-05-13                                                        ##
#                                                                                ##
###################################################################################

# SUMMARY:

# Script for QAQC of invertebrate collection data made from the B-236 Seaweed Gradients
# Antarctic cruise from 2019-04-12:2019-05-18 including expedition members from the 
# University of Oregon (OIMB), University of Alabama Birmingham, and Texas A&M.


# Required Files (check that script is loading latest version):
# B-236_Fatty-Acid_Collections.csv

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
# 2019-05-13 Changed name of source files

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

raw_inverts_0 <- read.csv("B-236_Fatty-Acid_Collections.csv")

str(raw_inverts_0)

# remove empty rows
raw_inverts_1 <- raw_inverts_0 %>%
  filter(SiteID != "")

# recode ice cover cats

raw_inverts_1$IceCoverCat <- raw_inverts_1$IceCoverCat %>%
  recode("<=60%" = "0.6",
         "<=80%" = "0.8",
         "<=90%" = "0.9",
         "50%" = "0.5",
         "80%" = "0.8",
         "<=40%" = "0.4")
raw_inverts_1$IceCoverCat <- as.numeric(as.character(raw_inverts_1$IceCoverCat))

# fix ice cover cat mistake in site 9
raw_inverts_1[679:680, 6] <- 0.6

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
# join unrecognized species
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
raw_inverts_2$genusSpecies <- raw_inverts_2$genusSpecies %>%
  recode("Neosmilaster georginus" = "Neosmilaster georgianus", 
         "Sterechinus neumeyeri" = "Sterechinus neumayeri",
         "Isotaelia antarctica" = "Isotealia antarctica",
         "Isotealia antactica" = "Isotealia antarctica",
         "Isotelia autare" = "Isotealia antarctica",
         "Heterocucumis cucumaris" = "Heterocucumis steineni",
         "Porania antarctica glabra" = "Glabraster antarctica",
         "Alcyonium antarctucum" = "Alcyonium antarcticum",
         "Metalaptamphopus pectinatus" = "Metaleptamphopus pectinatus",
         "Labidaster annulatus" = "Labidiaster annulatus",
         "Promachocrinus kerguelenensis" = "Promachocrinus kerguelensis",
         "Phorbus aerolatus" = "Phorbas areolatus",
         "Metalaptamphopus pectatus" = "Metaleptamphopus pectinatus",
         "Bovalia gigantea" = "Bovallia gigantea",
         "Austrodoris vergulensis" = "Doris kerguelenensis",
         "Austrodoris kerguelenensis" = "Doris kerguelenensis",
         "Pontogeneilla brevicornis" = "Prostebbingia brevicornis",
         "Panaceradocus miersi" = "Paraceradocus miersi",
         "Polynoidae so" = "Polynoidae sp.",
         "Desmerestia menziesii" = "Desmarestia menziesii",
         "Paraphimedia integricauda" = "Pariphimedia integricauda",
         "Oraderea bidentata" = "Oradarea bidentata"
         )



####################### make species list

species_list <-raw_inverts_2 %>% 
  group_by(genusSpecies) %>% 
  tally()

write_csv(species_list, "B236_fattyacid-species-list.csv")

####################### species found at all sites

Site1 <- raw_inverts_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "01")
Site2 <- raw_inverts_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "02")
Site3 <- raw_inverts_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "03")
Site4 <- raw_inverts_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "04")
Site5 <- raw_inverts_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "05")
Site6 <- raw_inverts_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "06")
Site7 <- raw_inverts_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "07")
Site8 <- raw_inverts_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "08")
Site9 <- raw_inverts_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "09")
Site10 <- raw_inverts_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "10")

Reduce(intersect, list(Site1$genusSpecies, Site2$genusSpecies, Site3$genusSpecies,
                       Site4$genusSpecies, Site5$genusSpecies, Site6$genusSpecies,
                       Site7$genusSpecies, Site8$genusSpecies, Site9$genusSpecies,
                       Site10$genusSpecies))

####################### extract algal species

specieslist1 <- species_list[1:50,1]
specieslist2 <- species_list[51:69,1]

alltaxa <- wm_records_taxamatch(specieslist1$genusSpecies)
alltaxa2 <- wm_records_taxamatch(specieslist2$genusSpecies)

taxatable <- alltaxa %>% 
  map_df(tibble::rownames_to_column, .id = 'kingdom')
taxatable <- taxatable[,2:29]

taxatable2 <- alltaxa2 %>%
  map_df(tibble::rownames_to_column, .id = 'kingdom')
taxatable2 <- taxatable2[,2:29]

finaltaxatable <- bind_rows(taxatable, taxatable2)

algaetable <- finaltaxatable %>%
  filter(kingdom != 'Animalia')

write.csv(algaetable[,4], "algaetable.csv")



###################################################################################
# COLLECTION VISUALIZATIONS                                                       #
###################################################################################

raw_inverts_f1 <- raw_inverts_2
raw_inverts_f1$count <- 1

# per individual
raw_inverts_f2 <- raw_inverts_f1 %>%
  group_by(IceCoverCat, SiteID) %>%
  summarise(sum(count))
# rename column
names(raw_inverts_f2)[names(raw_inverts_f2)=="sum(count)"] <- "totalOrganisms"

# per species
raw_inverts_f3 <- raw_inverts_f1 %>%
  group_by(IceCoverCat, genusSpecies, SiteID) %>% 
  summarise(sum(count))
# rename column
names(raw_inverts_f3)[names(raw_inverts_f3)=="sum(count)"] <- "totalSpecies"
raw_inverts_f3$totalSpecies <- 1


# total number of organisms collected for biomarkers per ice cover
ggplot(data = raw_inverts_f2, aes(x = IceCoverCat, y = totalOrganisms)) +
  geom_col() +
  theme_classic()

# total number of species collected for biomarkers per ice cover
ggplot(data = raw_inverts_f3, aes(x = IceCoverCat, y = totalSpecies)) +
  geom_col() +
  theme_classic()

# total number of organisms collected for biomarkers per site
ggplot(data = raw_inverts_f2, aes(x = SiteID, y = totalOrganisms)) +
  geom_col() +
  theme_classic()

# total number of species collected for biomarkers per site
ggplot(data = raw_inverts_f3, aes(x = SiteID, y = totalSpecies)) +
  geom_col() +
  theme_classic()

raw_inverts_f4 <- raw_inverts_f1 %>%
  group_by(SiteID, genusSpecies) %>%
  summarize(sum(count))
# rename column
names(raw_inverts_f4)[names(raw_inverts_f4)=="sum(count)"] <- "totalSpecies"

# rough plot of species by site
ggplot(data = raw_inverts_f3, aes(x = SiteID, y = genusSpecies)) +
  geom_count() 

# heatmap
ggplot(data = raw_inverts_f4, aes(x = SiteID, y = genusSpecies)) +
  geom_tile(aes(fill = totalSpecies)) +
  scale_fill_viridis(option = "B", begin = 0.3, end = 1) +
  theme_dark() +
  theme(panel.background = element_rect(fill = 'black')) 
  
str(raw_inverts_f4)




############### SUBSECTION HERE


#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#