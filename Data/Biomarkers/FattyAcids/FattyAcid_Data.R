###################################################################################
#                                                                                ##
# B-236 Fatty Acid Collections Data                                              ##
# Data are current as of 2019-05-19                                              ##
# Data source: B-236 Seaweed Gradients Cruise                                    ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2019-05-19                                                        ##
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
# COLLECTION VISUALIZATIONS                                                       #
# IDENTIFY TARGET SPECIES                                                         #
#                                                                                 #
###################################################################################

###################################################################################
# RECENT CHANGES TO SCRIPT                                                        #
###################################################################################

# 2019-05-12 Script created and added to git repository:
# https://github.com/whippo/SeaweedGradients.git
# 2019-05-13 Changed name of source files
# 2019-05-19 Added phylum extractions, worked on IDing common species

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

raw_biotaxa_0 <- read.csv("B-236_Fatty-Acid_Collections.csv")

str(raw_biotaxa_0)

# remove empty rows
raw_biotaxa_1 <- raw_biotaxa_0 %>%
  filter(SiteID != "")

# recode ice cover cats

raw_biotaxa_1$IceCoverCat <- raw_biotaxa_1$IceCoverCat %>%
  recode("<=60%" = "0.6",
         "<=80%" = "0.8",
         "<=90%" = "0.9",
         "50%" = "0.5",
         "80%" = "0.8",
         "<=40%" = "0.4",
         "<=70%" = "0.7")
raw_biotaxa_1$IceCoverCat <- as.numeric(as.character(raw_biotaxa_1$IceCoverCat))

# fix ice cover cat mistake in site 9
raw_biotaxa_1[679:680, 6] <- 0.6

###################################################################################
# SPECIES NAMES & LIST                                                            #
###################################################################################

################ identify incorrect species names

# create genusSpecies column
raw_biotaxa_1 <- raw_biotaxa_1 %>%
  unite(genusSpecies, Genus, species, sep = " ", remove = FALSE)

# check names against WoRMS database (exact match only)
taxon_names <- raw_biotaxa_1$genusSpecies

# start from here for recheck #

taxon_names <- unique(taxon_names)
# can only check 50 at a time
taxacheck1 <- wm_records_names(taxon_names[1:50])
taxacheck2 <- wm_records_names(taxon_names[51:100])
taxacheck3 <- wm_records_names(taxon_names[101:length(taxon_names)])

# join tibbles in table
taxatable1 <- taxacheck1 %>% 
  map_df(tibble::rownames_to_column, .id = 'kingdom')
taxatable1 <- taxatable1[,2:29]

taxatable2 <- taxacheck2 %>%
  map_df(tibble::rownames_to_column, .id = 'kingdom')
taxatable2 <- taxatable2[,2:29]

taxatable3 <- taxacheck3 %>%
  map_df(tibble::rownames_to_column, .id = 'kingdom')
taxatable3 <- taxatable3[,2:29]

checkedtaxatable <- bind_rows(taxatable1, taxatable2, taxatable3)

# ID species not returned and join in table
species_notreturned <- setdiff(taxon_names, checkedtaxatable$scientificname)
species_notreturned <- melt(species_notreturned)
species_notreturned$status <- "unaccepted"
names(species_notreturned)[names(species_notreturned)=="value"] <- "scientificname"
species_notreturned$scientificname <- as.character(species_notreturned$scientificname)

finaltaxacheck <- bind_rows(checkedtaxatable, species_notreturned, .id = "scientificname")

# extract species not found
unacceptedtaxa <- filter(finaltaxacheck[,2:length(finaltaxacheck)], status != "accepted")

# possible candidates for not returned values with fuzzy matching
possibletaxa <- wm_records_taxamatch(species_notreturned$scientificname)
possibletaxa <- possibletaxa %>% 
  map_df(tibble::rownames_to_column, .id = 'kingdom')
possibletaxa <- possibletaxa[,2:29]

####################### correct species names

# fix incorrect species names
raw_biotaxa_2 <- raw_biotaxa_1
raw_biotaxa_2$genusSpecies <- raw_biotaxa_2$genusSpecies %>%
  recode("Neosmilaster georginus" =        "Neosmilaster georgianus", 
         "Sterechinus neumeyeri" =         "Sterechinus neumayeri",
         "Isotaelia antarctica" =          "Isotealia antarctica",
         "Isotealia antactica" =           "Isotealia antarctica",
         "Isotelia autare" =               "Isotealia antarctica",
         "Heterocucumis cucumaris" =       "Heterocucumis steineni",
         "Porania antarctica glabra" =     "Glabraster antarctica",
         "Porania antarctica" =            "Glabraster antarctica",
         "Alcyonium antarctucum" =         "Alcyonium antarcticum",
         "Metalaptamphopus pectinatus" =   "Metaleptamphopus pectinatus",
         "Labidaster annulatus" =          "Labidiaster annulatus",
         "Promachocrinus kerguelenensis" = "Promachocrinus kerguelensis",
         "Phorbus aerolatus" =             "Phorbas bergmontae",
         "Metalaptamphopus pectatus" =     "Metaleptamphopus pectinatus",
         "Bovalia gigantea" =              "Bovallia gigantea",
         "Austrodoris vergulensis" =       "Doris kerguelenensis",
         "Austrodoris kerguelenensis" =    "Doris kerguelenensis",
         "Pontogeneilla brevicornis" =     "Prostebbingia brevicornis",
         "Panaceradocus miersi" =          "Paraceradocus miersi",
         "Polynoidae so" =                 "Polynoidae sp.",
         "Desmerestia menziesii" =         "Desmarestia menziesii",
         "Paraphimedia integricauda" =     "Pariphimedia integricauda",
         "Oraderea bidentata" =            "Oradarea bidentata",
         "Myriogramme mangini" =           "Myriogramme manginii",
         "Waldeckia obesa" =               "Charcotia obesa",
         "Piconiella plumosa" =            "Picconiella plumosa",
         "Rhodokrambe laingioides" =       "Rhodokrambe lanigioides",
         "Eusiris sp." =                   "Eusirus sp.",
         "Cysctosphaera jacquinotii" =     "Cystosphaera jacquinotii",
         "Cytstoclonium obtusangulum" =    "Cystoclonium obtusangulum",
         "Ballia callitricna" =            "Ballia callitricha",
         "Palmeria decipiens" =            "Palmaria decipiens",
         "Cordiea racovitzae" =            "Curdiea racovitzae",
         "Psilaster chascoti" =            "Psilaster charcoti",
         "Plocamium antarcticum" =         "Plocamium cartilagineum",
         "Astropugettia crassa" =          "Austropugetia crassa",
         "Pachymenia obicularis" =         "Pachymenia orbicularis",
         "Porphyra plocamienstris" =       "Porphyra plocamiestris",
         "Laevilacunana antarctica" =      "Laevilacunaria antarctica"
         
  )

# recheck all names
taxon_names <- unique(raw_biotaxa_2$genusSpecies)
# go back to recheck

####################### make species list

species_list <-raw_biotaxa_2 %>% 
  group_by(genusSpecies) %>% 
  tally()

write_csv(species_list, "B236_fattyacid-species-list.csv")

####################### extract algal species

algaetable <- checkedtaxatable %>%
  filter(kingdom != 'Animalia')
algaelist <- unique(algaetable$valid_name)

write.csv(algaelist, "algaetable.csv")


####################### species found at all sites

Site1 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "01")
Site2 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "02")
Site3 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "03")
Site4 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "04")
Site5 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "05")
Site6 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "06")
Site7 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "07")
Site8 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "08")
Site9 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "09")
Site10 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "10")
Site11 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "11")
Site12 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "12")
Site13 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "13")
Site14 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "14")
Site15 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "15")

Reduce(intersect, list(Site1$genusSpecies, Site2$genusSpecies, Site3$genusSpecies,
                       Site4$genusSpecies, Site5$genusSpecies, Site6$genusSpecies,
                       Site7$genusSpecies, Site8$genusSpecies, Site9$genusSpecies,
                       Site10$genusSpecies, Site11$genusSpecies, Site12$genusSpecies,
                       Site13$genusSpecies, Site14$genusSpecies, Site15$genusSpecies))


###################################################################################
# COLLECTION VISUALIZATIONS                                                       #
###################################################################################

raw_biotaxa_f1 <- raw_biotaxa_2
raw_biotaxa_f1$count <- 1

# per individual
raw_biotaxa_f2 <- raw_biotaxa_f1 %>%
  group_by(IceCoverCat, SiteID) %>%
  summarise(sum(count))
# rename column
names(raw_biotaxa_f2)[names(raw_biotaxa_f2)=="sum(count)"] <- "totalOrganisms"

# per species
raw_biotaxa_f3 <- raw_biotaxa_f1 %>%
  group_by(IceCoverCat, genusSpecies, SiteID) %>% 
  summarise(sum(count))
# rename column
names(raw_biotaxa_f3)[names(raw_biotaxa_f3)=="sum(count)"] <- "totalSpecies"
raw_biotaxa_f3$totalSpecies <- 1


# total number of organisms collected for biomarkers per ice cover
ggplot(data = raw_biotaxa_f2, aes(x = IceCoverCat, y = totalOrganisms)) +
  geom_col() +
  theme_classic()

# total number of species collected for biomarkers per ice cover
ggplot(data = raw_biotaxa_f3, aes(x = IceCoverCat, y = totalSpecies)) +
  geom_col() +
  theme_classic()

# total number of organisms collected for biomarkers per site
ggplot(data = raw_biotaxa_f2, aes(x = SiteID, y = totalOrganisms)) +
  geom_col() +
  theme_classic()

# total number of species collected for biomarkers per site
ggplot(data = raw_biotaxa_f3, aes(x = SiteID, y = totalSpecies)) +
  geom_col() +
  theme_classic()



raw_biotaxa_f4 <- raw_biotaxa_f1 %>%
  group_by(SiteID, genusSpecies) %>%
  summarize(sum(count))
# rename column
names(raw_biotaxa_f4)[names(raw_biotaxa_f4)=="sum(count)"] <- "totalSpecies"


# split into animal/algae
names(checkedtaxatable)[names(checkedtaxatable)=="scientificname"] <- "genusSpecies"
raw_biotaxa_f4 <- raw_biotaxa_f4 %>%
  left_join(checkedtaxatable, by = "genusSpecies")
# add unknown phyla manually
raw_biotaxa_f4 <- raw_biotaxa_f4 %>%
  mutate(phylum1 = case_when(genusSpecies == "Cnemidocarpa sp." ~ "Chordata",
                            genusSpecies == "Polynoidae sp." ~ "Annelida",
                            genusSpecies == "Nematoflustra sp." ~ "Bryozoa",
                            genusSpecies == "Pycnogonid sp." ~ "Arthropoda",
                            genusSpecies == "Terebellidae sp." ~ "Annelida",
                            genusSpecies == "Munnid isopod" ~ "Arthropoda",
                            genusSpecies == "Oraderea sp." ~ "Arthropoda",
                            genusSpecies == "Flabelligera cf. mundata" ~ "Annelida",
                            genusSpecies == "Flabelligera sp." ~ "Annelida",
                            genusSpecies == "Holothuroidea sp." ~ "Echinodermata",
                            genusSpecies == "Nymphon sp." ~ "Arthropoda",
                            genusSpecies == "Heterocucumis sp." ~ "Echinodermata",
                            genusSpecies == "Arctarid sp." ~ "Arthropoda",
                            genusSpecies == "Hymenocladiopsis sp." ~ "Rhodophyta",
                            genusSpecies == "Camptoplites sp." ~ "Bryozoa",
                            genusSpecies == "Unknown amphipod 1 sp." ~ "Arthropoda",
                            genusSpecies == "Acodontaster sp." ~ "Echinodermata",
                            genusSpecies == "Benthic diatoms" ~ "Ochrophyta",
                            genusSpecies == "Lysianassid sp." ~ "Arthropoda",
                            genusSpecies == "Eusirus sp." ~ "Arthropoda",
                            genusSpecies == "Octopus/Stubby Squid sp." ~ "Mollusca",
                            genusSpecies == "Diatom benthic" ~ "Ochrophyta",
                            genusSpecies == "Benthic Diatoms sp." ~ "Ochrophyta",
                            genusSpecies == "Diatoms sp." ~ "Ochrophyta",
                            genusSpecies == "Perknaster sp." ~ "Echinodermata",
                            genusSpecies == "Oraderea cf. bidentata" ~ "Arthropoda",
                            genusSpecies == "Liothyrella sp." ~ "Brachiopoda",
                            genusSpecies == "Serolis sp." ~ "Arthropoda"
                            ))
raw_biotaxa_f4 <- raw_biotaxa_f4 %>% 
  mutate(phylum = coalesce(phylum, phylum1)) 


# rough plot of species by site
ggplot(data = raw_biotaxa_f3, aes(x = SiteID, y = genusSpecies)) +
  geom_count() 

# heatmap
ggplot(data = raw_biotaxa_f4, aes(x = SiteID, y = genusSpecies)) +
  geom_tile(aes(fill = totalSpecies)) +
  facet_grid(.~phylum) +
  scale_fill_viridis(option = "B", begin = 0.3, end = 1) +
  theme_dark() +
  theme(panel.background = element_rect(fill = 'black')) 
  
str(raw_biotaxa_f4)

###################################################################################
# IDENTIFY TARGET SPECIES                                                         #
###################################################################################


####################### species found at all sites

Site1 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "01")
Site2 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "02")
Site3 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "03")
Site4 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "04")
Site5 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "05")
Site6 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "06")
Site7 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "07")
Site8 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "08")
Site9 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "09")
Site10 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "10")
Site11 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "11")
Site12 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "12")
Site13 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "13")
Site14 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "14")
Site15 <- raw_biotaxa_2 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "15")

Reduce(intersect, list(Site1$genusSpecies, Site2$genusSpecies, Site3$genusSpecies,
                       Site4$genusSpecies, Site5$genusSpecies, Site6$genusSpecies,
                       Site7$genusSpecies, Site8$genusSpecies, Site9$genusSpecies,
                       Site10$genusSpecies, Site11$genusSpecies, Site12$genusSpecies,
                       Site13$genusSpecies, Site14$genusSpecies, Site15$genusSpecies))

# Total number of each species replicates per site
species_replicates <- raw_biotaxa_f4

reps <- 3
sites <- 8


reptable <- species_replicates %>%
  filter(totalSpecies >= reps)

sitestable <- reptable %>% 
  group_by(genusSpecies) %>% 
  mutate(count = n())
sitestable <- sitestable %>%
  filter(count >= sites)

# core as defined by at least 3 reps, at at least 8 sites
core_species <- unique(sitestable$genusSpecies)
core_table <- melt(core_species)

#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#