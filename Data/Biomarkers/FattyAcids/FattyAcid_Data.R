###################################################################################
#                                                                                ##
# B-236 Fatty Acid Collections Data                                              ##
# Data are current as of 2019-05-19                                              ##
# Data source: B-236 Seaweed Gradients Cruise                                    ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2019-11-13                                                        ##
#                                                                                ##
###################################################################################

# SUMMARY:

# Script for QAQC of invertebrate collection data made from the B-236 Seaweed Gradients
# Antarctic cruise from 2019-04-12:2019-05-18 including expedition members from the 
# University of Oregon (OIMB), University of Alabama Birmingham, and Texas A&M.


# Required Files (check that script is loading latest version):
# B-236_Fatty-Acid_Collections.csv

# Associated Scripts:
# GradientsBiomarkerSamples.Rmd

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
# 2019-10-17 Created QAQC .csv export for full data set for markdown summary script
# 2019-11-14 Fixed incorrect kingdom assignments of algae

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

# NEED TO CONFIRM THIS IS THE CORRECT ICE COVER VALUE FOR SITE 14, and XX values
raw_biotaxa_1 <- raw_biotaxa_1 %>%
  mutate(IceCoverCat = if_else(is.na(IceCoverCat), 0.5, IceCoverCat))

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
checkedtaxatable <- checkedtaxatable %>%
  filter(status == "accepted")

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
         "Laevilacunana antarctica" =      "Laevilacunaria antarctica",
         "Diatom benthic" =                "Benthic diatoms",
         "Benthic Diatoms sp." =           "Benthic diatoms",
         "Diatoms sp." =                   "Benthic diatoms"
         
  )

# recheck all names
taxon_names <- unique(raw_biotaxa_2$genusSpecies)
# go back to recheck, DO A COMPLETE RERUN

####################### add complete phylogeny and save full corrected data

checkedtaxatable_2 <- checkedtaxatable 
names(checkedtaxatable_2)[names(checkedtaxatable_2)=="scientificname"] <- "genusSpecies"

raw_biotaxa_3 <- raw_biotaxa_2 %>%
  left_join(checkedtaxatable_2[c(4,14:18)], by = 'genusSpecies')

raw_biotaxa_3$phylum1 <- raw_biotaxa_3$phylum

# add unknown phyla manually
raw_biotaxa_4 <- raw_biotaxa_3 %>%
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
                             genusSpecies == "Serolis sp." ~ "Arthropoda",
                             genusSpecies == "Doris kerguelenensis" ~ "Mollusca",
                             genusSpecies == "Heterocucumis steineni" ~ "Echinodermata",
                             genusSpecies == "Rhodokrambe lanigioides" ~ "Rhodophyta",
                             genusSpecies == "Pariphimedia integricauda" ~ "Arthropoda",
                             genusSpecies == "Curdiea racovitzae" ~ "Rhodophyta",
                             genusSpecies == "Psilaster charcoti" ~ "Echinodermata",
                             genusSpecies == "Glabraster antarctica" ~ "Echinodermata",
                             genusSpecies == "Oradarea bidentata" ~ "Arthropoda",
                             genusSpecies == "Charcotia obesa" ~ "Arthropoda",
                             genusSpecies == "Austropugetia crassa" ~ "Rhodophyta",
                             genusSpecies == "Phorbas bergmontae" ~ "Porifera"
  ))
raw_biotaxa_5 <- raw_biotaxa_4 %>% 
  mutate(phylum = coalesce(phylum, phylum1)) %>%
  select( -phylum1)

# add missing kingdoms manually
raw_biotaxa_6 <- raw_biotaxa_5 %>%
  mutate(kingdom1 = case_when(phylum == "Rhodophyta" ~ "Plantae",
                             phylum == "Ochrophyta" ~ "Chromista")
                           )
raw_biotaxa_7 <- raw_biotaxa_6 %>% 
  mutate(kingdom = coalesce(kingdom, kingdom1)) %>%
  select( -kingdom1)
raw_biotaxa_7$kingdom <- raw_biotaxa_7$kingdom %>%
  replace_na("Animalia")
# add gps collection location for each sample

ant_gps <- read.csv("FinalSiteLocations.csv", colClasses=c("SiteID"="character"))
ant_gps$SiteID <- as.factor(ant_gps$SiteID)

# harmonize SiteID names

ant_gps$SiteID <- ant_gps$SiteID %>%
  recode("1" =        "01", 
         "2" =        "02",
         "3" =        "03",
         "4" =        "04",
         "5" =        "05",
         "6" =        "06",
         "7" =        "07",
         "8" =        "08",
         "9" =        "09"
  )

raw_biotaxa_8 <- raw_biotaxa_7 %>%
  left_join(ant_gps, by = 'SiteID')



# write_csv(raw_biotaxa_8, "B-236_Fatty-Acid_Collections_QAQC.csv")

####################### make species list

species_list <-raw_biotaxa_8 %>% 
  group_by(genusSpecies) %>% 
  tally()

# write_csv(species_list, "B236_fattyacid-species-list.csv")

####################### extract algal species

algaetable <- checkedtaxatable %>%
  filter(kingdom != 'Animalia')
algaelist <- unique(algaetable$valid_name)

# write.csv(algaelist, "algaetable.csv")


####################### species found at all sites

Site1 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "01")
Site2 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "02")
Site3 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "03")
Site4 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "04")
Site5 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "05")
Site6 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "06")
Site7 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "07")
Site8 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "08")
Site9 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "09")
Site10 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "10")
Site11 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "11")
Site12 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "12")
Site13 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "13")
Site14 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "14")
Site15 <- raw_biotaxa_8 %>%
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

raw_biotaxa_f1 <- raw_biotaxa_8
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
                            genusSpecies == "Serolis sp." ~ "Arthropoda",
                            genusSpecies == "Doris kerguelenensis" ~ "Mollusca",
                            genusSpecies == "Heterocucumis steineni" ~ "Echinodermata",
                            genusSpecies == "Rhodokrambe lanigioides" ~ "Rhodophyta",
                            genusSpecies == "Pariphimedia integricauda" ~ "Arthropoda",
                            genusSpecies == "Curdiea racovitzae" ~ "Rhodophyta",
                            genusSpecies == "Psilaster charcoti" ~ "Echinodermata",
                            genusSpecies == "Glabraster antarctica" ~ "Echinodermata",
                            genusSpecies == "Oradarea bidentata" ~ "Arthropoda",
                            genusSpecies == "Charcotia obesa" ~ "Arthropoda",
                            genusSpecies == "Austropugetia crassa" ~ "Rhodophyta",
                            genusSpecies == "Phorbas bergmontae" ~ "Porifera"
                            ))
raw_biotaxa_f4 <- raw_biotaxa_f4 %>% 
  mutate(phylum = coalesce(phylum, phylum1)) 

# for latitude tables
raw_biotaxa_f5 <- raw_biotaxa_f1 %>%
  group_by(LAT_DD.1, genusSpecies) %>%
  summarize(sum(count))
# rename column
names(raw_biotaxa_f5)[names(raw_biotaxa_f5)=="sum(count)"] <- "totalSpecies"
raw_biotaxa_f5$LAT_DD.1 <- as.character(raw_biotaxa_f5$LAT_DD.1)


# rough plot of species by site
ggplot(data = raw_biotaxa_f3, aes(x = SiteID, y = genusSpecies)) +
  geom_count() 

# stacked bar graph of phyla collected
ggplot(data = raw_biotaxa_f4, aes(SiteID)) +
  geom_bar(aes(fill = phylum)) +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal()

# heatmap
ggplot(data = raw_biotaxa_f4, aes(x = SiteID, y = genusSpecies)) +
  geom_tile(aes(fill = totalSpecies)) +
  scale_fill_viridis(option = "B", begin = 0.3, end = 1) +
  theme_dark() +
  theme(panel.background = element_rect(fill = 'black')) +
  geom_text(show.legend = FALSE, colour = "white", aes(label = totalSpecies)) 

# heatmap by phylum
ggplot(data = raw_biotaxa_f4, aes(x = SiteID, y = genusSpecies)) +
  geom_tile(aes(fill = totalSpecies)) +
  facet_grid(.~phylum) +
  scale_fill_viridis(option = "B", begin = 0.3, end = 1) +
  theme_dark() +
  theme(panel.background = element_rect(fill = 'black')) +
  geom_text(show.legend = FALSE, colour = "white", aes(label = totalSpecies)) 
  

str(raw_biotaxa_f4)

# heatmap by latitude
ggplot(data = raw_biotaxa_f5, aes(x = LAT_DD.1, y = genusSpecies)) +
  geom_tile(aes(fill = totalSpecies)) +
  scale_fill_viridis(option = "B", begin = 0.3, end = 1) +
  theme_dark() +
  theme(panel.background = element_rect(fill = 'black')) +
  geom_text(show.legend = FALSE, colour = "white", aes(label = totalSpecies)) 

###################################################################################
# IDENTIFY TARGET SPECIES                                                         #
###################################################################################


####################### species found at all sites

Site1 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "01")
Site2 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "02")
Site3 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "03")
Site4 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "04")
Site5 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "05")
Site6 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "06")
Site7 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "07")
Site8 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "08")
Site9 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "09")
Site10 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "10")
Site11 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "11")
Site12 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "12")
Site13 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "13")
Site14 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "14")
Site15 <- raw_biotaxa_8 %>%
  group_by(SiteID, genusSpecies) %>%
  filter(SiteID == "15")

Reduce(intersect, list(Site1$genusSpecies, Site2$genusSpecies, Site3$genusSpecies,
                       Site4$genusSpecies, Site5$genusSpecies, Site6$genusSpecies,
                       Site7$genusSpecies, Site8$genusSpecies, Site9$genusSpecies,
                       Site10$genusSpecies, Site11$genusSpecies, Site12$genusSpecies,
                       Site13$genusSpecies, Site14$genusSpecies, Site15$genusSpecies))

# Total number of each species replicates per site
species_replicates <- raw_biotaxa_f4

# change these numbers and rerun to make table below
reps <- 2
sites <- 14

reptable <- species_replicates %>%
  filter(totalSpecies >= reps)

sitestable <- reptable %>% 
  group_by(genusSpecies) %>% 
  mutate(count = n())
sitestable <- sitestable %>%
  filter(count >= sites)

# core as defined by at least 3 reps, at at least 8 sites
core_species <- unique(sitestable$genusSpecies)
core_table <- melt(sort(core_species))

core_3_6 <- core_table 
core_3_8 <- core_table
core_3_10 <- core_table
core_3_14 <- core_table
core_2_10 <- core_table
core_2_14 <- core_table

core_3_6$reps <- 3
core_3_6$sites <- 6
core_3_8$reps <- 3
core_3_8$sites <- 8
core_3_10$reps <- 3
core_3_10$sites <- 10
core_3_14$reps <- 3
core_3_14$sites <- 14
core_2_10$reps <- 2
core_2_10$sites <- 10
core_2_14$reps <- 2
core_2_14$sites <- 14

core_siterep <- bind_rows(core_3_6, core_3_8, core_3_10, core_3_14,
                          core_2_10, core_2_14)

# write.csv(core_siterep, "core_species_siterep.csv")
#####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#


######### SCRATCH PAD

invert_heat_lat <- invert_heat %>%
  group_by(LAT_DD.1, genusSpecies) %>%
  summarize(sum(count))
