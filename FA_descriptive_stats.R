#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                                                                ##
# Antarctica Algae Descriptive Stats                                             ##
# Data are current as of 2021-04-26                                              ##
# Data source: Antarctic Gradients 2019                                          ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2021-05-14                                                        ##
#                                                                                ##
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# SUMMARY:

# Script to generate descriptive statistical summaries of fatty acid data collected
# from Antarctic Gradients 2019 project. 


# Required Files (check that script is loading latest version):
# dummy_sample_values.csv
# Gradients19_FA_Concs.csv
# Gradients_FA_Concs_Insight.csv
# Whippo_FA_extraction_log.csv


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
# 2021-05-14 Updated metadata, import and treatment of Insight data

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

# Browser data
Gradients19_FA_Concs <- read_csv("Data/Biomarkers/FattyAcids/Gradients19_FA_Concs.csv", 
                                                 col_types = cols(Conc = col_double(), 
                                                 Date.anal = col_character(), Notes = col_character()))

# Insight Data
Insight_concs <- read_csv("Data/Biomarkers/FattyAcids/Gradients_FA_Concs_Insight.csv", 
                               col_types = cols(Area = col_number(), 
                                                Conc = col_number(),
                                                          Date_Insight_Anal = col_character(), 
                                                          Number_ID = col_character()))
Insight_concs <- Insight_concs %>%
  mutate(across(where(is.character), ~na_if(., "----")))

Whippo_FA_extraction_log <- read_csv("~/Dropbox/OSF/Fatty Acid Extractions/Whippo_FA_extraction_log.csv")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# MANIPULATE DATA                                                              ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

############### Dummy Data

# convert data to long form
dummy_trans <- t(dummy_data)
dummy_trans <- as.data.frame(t(dummy_data[,-1]))
colnames(dummy_trans) <- dummy_data$X1
dummy_trans <- as.numeric(dummy_trans)
dummy_nums <- mutate_all(dummy_trans, function(x) as.numeric(as.character(x)))
dummy_nums <- dummy_nums %>% 
  mutate(
    across(everything(), ~replace_na(.x, 0))
  )



#long_dummy <- dummy_trans %>%
#  pivot_longer(names_to = "")



############### Actual data BROWSER

# create species-genus ID column (BROWSER)
grad_conc <- Gradients19_FA_Concs %>%
  mutate(species = str_replace_all(Sample, c("HIGR_10F0836_0130_5252020_22" = "H. grandifolius",
                                            "MYMA_09F0778_0166_5252020_7" = "M. manginii",
                                            "LAAN_14F1207_0194_5252020_14" = "L. antarctica")))

grad_conc <- grad_conc %>% 
  drop_na(Conc)

# replace NAs with 0 in Conc
grad_conc$Conc[is.na(grad_conc$Conc)] <- 0

# extract FAnumber from sample name to join with weights
grad_conc <- grad_conc %>%
  mutate(FAnumber = substr(Sample, 14, 17))

# Get original sample weights to standardize concentrations
weight_ID <- Whippo_FA_extraction_log %>%
  select(FAnumber, boatSampleWeight, boatAfterWeighing, projectID, evapVol) %>%
  filter(projectID %in% c("Gradients2019", "Gradients2019 - B sides"))

weight_net <- weight_ID %>%
  mutate(weight = boatSampleWeight - boatAfterWeighing)

conc_weight <- grad_conc %>%
  left_join(weight_net, by = "FAnumber")

conc_weight$evapVol <- recode(conc_weight$evapVol, "nd" = "1.5")
conc_weight$evapVol <- as.numeric(conc_weight$evapVol)

final_concs <- conc_weight %>%
  mutate(stand_conc = ((Conc*1000)*evapVol)/weight)
summed_FA <- final_concs %>%
  group_by(Sample) %>%
  summarise(sum(stand_conc))
final_concs <- final_concs %>%
  left_join(summed_FA, by = "Sample")
# rename column
final_concs <- final_concs %>%
  rename(summed_FA = "sum(stand_conc)")
# calculate percent of each FA per total
final_concs$FA_percent <- final_concs$stand_conc/final_concs$summed_FA





############### Actual data INSIGHT


# create species-genus ID column (INSIGHT)
grad_conc_insight <- Insight_concs %>%
  separate(Sample, into = c("speciesAbv", "sampleID", "FAnumber", "extractDate", "GCrunNumber"), sep = "_") %>%
  mutate(species = str_replace_all(speciesAbv, c("HIGR" = "H. grandifolius",
                                             "MYMA" = "M. manginii",
                                             "LAAN" = "L. antarctica")))

grad_conc_insight <- grad_conc_insight %>% 
  drop_na(Conc)

# replace NAs with 0 in Conc
grad_conc_insight$Conc[is.na(grad_conc_insight$Conc)] <- 0


# Get original sample weights to standardize concentrations
weight_ID <- Whippo_FA_extraction_log %>%
  select(FAnumber, boatSampleWeight, boatAfterWeighing, projectID, evapVol) %>%
  filter(projectID %in% c("Gradients2019", "Gradients2019 - B sides"))

weight_net <- weight_ID %>%
  mutate(weight = boatSampleWeight - boatAfterWeighing)

conc_weight_insight <- grad_conc_insight %>%
  left_join(weight_net, by = "FAnumber")

conc_weight_insight$evapVol <- recode(conc_weight_insight$evapVol, "nd" = "1.5")
conc_weight_insight$evapVol <- as.numeric(conc_weight_insight$evapVol)

final_concs_insight <- conc_weight_insight %>%
  mutate(stand_conc = ((Conc*1000)*evapVol)/weight)
summed_FA <- final_concs_insight %>%
  group_by(FAnumber) %>%
  summarise(sum(stand_conc))
final_concs_insight <- final_concs_insight %>%
  left_join(summed_FA, by = "FAnumber")
# rename column
final_concs_insight <- final_concs_insight %>%
  rename(summed_FA = "sum(stand_conc)")
# calculate percent of each FA per total
final_concs_insight$FA_proportion <- final_concs_insight$stand_conc/final_concs_insight$summed_FA




### MDS

# dummy data
dummy_MDS <- metaMDS(dummy_nums[2:18])
plot(dummy_MDS, type = "t")
dummy_MDS_points <- dummy_MDS$points
dummy_MDS_points <- data.frame(dummy_MDS_points)
plot_data_dummy <- data.frame(dummy_MDS_points, dummy_nums[,1])
ggplot(plot_data_dummy, aes(x=MDS1, y=MDS2)) +  
  theme_minimal() +
  geom_point(size = 4)

# pivot data wide for mds
grad_conc_wide <- grad_conc %>%
  select(Name, Conc, species) %>%
  pivot_wider(names_from = Name, values_from = Conc, values_fill = 0)

metaMDS(grad_conc_wide[2:14])

### Figure  BROWSER

# dotplot of final concentration by species by FA
ggplot(filter(final_concs, stand_conc > 250), aes(x = Name, y = stand_conc, colour = species)) +
  geom_point(size = 4) +
  theme_classic() +
  scale_colour_viridis(discrete = TRUE, end = 0.9) +
  labs(x = "Fatty Acid", y = "Concentration (ng/mg)") +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.1))

# stacked barplot of percent contribution of each FA to total FA (dominant)
ggplot(filter(final_concs, stand_conc > 250 & Name != 'c19.0'), aes(x = FA_proportion, y = species, fill = Name)) +
  geom_col(position = "stack") +
  scale_fill_viridis(discrete = TRUE) +
  theme_classic() +
  labs(fill = "Fatty Acid") +
  xlab("Proportion of Total FA Content") +
  ylab("Species")

# stacked barplot of percent contribution of each FA to total FA (minimal)
ggplot(filter(final_concs, stand_conc > 25 & stand_conc < 250 & Name != 'c19.0'), aes(x = FA_percent, y = species, fill = Name)) +
  geom_col(position = "stack") +
  scale_fill_viridis(discrete = TRUE) +
  labs(fill = "Fatty Acid") +
  xlab("Proportion of Total FA Content") +
  ylab("Species")

# stacked barplot of percent contribution of each FA to total FA (trace)
ggplot(filter(final_concs, stand_conc < 25 & Name != 'c19.0'), aes(x = FA_percent, y = species, fill = Name)) +
  geom_col(position = "stack") +
  scale_fill_viridis(discrete = TRUE) +
  labs(fill = "Fatty Acid") +
  xlab("Proportion of Total FA Content") +
  ylab("Species")

# LIN, ALA, SDA, ARA, EPA
ggplot(filter(final_concs, Name %in% c("c18.2n6c", "c18.3n3", "c18.4n3", "c20.4n6", "c20.5n3")), aes(x = Name, y = stand_conc, fill = species)) +
  geom_col(position = "dodge") +
  theme_classic() +
  scale_fill_viridis(discrete = TRUE) +
  scale_x_discrete(labels=c("c18.2n6c" = "LIN (c18.2n6c)", "c18.3n3" = "ALA (c18.3n3)",
                            "c20.4n6" = "ARA (c20.4n6)", "c20.5n3" = "EPA (c20.5n3)")) +
  labs(fill = "Species") +
  xlab("Fatty Acid") +
  ylab("Standardized Concentration (ng/ul)")


### Figure  INSIGHT

# dotplot of final concentration by species by FA
ggplot(filter(final_concs_insight, stand_conc > 250), aes(x = Name, y = stand_conc, colour = species)) +
  geom_point(size = 4) +
  theme_classic() +
  scale_colour_viridis(discrete = TRUE, end = 0.9) +
  labs(x = "Fatty Acid", y = "Concentration (ng/mg)") +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.1))

# boxplot of final concentration by species by FA
ggplot(filter(final_concs_insight, stand_conc > 250), aes(x = Name, y = stand_conc, colour = species)) +
  geom_boxplot() +
  theme_classic() +
  scale_colour_viridis(discrete = TRUE, end = 0.9) +
  labs(x = "Fatty Acid", y = "Concentration (ng/mg)") +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.1))

# stacked barplot of percent contribution of each FA to total FA (dominant)
ggplot(filter(final_concs_insight, stand_conc > 250 & Name != 'c19.0'), aes(x = FA_proportion, y = sampleID, fill = Name, group = as.factor(species))) +
  geom_col(position = "stack") +
  scale_fill_viridis(discrete = TRUE) +
  theme_classic() +
  labs(fill = "Fatty Acid") +
  xlab("Proportion of Total FA Content") +
  ylab("Samples") +
  facet_wrap(~ species)

  # stacked barplot of percent contribution of each FA to total FA (minimal)
ggplot(filter(final_concs_insight, stand_conc > 25 & stand_conc < 250 & Name != 'c19.0'), aes(x = FA_percent, y = species, fill = Name)) +
  geom_col(position = "stack") +
  scale_fill_viridis(discrete = TRUE) +
  labs(fill = "Fatty Acid") +
  xlab("Proportion of Total FA Content") +
  ylab("Species")

# stacked barplot of percent contribution of each FA to total FA (trace)
ggplot(filter(final_concs_insight, stand_conc < 25 & Name != 'c19.0'), aes(x = FA_percent, y = species, fill = Name)) +
  geom_col(position = "stack") +
  scale_fill_viridis(discrete = TRUE) +
  labs(fill = "Fatty Acid") +
  xlab("Proportion of Total FA Content") +
  ylab("Species")

# LIN, ALA, SDA, ARA, EPA
ggplot(filter(final_concs_insight, Name %in% c("c18.2n6c", "c18.3n3", "c18.4n3", "c20.4n6", "c20.5n3")), aes(x = Name, y = stand_conc, fill = species)) +
  geom_boxplot(position = "dodge") +
  theme_classic() +
  scale_fill_viridis(discrete = TRUE) +
  scale_x_discrete(labels=c("c18.2n6c" = "LIN (c18.2n6c)", "c18.3n3" = "ALA (c18.3n3)",
                            "c20.4n6" = "ARA (c20.4n6)", "c20.5n3" = "EPA (c20.5n3)")) +
  labs(fill = "Species") +
  xlab("Fatty Acid") +
  ylab("Standardized Concentration (ng/ul)")

  


  ####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#

# SCRATCH PAD ####