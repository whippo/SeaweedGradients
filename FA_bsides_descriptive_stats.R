#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                                                                ##
# Antarctica Algae FA Bsides SI Analyses                                         ##
# Script Created 2023-03-21                                                      ##
# Data source: Antarctic Gradients 2019                                          ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2023-04-01                                                        ##
#                                                                                ##
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# SUMMARY:

# Script to generate descriptive statistical summaries of bside fatty acid and
# stable isotope data from Antarctic Gradients 2019 project. 


# Required Files (check that script is loading latest version):
# Whippo_FA_extraction_log.csv
# gradients2019_bsides_FASI_QAQC.csv


# Associated Scripts:
# FA_bsides_SI_data_pipeline.R

# TO DO 
# 2. Run PERMANOVA and PCA for SI and FA separately
#     - samples in common
#     - FA with additional samples
#     - 'All' samples
# 3. Fill in all values in paper table for FA and SI


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# TABLE OF CONTENTS                                                            ####
#                                                                                 +
# RECENT CHANGES TO SCRIPT                                                        +
# LOAD PACKAGES                                                                   +
# READ IN AND PREPARE DATA                                                        +
# DATA SUMMARY                                                                    +
#                                                                                 +
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# RECENT CHANGES TO SCRIPT                                                     ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 2023-03-21 Script created from FA_descriptive_stats.R
# 2023-03-25 had to remove all c19 standard from samples

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LOAD PACKAGES                                                                ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(tidyverse) # data cleanup
library(vegan) # 'community' analyses
library(viridis) # color palette
library(psych) # pairs panel
library(ggfortify) # PCA visualizations
library(stringi) # order FA's in columns
library(factoextra) # clustering dendrogram
library(ggpubr)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LOAD FUNCTIONS                                                               ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# fuction for "%notin%
`%notin%` <- Negate(`%in%`)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# READ IN AND PREPARE DATA                                                     ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Read in all core species from data pipeline and remove duplicated 'all' FA
FASI_QAQC <- read_csv("Data/Biomarkers/FattyAcids/gradients2019_bsides_FASI_QAQC.csv")
all_species <- FASI_QAQC %>%
  select(!`19:0`) %>%
  select_if(~ !is.numeric(.) || sum(., na.rm = TRUE) != 0)

# Create simplified long dataset for analysis, remove non-overlapping samples
long_species <- all_species %>%
  select(ProjID, siteName, revisedSpecies, `Ice cover (NIC-Midpoint-Annual)`,
         targetFA, `CN ratio`:`24:1w9`) %>%
  filter(targetFA == "standards", !is.na(`CN ratio`)) %>%
  pivot_longer(cols = `CN ratio`:`24:1w9`, names_to = 'marker', values_to = 'value')

# create wide dataset, remove non-overlapping samples
overlap_species <- all_species %>%
  filter(targetFA == "standards", !is.na(`CN ratio`)) %>%
  select_if(~ !is.numeric(.) || sum(., na.rm = TRUE) != 0)

# SI only wide dataset
SI_wide <- all_species %>%
  filter(!is.na(`CN ratio`))

# FA only wide dataset
FA_wide <- all_species %>%
  filter(!is.na(`8:0`))




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# DATA SUMMARY                                                                 ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

###### BIOMARKER VALUES FOR TABLE

# calc mean and sd of each FA for each sp
FA_means <- all_species %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  mutate(across(c(`8:0`:`24:1w9`), function(x) x*100)) 
FA_means <- FA_means %>%
  select(revisedSpecies, phylum, `8:0`:`24:1w9`) %>%
  group_by(phylum, revisedSpecies) %>%
  filter(!is.na(`8:0`)) %>%
  summarise(across(`8:0`:`24:1w9`, list(mean = mean, sd = sd))) 
FA_means <- FA_means %>%
  mutate(across(where(is.numeric), round, 3))
FA_means <- as.data.frame(t(FA_means)) 
colnames(FA_means) <- FA_means[2,]


# calc mean and sd of each SI for each sp

SI_means <- all_species %>%
  select(revisedSpecies, phylum, `CN ratio`:`d13C`) %>%
  filter(!is.na(`CN ratio`)) %>%
  group_by(phylum, revisedSpecies) %>%
  summarise(across(`CN ratio`:`d13C`, list(mean = mean, sd = sd)))
SI_means <- as.data.frame(t(SI_means)) 
colnames(SI_means) <- SI_means[2,]

marker_means <- FA_means %>%
  bind_rows(SI_means[3:8,]) 
marker_means <- marker_means %>% 
  replace(is.na(.), "-") %>%
  rownames_to_column()

# write_csv(marker_means, "marker_means.csv")


###### OVERLAPPING SAMPLES no diatoms



### PERMANOVA 

# algal FA for adonis
marker_only <- overlap_species %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`CN ratio`:`24:1w9`) 

adonis2(abs(marker_only) ~ revisedSpecies, data = filter(overlap_species, revisedSpecies != "Benthic diatoms"), method = 'bray', na.rm = TRUE)

# run PCA
PCA_results <-  rda(marker_only, scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
ev <- PCA_results$CA$eig
ev>mean(ev)
# proportion explained
barplot(ev, main="eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(filter(overlap_species, revisedSpecies != "Benthic diatoms")), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make final ggplot figure
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_text(data = vscores, aes(x = -PC1, y = PC2, label = rownames(vscores)), col = 'red') +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = -PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  geom_point(aes(x = -PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies,
                 shape = phylum), size = 4) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
     y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))



##### WITHOUT DIATOMS

### PERMANOVA 

# algal FA for adonis
marker_only <- overlap_species %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`CN ratio`:`24:1w9`) 

adonis2(abs(marker_only) ~ revisedSpecies, data = filter(overlap_species, revisedSpecies != "Benthic diatoms"), method = 'bray', na.rm = TRUE)

# run PCA
PCA_results <-  rda(marker_only, scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
ev <- PCA_results$CA$eig
ev>mean(ev)
# proportion explained
barplot(ev, main="eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(filter(overlap_species, revisedSpecies != "Benthic diatoms")), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make final ggplot figure
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_text(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), col = 'red') +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  geom_point(aes(x = PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies,
                 shape = phylum), size = 4) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))






###### SI VALUES ONLY

### PERMANOVA 

# algal SI for adonis
SI_only <- SI_wide %>%
  select(`CN ratio`:d13C)

adonis2(abs(SI_only) ~ revisedSpecies, data = SI_wide, method = 'bray', na.rm = TRUE)

# run PCA
PCA_results <-  rda(SI_wide[,c(21:23)], scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
ev <- PCA_results$CA$eig
ev>mean(ev)
# proportion explained
barplot(ev, main="eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(SI_wide), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make final ggplot figure (Points scaled by 1.5)
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_text(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), col = 'red') +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  geom_point(aes(x = PC1*1.5, y = PC2*1.5, fill = revisedSpecies, color = revisedSpecies,
                 shape = phylum), size = 4) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))



###### SI VALUES ONLY no diatoms

### PERMANOVA 

# algal SI for adonis
SI_only <- SI_wide %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`CN ratio`:d13C) 

adonis2(abs(SI_only) ~ revisedSpecies, data = filter(SI_wide, revisedSpecies != "Benthic diatoms"), method = 'bray', na.rm = TRUE)

# run PCA
PCA_results <-  rda(SI_wide[,c(21:23)], scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
ev <- PCA_results$CA$eig
ev>mean(ev)
# proportion explained
barplot(ev, main="eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(filter(SI_wide, revisedSpecies != "Benthic diatoms")), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make final ggplot figure (Points scaled by 1.5)
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_text(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), col = 'red') +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  geom_point(aes(x = PC1*1.5, y = PC2*1.5, fill = revisedSpecies, color = revisedSpecies,
                 shape = phylum), size = 4) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))







###### FA VALUES ONLY



### PERMANOVA 

# algal FA for adonis
FA_only <- FA_wide %>%
  select(`8:0`:`24:1w9`) 

adonis2(abs(FA_only) ~ revisedSpecies, data = FA_wide, method = 'bray', na.rm = TRUE)

# PCA


# run PCA
PCA_results <-  rda(FA_wide[,c(24:62)], scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
ev <- PCA_results$CA$eig
ev>mean(ev)
# proportion explained
barplot(ev, main="eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(all_species), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make final ggplot figure
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_text(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), col = 'red') +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  geom_point(aes(x = PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies,
                 shape = phylum), size = 4) +
 #   geom_text(
#      aes(x = PC1, y = PC2),
#      label=uscores1$revisedSpecies,
#      check_overlap=T
#    ) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))



###### FA VALUES ONLY NO DIATOMS



### PERMANOVA 

# algal FA for adonis
FA_only <- FA_wide %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`8:0`:`24:1w9`) 

adonis2(abs(FA_only) ~ revisedSpecies, data = filter(FA_wide, revisedSpecies != "Benthic diatoms"), method = 'bray', na.rm = TRUE)

# PCA


# run PCA
PCA_results <-  rda(FA_only, scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
ev <- PCA_results$CA$eig
ev>mean(ev)
# proportion explained
barplot(ev, main="eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(filter(FA_wide, revisedSpecies != "Benthic diatoms")), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make final ggplot figure
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_text(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), col = 'red') +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  geom_point(aes(x = PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies,
                 shape = phylum), size = 4) +
  #   geom_text(
  #      aes(x = PC1, y = PC2),
  #      label=uscores1$revisedSpecies,
  #      check_overlap=T
  #    ) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))



###### FA VALUES ONLY NO DIATOMS REDUCED FA



### PERMANOVA 

# algal FA for adonis
FA_only <- FA_wide %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`20:2w6`:`20:4w6`, `16:0`, `18:3w3`) 

adonis2(abs(FA_only) ~ revisedSpecies, data = filter(FA_wide, revisedSpecies != "Benthic diatoms"), method = 'bray', na.rm = TRUE)

# PCA


# run PCA
PCA_results <-  rda(FA_only, scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
ev <- PCA_results$CA$eig
ev>mean(ev)
# proportion explained
barplot(ev, main="eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(filter(FA_wide, revisedSpecies != "Benthic diatoms")), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make final ggplot figure
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_text(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), col = 'red') +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  geom_point(aes(x = PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies,
                 shape = phylum), size = 4) +
  #   geom_text(
  #      aes(x = PC1, y = PC2),
  #      label=uscores1$revisedSpecies,
  #      check_overlap=T
  #    ) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))



###### FA VALUES ONLY NO DIATOMS REDUCED FA



### PERMANOVA 

# algal FA for adonis
FA_only <- FA_wide %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`20:5w3`:`20:4w6`, `18:2w6c`, `18:3w3`) 

adonis2(abs(FA_only) ~ revisedSpecies, data = filter(FA_wide, revisedSpecies != "Benthic diatoms"), method = 'bray', na.rm = TRUE)

# PCA


# run PCA
PCA_results <-  rda(FA_only, scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
ev <- PCA_results$CA$eig
ev>mean(ev)
# proportion explained
barplot(ev, main="eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(filter(FA_wide, revisedSpecies != "Benthic diatoms")), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make final ggplot figure
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_text(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), col = 'red') +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  geom_point(aes(x = PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies,
                 shape = phylum), size = 4) +
  #   geom_text(
  #      aes(x = PC1, y = PC2),
  #      label=uscores1$revisedSpecies,
  #      check_overlap=T
  #    ) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))



###### FA VALUES ONLY NO DIATOMS REDUCED FA



### PERMANOVA 

# algal FA for adonis
FA_only <- FA_wide %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`16:0`, `16:1w7c`, `18:3w3`, `18:1w9c`, `20:4w6`,`20:5w3`) 

adonis2(abs(FA_only) ~ revisedSpecies, data = filter(FA_wide, revisedSpecies != "Benthic diatoms"), method = 'bray', na.rm = TRUE)

# PCA


# run PCA
PCA_results <-  rda(FA_only, scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
ev <- PCA_results$CA$eig
ev>mean(ev)
# proportion explained
barplot(ev, main="eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(filter(FA_wide, revisedSpecies != "Benthic diatoms")), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make final ggplot figure
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_text(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), col = 'red') +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  geom_point(aes(x = PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies,
                 shape = phylum), size = 4) +
  #   geom_text(
  #      aes(x = PC1, y = PC2),
  #      label=uscores1$revisedSpecies,
  #      check_overlap=T
  #    ) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))



###### FA VALUES and SI NO DIATOMS REDUCED FA



### PERMANOVA 

# algal FA for adonis
FA_only <- overlap_species %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`CN ratio`, `16:0`, `20:2w6`, `18:3w3`, `20:4w6`, `18:1w9c`, `20:3w3`, `20:2w6`) 

adonis2(abs(FA_only) ~ revisedSpecies, data = filter(FA_wide, revisedSpecies != "Benthic diatoms"), method = 'bray', na.rm = TRUE)

# PCA


# run PCA
PCA_results <-  rda(FA_only, scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
ev <- PCA_results$CA$eig
ev>mean(ev)
# proportion explained
barplot(ev, main="eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(filter(overlap_species, revisedSpecies != "Benthic diatoms")), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make final ggplot figure
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_text(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), col = 'red') +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  geom_point(aes(x = PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies,
                 shape = phylum), size = 4) +
  #   geom_text(
  #      aes(x = PC1, y = PC2),
  #      label=uscores1$revisedSpecies,
  #      check_overlap=T
  #    ) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))




###### TESTING MINIMAL FA BY LOOKING AT VECTORS IN PCA (`16:0`, `18:3w3`, `20:4w6`, `22:5w6`)

### PERMANOVA 

# algal min FA for adonis
minFA_only <- overlap_species %>%
  select(`16:0`, `18:3w3`, `20:4w6`, `22:5w6`) 

adonis2(abs(minFA_only) ~ revisedSpecies, data = overlap_species, method = 'bray', na.rm = TRUE)


# run PCA
PCA_results <-  rda(minFA_only, scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
ev <- PCA_results$CA$eig
ev>mean(ev)
# proportion explained
barplot(ev, main="eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(overlap_species), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make final ggplot figure (scale on points changed to highlight spp differences)
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_text(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), col = 'red') +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  geom_point(aes(x = PC1*2, y = PC2*2, fill = revisedSpecies, color = revisedSpecies,
                 shape = phylum), size = 4) +
 #   geom_text(
#      aes(x = PC1, y = PC2),
#      label=uscores1$revisedSpecies,
#      check_overlap=T
#    ) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))





###### TESTING MINIMAL FA BY LOOKING AT VECTORS IN PCA (`16:0`, `18:3w3`, `20:4w6`, `22:5w6`)

### PERMANOVA 

# algal min FA for adonis
minFA_only <- overlap_species %>%
  select(`16:0`, `18:2w6c`, `16:1w7c`, `17:0`) 

adonis2(abs(minFA_only) ~ revisedSpecies, data = overlap_species, method = 'bray', na.rm = TRUE)


# run PCA
PCA_results <-  rda(minFA_only, scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
ev <- PCA_results$CA$eig
ev>mean(ev)
# proportion explained
barplot(ev, main="eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(overlap_species), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make final ggplot figure (scale on points changed to highlight spp differences)
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_text(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), col = 'red') +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  geom_point(aes(x = PC1*2, y = PC2*2, fill = revisedSpecies, color = revisedSpecies,
                 shape = phylum), size = 4) +
  #   geom_text(
  #      aes(x = PC1, y = PC2),
  #      label=uscores1$revisedSpecies,
  #      check_overlap=T
  #    ) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))



###### TESTING MINIMAL FA BY LOOKING AT VECTORS IN PCA (Just EFA)

### PERMANOVA 

# algal min FA for adonis
minFA_only <- overlap_species %>%
  select(`18:2w6c`, `18:3w3`, `20:4w6`, `20:5w3`, `22:6w3`) 

adonis2(abs(minFA_only) ~ revisedSpecies, data = overlap_species, method = 'bray', na.rm = TRUE)


# run PCA
PCA_results <-  rda(minFA_only, scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
ev <- PCA_results$CA$eig
ev>mean(ev)
# proportion explained
barplot(ev, main="eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(overlap_species), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make final ggplot figure (scale on points changed to highlight spp differences)
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_text(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), col = 'red') +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  geom_point(aes(x = PC1*2, y = PC2*2, fill = revisedSpecies, color = revisedSpecies,
                 shape = phylum), size = 4) +
  #   geom_text(
  #      aes(x = PC1, y = PC2),
  #      label=uscores1$revisedSpecies,
  #      check_overlap=T
  #    ) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))




###### TESTING MINIMAL FA BY LOOKING AT SIMPER VALUES NO DIATOMS

### PERMANOVA 

# algal min FA for adonis
minFA_only <- FA_wide %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`18:1w9c`, `18:3w3`, `20:4w6`, `20:5w3`, `16:0`, `16:1w7c`) 

adonis2(abs(minFA_only) ~ revisedSpecies, data = overlap_species, method = 'bray', na.rm = TRUE)


# run PCA
PCA_results <-  rda(minFA_only, scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
ev <- PCA_results$CA$eig
ev>mean(ev)
# proportion explained
barplot(ev, main="eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(filter(FA_wide, revisedSpecies != "Benthic diatoms")), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make final ggplot figure (scale on points changed to highlight spp differences)
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_text(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), col = 'red') +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  geom_point(aes(x = PC1*2, y = PC2*2, fill = revisedSpecies, color = revisedSpecies,
                 shape = phylum), size = 4) +
  #   geom_text(
  #      aes(x = PC1, y = PC2),
  #      label=uscores1$revisedSpecies,
  #      check_overlap=T
  #    ) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))






###### TESTING MINIMAL FA BY LOOKING AT VECTORS IN PCA (`16:0`, `18:3w3`, `20:4w6`, `22:5w6`)
###### WITH SI ADDED

### PERMANOVA 

# algal min FA for adonis
minFA_only <- overlap_species %>%
  select(`CN ratio`, `d15N`, `d13C`, `16:0`, `16:1w7c`, `20:5w3`) 

adonis2(abs(minFA_only) ~ revisedSpecies, data = overlap_species, method = 'bray', na.rm = TRUE)


# run PCA
PCA_results <-  rda(minFA_only, scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
ev <- PCA_results$CA$eig
ev>mean(ev)
# proportion explained
barplot(ev, main="eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(overlap_species), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make final ggplot figure (scale on points changed to highlight spp differences)
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_text(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), col = 'red') +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  geom_point(aes(x = PC1*2, y = PC2*2, fill = revisedSpecies, color = revisedSpecies,
                 shape = phylum), size = 4) +
  #   geom_text(
  #      aes(x = PC1, y = PC2),
  #      label=uscores1$revisedSpecies,
  #      check_overlap=T
  #    ) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))



###### TESTING REDS ONLY ON ALL FA SAMPLES

### PERMANOVA 

# algal red FA for adonis
redFA_only <- FA_wide %>%
  filter(phylum == "Rhodophyta") %>%
  select(`8:0`:`24:1w9`)

adonis2(abs(redFA_only) ~ revisedSpecies, data = filter(FA_wide, phylum == "Rhodophyta"), method = 'bray', na.rm = TRUE)


# run PCA
PCA_results <-  rda(redFA_only, scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
ev <- PCA_results$CA$eig
ev>mean(ev)
# proportion explained
barplot(ev, main="eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(filter(FA_wide, phylum == "Rhodophyta")), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make final ggplot figure (scale on points changed to highlight spp differences)
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_text(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), col = 'red') +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  geom_point(aes(x = PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies,
                 shape = phylum), size = 4) +
  #   geom_text(
  #      aes(x = PC1, y = PC2),
  #      label=uscores1$revisedSpecies,
  #      check_overlap=T
  #    ) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))



###### TESTING REDS ONLY ON ALL FA SAMPLES minimum data

### PERMANOVA 

# algal red FA for adonis
redminFA_only <- FA_wide %>%
  filter(phylum == "Rhodophyta") %>%
  select(`16:0`, `20:3w3`, `22:3w6`, `22:5w3`, `13:0`, `18:1w7c`)

adonis2(abs(redminFA_only) ~ revisedSpecies, data = filter(FA_wide, phylum == "Rhodophyta"), method = 'bray', na.rm = TRUE)


# run PCA
PCA_results <-  rda(redminFA_only, scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
ev <- PCA_results$CA$eig
ev>mean(ev)
# proportion explained
barplot(ev, main="eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(filter(FA_wide, phylum == "Rhodophyta")), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make final ggplot figure (scale on points changed to highlight spp differences)
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_text(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), col = 'red') +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  geom_point(aes(x = PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies,
                 shape = phylum), size = 4) +
  #   geom_text(
  #      aes(x = PC1, y = PC2),
  #      label=uscores1$revisedSpecies,
  #      check_overlap=T
  #    ) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))


###### TESTING REDS ONLY ON ALL FA SAMPLES minimum data and SI

### PERMANOVA 

# algal red FA for adonis
redminFASI_only <- overlap_species %>%
  filter(phylum == "Rhodophyta") %>%
  select(`CN ratio`, `d15N`, `d13C`, `16:0`, `20:3w3`, `22:3w6`, `22:5w3`, `13:0`, `18:1w7c`)

adonis2(abs(redminFASI_only) ~ revisedSpecies, data = filter(overlap_species, phylum == "Rhodophyta"), method = 'bray', na.rm = TRUE)


# run PCA
PCA_results <-  rda(redminFASI_only, scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
ev <- PCA_results$CA$eig
ev>mean(ev)
# proportion explained
barplot(ev, main="eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(filter(overlap_species, phylum == "Rhodophyta")), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make final ggplot figure (scale on points changed to highlight spp differences)
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_text(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), col = 'red') +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  geom_point(aes(x = PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies,
                 shape = phylum), size = 4) +
  #   geom_text(
  #      aes(x = PC1, y = PC2),
  #      label=uscores1$revisedSpecies,
  #      check_overlap=T
  #    ) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))



###### TESTING REDS ONLY ON ALL SI SAMPLES 

### PERMANOVA 

# algal red SI for adonis
redminSI_only <- overlap_species %>%
  filter(phylum == "Rhodophyta") %>%
  select(`CN ratio`:`d13C`)

adonis2(abs(redminSI_only) ~ revisedSpecies, data = filter(overlap_species, phylum == "Rhodophyta"), method = 'bray', na.rm = TRUE)


# run PCA
PCA_results <-  rda(redminSI_only, scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
ev <- PCA_results$CA$eig
ev>mean(ev)
# proportion explained
barplot(ev, main="eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(filter(overlap_species, phylum == "Rhodophyta")), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make final ggplot figure (scale on points changed to highlight spp differences)
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_text(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), col = 'red') +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  geom_point(aes(x = PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies,
                 shape = phylum), size = 4) +
  #   geom_text(
  #      aes(x = PC1, y = PC2),
  #      label=uscores1$revisedSpecies,
  #      check_overlap=T
  #    ) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))


# Proportional stacked barplot of FA composition for each species

# calc mean of each FA for each sp
FA_quartile <- long_species %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  filter(marker %notin% c("d13C", "d15N", "CN ratio")) %>%
  filter(value != 0)
summary(FA_quartile$value)

FA_means <- FA_wide %>%
  select(revisedSpecies, phylum, `8:0`:`24:1w9`) %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  group_by(phylum, revisedSpecies) %>%
  summarise(across(`8:0`:`24:1w9`, mean)) %>%
  pivot_longer(cols = `8:0`:`24:1w9`, names_to = 'Fatty Acid', values_to = 'value') %>%
  filter(value >= 0.05) %>% # 1st quartile of all FA values in dataset
  mutate(phylum = case_when(phylum == "Chlorophyta" ~ "",
                            phylum == "Ochrophyta" ~ "Ochrophyta",
                            phylum == "Rhodophyta" ~ "Rhodophyta"))


FA_means %>%
  ggplot() +
  geom_col(aes(x = revisedSpecies, y = value, fill = `Fatty Acid`), position = "fill") +
  scale_fill_viridis(discrete = TRUE, option = 6) +
  theme_bw() +
  labs(x = "Species", y = "Mean Proportional Composition") +
  guides(fill = guide_legend(title = "Fatty Acids")) +
  facet_grid(cols = vars(phylum), scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.25)) +
  theme(axis.title.x = element_text(vjust = -1)) +
  theme(axis.title.y = element_text(vjust = 2.5))



# MDS of community
FA_mds <- metaMDS(FA_only)
FA_mds_points <- FA_mds$points
FA_mds_points <- data.frame(FA_mds_points)
plot_data_tax <- data.frame(FA_mds_points, FA_mat[,1:3])
library(plyr)
chulls_tax <- ddply(plot_data_tax, .(Habitat), function(df) df[chull(df$MDS1, df$MDS2), ])
detach(package:plyr)
plot(FA_mds)

simper_FA <- FA_wide %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`8:0`:`24:1w9`) 

full_algal_simper <- simper(simper_FA)
summary(full_algal_simper)


# cluster analysis of FA only

no_diatoms_meanFA <- FA_wide %>%
  select(revisedSpecies, phylum, `8:0`:`24:1w9`) %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  group_by(phylum, revisedSpecies) %>%
  summarise(across(`8:0`:`24:1w9`, mean)) 


Alg_dist <- vegdist(no_diatoms_meanFA[,3:41])
Alg_clust <- hclust(Alg_dist, method="ward.D2")
Alg_order <- data.frame(no_diatoms_meanFA$revisedSpecies, Alg_clust$order)
Alg_order <- Alg_order[order(Alg_order$Alg_clust.order), ]

plot(Alg_clust, las = 1, 
     main="Cluster diagram of algae", 
     xlab="Sample", 
     ylab="Euclidean distance",
     label = Alg_order$FA_wide.revisedSpecies)

Alg_clust$labels <- no_diatoms_meanFA$revisedSpecies
clust_col <- viridis(4, alpha = 1, begin = 0.2, end = 0.8, direction = 1, option = "C")

fviz_dend(x = Alg_clust, cex = 0.8, lwd = 0.8, k = 4, 
          k_colors = clust_col,
          rect = TRUE, 
          rect_border = "jco", 
          rect_fill = TRUE,
          type = "circular",
          ggtheme = theme_bw())

# cluster analysis of SI only

no_diatoms_meanSI <- SI_wide %>%
  select(revisedSpecies, phylum, `CN ratio`:`d13C`) %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  group_by(phylum, revisedSpecies) %>%
  summarise(across(`CN ratio`:`d13C`, mean)) 

Alg_dist <- vegdist(abs(no_diatoms_meanSI[,3:5]))
Alg_clust <- hclust(Alg_dist, method="ward.D2")
Alg_order <- data.frame(SI_wide$revisedSpecies, Alg_clust$order)
Alg_order <- Alg_order[order(Alg_order$Alg_clust.order), ]

plot(Alg_clust, las = 1, 
     main="Cluster diagram of algae", 
     xlab="Sample", 
     ylab="Euclidean distance",
     label = Alg_order$FA_wide.revisedSpecies)

Alg_clust$labels <- no_diatoms_meanSI$revisedSpecies
clust_col <- viridis(4, alpha = 1, begin = 0.2, end = 0.8, direction = 1, option = "C")

fviz_dend(x = Alg_clust, cex = 0.8, lwd = 0.8, k = 4, 
          k_colors = clust_col,
          rect = TRUE, 
          rect_border = "jco", 
          rect_fill = TRUE,
          type = "circular",
          ggtheme = theme_bw())



##### What are the diffs if any between published and non-published FA per phylum?


  
  # calc mean and sd of each FA for each sp
  FA_means <- all_species %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
    mutate(published = case_when(revisedSpecies == "Lambia antarctica" ~ "published",
                                 revisedSpecies == "Ascoseira mirabilis" ~ "published",
                                 revisedSpecies == "Desmarestia anceps" ~ "published",
                                 revisedSpecies == "Desmarestia antarctica" ~ "published",
                                 revisedSpecies == "Desmarestia menziesii" ~ "published",
                                 revisedSpecies == "Himantothallus grandifolius" ~ "published",
                                 revisedSpecies == "Adenocystis utricularis" ~ "published",
                                 revisedSpecies == "Delisea pulchra" ~ "published",
                                 revisedSpecies == "Georgiella confluens" ~ "published",
                                 revisedSpecies == "Myriogramme smithii" ~ "published",
                                 revisedSpecies == "Myriogramme manginii" ~ "published",
                                 revisedSpecies == "Pantoneura plocamioides" ~ "published",
                                 revisedSpecies == "Sarcopeltis antarctica" ~ "published",
                                 revisedSpecies == "Iridaea cordata" ~ "published",
                                 revisedSpecies == "Curdiea racovitzae" ~ "published",
                                 revisedSpecies == "Palmaria decipiens" ~ "published",
                                 revisedSpecies == "Plocamium sp." ~ "published",
                                 revisedSpecies == "Hymenocladiopsis sp." ~ "published",
                                 revisedSpecies == "Cystosphaera jacquinotii" ~ "unpublished",
                                 revisedSpecies == "Microzonia australe" ~ "unpublished",
                                 revisedSpecies == "Ballia callitricha" ~ "unpublished",
                                 revisedSpecies == "Porphyra plocamiestris" ~ "unpublished",
                                 revisedSpecies == "Paraglossum salicifolium" ~ "unpublished",
                                 revisedSpecies == "Picconiella plumosa" ~ "unpublished",
                                 revisedSpecies == "Meridionella antarctica" ~ "unpublished",
                                 revisedSpecies == "Austropugetia crassa" ~ "unpublished",
                                 revisedSpecies == "Callophyllis atrosanguinea" ~ "unpublished",
                                 revisedSpecies == "Gymnogongrus antarcticus" ~ "unpublished",
                                 revisedSpecies == "Phyllophora antarctica" ~ "unpublished",
                                 revisedSpecies == "Pachymenia orbicularis" ~ "unpublished",
                                 revisedSpecies == "Trematocarpus antarcticus" ~ "unpublished")) %>%
  mutate(across(c(`8:0`:`24:1w9`), function(x) x*100))
FA_means <- FA_means %>%
  select(revisedSpecies, phylum, published, `8:0`:`24:1w9`) %>%
  group_by(phylum, published) %>%
  filter(!is.na(`8:0`)) %>%
  summarise(across(`8:0`:`24:1w9`, list(mean = mean))) 
top_FA <- FA_means %>%
  pivot_longer(`8:0_mean`:`24:1w9_mean`, names_to = "FA", values_to = "Percent") %>%
  arrange(desc(Percent)) %>% 
  group_by(phylum, published) %>%
  slice(1:5) %>%
  mutate(rank = c(1:5)) %>%
  pivot_wider(names_from = "rank", values_from = c("FA", "Percent"))

# grand means

FA_means <- FA_means %>%
  select(revisedSpecies, phylum, `8:0`:`24:1w9`) %>%
  group_by(phylum) %>%
  filter(!is.na(`8:0`)) %>%
  summarise(across(`8:0`:`24:1w9`, list(mean = mean))) 
top_FA <- FA_means %>%
  pivot_longer(`8:0_mean`:`24:1w9_mean`, names_to = "FA", values_to = "Percent") %>%
  arrange(desc(Percent)) %>% 
  group_by(phylum) %>%
  slice(1:10) %>%
  mutate(rank = c(1:10)) %>%
  pivot_wider(names_from = "rank", values_from = c("FA", "Percent"))




###### FA VALUES ONLY NO DIATOMS, ORDER LEVEL ANALYSIS
# algal FA for adonis
FA_only <- FA_wide %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`8:0`:`24:1w9`) 
FA_tax <- FA_wide %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  mutate(order = case_when(revisedSpecies == "Lambia antarctica" ~ "Bryopsidales",
                           revisedSpecies == "Ascoseira mirabilis" ~ "Ascoseirales",
                           revisedSpecies == "Desmarestia anceps" ~ "Desmarestiales",
                           revisedSpecies == "Desmarestia antarctica" ~ "Desmarestiales",
                           revisedSpecies == "Desmarestia menziesii" ~ "Desmarestiales",
                           revisedSpecies == "Himantothallus grandifolius" ~ "Desmarestiales",
                           revisedSpecies == "Adenocystis utricularis" ~ "Ectocarpales",
                           revisedSpecies == "Cystosphaera jacquinotii" ~ "Fucales",
                           revisedSpecies == "Microzonia australe" ~ "Syringodermatales",
                           revisedSpecies == "Ballia callitricha" ~ "Balliales",
                           revisedSpecies == "Porphyra plocamiestris" ~ "Bangiales",
                           revisedSpecies == "Delisea pulchra" ~ "Bonnemaisoniales",
                           revisedSpecies == "Georgiella confluens" ~ "Ceramiales",
                           revisedSpecies == "Myriogramme smithii" ~ "Ceramiales",
                           revisedSpecies == "Myriogramme manginii" ~ "Ceramiales",
                           revisedSpecies == "Pantoneura plocamioides" ~ "Ceramiales",
                           revisedSpecies == "Paraglossum salicifolium" ~ "Ceramiales",
                           revisedSpecies == "Picconiella plumosa" ~ "Ceramiales",
                           revisedSpecies == "Meridionella antarctica" ~ "Gigartinales",
                           revisedSpecies == "Sarcopeltis antarctica" ~ "Gigartinales",
                           revisedSpecies == "Iridaea cordata" ~ "Gigartinales",
                           revisedSpecies == "Austropugetia crassa" ~ "Gigartinales",
                           revisedSpecies == "Callophyllis atrosanguinea" ~ "Gigartinales",
                           revisedSpecies == "Gymnogongrus antarcticus" ~ "Gigartinales",
                           revisedSpecies == "Phyllophora antarctica" ~ "Gigartinales",
                           revisedSpecies == "Curdiea racovitzae" ~ "Gracilariales",
                           revisedSpecies == "Pachymenia orbicularis" ~ "Halymeniales",
                           revisedSpecies == "Palmaria decipiens" ~ "Palmariales",
                           revisedSpecies == "Plocamium sp." ~ "Plocamiales",
                           revisedSpecies == "Trematocarpus antarcticus" ~ "Plocamiales",
                           revisedSpecies == "Hymenocladiopsis sp." ~ "Rhodymeniales")) %>%
  mutate(family = case_when(revisedSpecies == "Lambia antarctica" ~ "Bryopsidaceae",
                           revisedSpecies == "Ascoseira mirabilis" ~ "Ascoseiraceae",
                           revisedSpecies == "Desmarestia anceps" ~ "Desmarestiaceae",
                           revisedSpecies == "Desmarestia antarctica" ~ "Desmarestiaceae",
                           revisedSpecies == "Desmarestia menziesii" ~ "Desmarestiaceae",
                           revisedSpecies == "Himantothallus grandifolius" ~ "Desmarestiaceae",
                           revisedSpecies == "Adenocystis utricularis" ~ "Adenocystaceae",
                           revisedSpecies == "Cystosphaera jacquinotii" ~ "Seirococcaceae",
                           revisedSpecies == "Microzonia australe" ~ "Syringodermataceae",
                           revisedSpecies == "Ballia callitricha" ~ "Balliaceae",
                           revisedSpecies == "Porphyra plocamiestris" ~ "Bangiaceae",
                           revisedSpecies == "Delisea pulchra" ~ "Bonnemaisoniaceae",
                           revisedSpecies == "Georgiella confluens" ~ "Callithamniaceae",
                           revisedSpecies == "Myriogramme smithii" ~ "Delesseriaceae",
                           revisedSpecies == "Myriogramme manginii" ~ "Delesseriaceae",
                           revisedSpecies == "Pantoneura plocamioides" ~ "Delesseriaceae",
                           revisedSpecies == "Paraglossum salicifolium" ~ "Delesseriaceae",
                           revisedSpecies == "Picconiella plumosa" ~ "Rhodomelaceae",
                           revisedSpecies == "Meridionella antarctica" ~ "Cystocloniaceae",
                           revisedSpecies == "Sarcopeltis antarctica" ~ "Gigartinaceae",
                           revisedSpecies == "Iridaea cordata" ~ "Gigartinaceae",
                           revisedSpecies == "Austropugetia crassa" ~ "Kallymeniaceae",
                           revisedSpecies == "Callophyllis atrosanguinea" ~ "Kallymeniaceae",
                           revisedSpecies == "Gymnogongrus antarcticus" ~ "Phyllophoraceae",
                           revisedSpecies == "Phyllophora antarctica" ~ "Phyllophoraceae",
                           revisedSpecies == "Curdiea racovitzae" ~ "Gracilariaceae",
                           revisedSpecies == "Pachymenia orbicularis" ~ "Halymeniaceae",
                           revisedSpecies == "Palmaria decipiens" ~ "Palmariaceae",
                           revisedSpecies == "Plocamium sp." ~ "Plocamiaceae",
                           revisedSpecies == "Trematocarpus antarcticus" ~ "Sarcodiaceae",
                           revisedSpecies == "Hymenocladiopsis sp." ~ "Fryeellaceae")) 
  

### PERMANOVA 


adonis2(abs(FA_only) ~ revisedSpecies, data = filter(FA_wide, revisedSpecies != "Benthic diatoms"), method = 'bray', na.rm = TRUE)

# PCA


# run PCA
PCA_results <-  rda(FA_only, scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
ev <- PCA_results$CA$eig
ev>mean(ev)
# proportion explained
barplot(ev, main="eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(FA_tax), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make final ggplot figure
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  geom_text(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)), col = 'red') +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  geom_point(aes(x = PC1, y = PC2, fill = family, color = family,
                 shape = phylum), size = 4) +
  #   geom_text(
  #      aes(x = PC1, y = PC2),
  #      label=uscores1$revisedSpecies,
  #      check_overlap=T
  #    ) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))







######## nMDS

FA_matrix <- FA_only

# run the nMDS
FA_mds <- metaMDS(FA_matrix)
# extract the 'points' from the nMDS that you will plot in ggplot2
FA_mds_points <- FA_mds$points
# turn those plot points into a dataframe that ggplot2 can read
FA_mds_points <- data.frame(FA_mds_points)
# join your plot points with your summed species observations from each habitat type
plot_data_tax <- data.frame(FA_mds_points, FA_tax[,c(7,8,63,64)])
plot_data_tax <- plot_data_tax %>%
  rename("division" = "phylum")
# IF you want to add hulls around your points (totally optional), use this code
# first you have to load plyr (DO NOT LOAD PRIOR TO THIS. It tends to mess with tidyverse functions)
library(plyr)
# create the list of points that will connect the 'hulls' together from your nMDS point data
chulls_tax <- ddply(plot_data_tax, .(Habitat), function(df) df[chull(df$MDS1, df$MDS2), ])
# DETACH PLYR so it won't mess with anything!
detach(package:plyr)

# run the ggplot
phylum_plot <- 
  ggplot(plot_data_tax, aes(x=MDS1, y=MDS2, # note pch and color are linked to *factors*
                          color = division)) +  
  theme_classic() + # optional, I just like this theme
  geom_point(size = 2) + # set size of points to whatever you want
 # geom_polygon(data=chulls_tax, aes(x=MDS1, y=MDS2, group=Habitat), fill=NA) + # optional 'hulls' around points
  scale_color_viridis(discrete = TRUE, begin = 0.2, end = 0.9, option = "G") # my favorite color-blind and b&w friendly palette, look at the viridis package for more details 
phylum_leg <- as_ggplot(get_legend(phylum_plot))

order_plot <- 
  ggplot(plot_data_tax, aes(x=MDS1, y=MDS2, # note pch and color are linked to *factors*
                                         color = order)) +  
  theme_classic() + # optional, I just like this theme
  geom_point(size = 2) + # set size of points to whatever you want
    guides(color=guide_legend(ncol=2)) +
  # geom_polygon(data=chulls_tax, aes(x=MDS1, y=MDS2, group=Habitat), fill=NA) + # optional 'hulls' around points
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "B") # my favorite color-blind and b&w friendly palette, look at the viridis package for more details 
order_leg <- as_ggplot(get_legend(order_plot))


family_plot <- 
  ggplot(plot_data_tax, aes(x=MDS1, y=MDS2, # note pch and color are linked to *factors*
                                         color = family)) +  
  theme_classic() + # optional, I just like this theme
  geom_point(size = 2) + # set size of points to whatever you want
  # geom_polygon(data=chulls_tax, aes(x=MDS1, y=MDS2, group=Habitat), fill=NA) + # optional 'hulls' around points
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "H") # my favorite color-blind and b&w friendly palette, look at the viridis package for more details 
family_leg <- as_ggplot(get_legend(family_plot))

FigureMDS <- ggarrange(ggarrange(phylum_plot, order_plot, family_plot,
                     labels = c("A", "B", "C"),
                     ncol = 1, nrow = 3,
                     legend = "none"), 
                     ggarrange(phylum_leg, order_leg, family_leg,
                               ncol = 1, nrow = 3, align = "v",
                               legend = "none"), 
                     ncol = 2, nrow = 1, legend = "none")
FigureMDS

# 800 x 1200

annotate_figure(FigureMDS, top = text_grob("2D stress = 0.12", size = 10))

# best size: ~630x700




####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#

# SCRATCH PAD ####
algal_grad <- grad_conc_wide %>%
  select(genusSpecies, `NIC-Klein-Midpoint-Annual`, '14:0', '16:0', '16:1w7c', '18:0', '18:4w3c', '20:4w6', '20:5w3') %>%
  filter(genusSpecies %in% c("Desmarestia menziesii", "Phyllophora antarctica", "Himantothallus grandifolius"))
algal_grad_long <- algal_grad %>%
  pivot_longer(cols = `14:0`:`20:5w3`, names_to = 'FA', values_to = 'proportion')


algal_grad_long %>%
  ggplot(aes(`NIC-Klein-Midpoint-Annual`, proportion)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_wrap(FA ~ genusSpecies, scales = 'free', ncol = 3)



# MDS ALGAE


# select EPA, ARA, SDA, PAL, OLE, LIN, VAC, and dominant sats (16, 18) 
sub_wide <- grad_conc_wide %>% # fix 16:0 error in inverts
  select(genusSpecies, FAsampleName, "14:0", "16:0", "16:1w7c", "18:0", "18:1w7c", "18:3w3", "18:4w3c", "18:1w9c", "20:4w6", "20:5w3")
sub_wide_trans <- (sub_wide[,3:12])
sub_wide_trans <- bind_cols(sub_wide[1:2], sub_wide_trans)

batch_1_2_MDS <- metaMDS(sub_wide_trans[3:12], autotransform = TRUE, distance = "clark")
batch_1_2_MDS_points <- batch_1_2_MDS$points
batch_1_2_MDS_points <- data.frame(batch_1_2_MDS_points)
plot_data_batch_1_2 <- data.frame(batch_1_2_MDS_points, sub_wide[,1])

library(plyr)
# create the list of points that will connect the 'hulls' together from your nMDS point data
chulls_tax <- ddply(plot_data_batch_1_2, .(genusSpecies), function(df) df[chull(df$MDS1, df$MDS2), ])
# DETACH PLYR so it won't mess with anything!
detach(package:plyr)

# create vectors to plot over MDS
scrs <- as.data.frame(scores(batch_1_2_MDS, display = "sites"))
scrs <- cbind(scrs, genusSpecies = sub_wide_trans$genusSpecies)

vf <- envfit(batch_1_2_MDS, sub_wide_trans[3:12], perm = 999)

spp.scrs <- as.data.frame(scores(vf, display = "vectors"))
spp.scrs <- cbind(spp.scrs, FA = rownames(spp.scrs))
spp.scrs <- spp.scrs %>%
  rename(MDS1 = NMDS1,
         MDS2 = NMDS2)

# MDS FOR ITRS - Aaron

ggplot(plot_data_batch_1_2, aes(x=MDS1, y=MDS2)) +
  coord_fixed() +
  geom_point(size = 4, aes(color = genusSpecies)) + # set size of points to whatever you want
  scale_color_viridis(discrete = TRUE, end = 0.9, name = "Algal Species") + # my favorite color-blind and b&w friendly palette, look at the viridis package for more details
  # geom_segment(data = spp.scrs,
  #            aes(x = 0, xend = MDS1, y = 0, yend = MDS2),
  #            arrow = arrow(length = unit(0.25, "cm")), color = "grey") +
  # geom_text(data = spp.scrs, aes(label = FA), 
  #         size = 3) +
  # geom_polygon(data = chulls_tax,
  #            aes(x = MDS1, y = MDS2, color = genusSpecies), 
  #            fill = NA) + # optional 'hulls' around points
  theme_classic() + # optional, I just like this theme
  annotate("text", x = 1, y = 1, label = "3D Stress = 0.17") 



# MDS INVERTS

# pivot data wide for mds
grad_conc_wide_invert <- long_inverts %>%
  select(FA, genusSpecies, proportion, FAsampleName) %>%
  pivot_wider(names_from = FA, values_from = proportion, values_fill = 0)

# select EPA, ARA, PAL, OLE, LIN, VAC, and dominant sats (16, 18), reduce inverts included
sub_wide_invert <- grad_conc_wide_invert %>% # fix 16:0 error in inverts
  select(genusSpecies, FAsampleName, "14:0", "16:0", "16:1w7c", "18:0", "18:1w7c", "18:3w3", "18:1w9c", "20:4w6", "20:5w3")
sub_wide_invert <- sub_wide_invert %>%
  filter(genusSpecies %in% c("Odontaster validus", 
                             "Cnemidocarpa sp.",
                             "Nacella concinna",
                             "Sterechinus neumayeri",
                             "Prostebbingia gracilis",
                             "Gondogeneia antarctica"))
sub_wide_trans_invert <- (sub_wide_invert[,3:11])
sub_wide_trans_invert <- bind_cols(sub_wide_invert[1:2], sub_wide_trans_invert)


batch_1_2_MDS_invert <- metaMDS(sqrt(sub_wide_trans_invert[3:11]), autotransform = TRUE, distance = "manhattan")
batch_1_2_MDS_points_invert <- batch_1_2_MDS_invert$points
batch_1_2_MDS_points_invert <- data.frame(batch_1_2_MDS_points_invert)
plot_data_batch_1_2_invert <- data.frame(batch_1_2_MDS_points_invert, sub_wide_invert[,1])

library(plyr)
# create the list of points that will connect the 'hulls' together from your nMDS point data
chulls_tax_invert <- ddply(plot_data_batch_1_2_invert, .(genusSpecies), function(df) df[chull(df$MDS1, df$MDS2), ])
# DETACH PLYR so it won't mess with anything!
detach(package:plyr)

# create vectors to plot over MDS
scrs_invert <- as.data.frame(scores(batch_1_2_MDS_invert, display = "sites"))
scrs_invert <- cbind(scrs_invert, genusSpecies = sub_wide_trans_invert$genusSpecies)

vf_invert <- envfit(batch_1_2_MDS_invert, sub_wide_trans_invert[3:11], perm = 999)

spp.scrs_invert <- as.data.frame(scores(vf_invert, display = "vectors"))
spp.scrs_invert <- cbind(spp.scrs_invert, FA = rownames(spp.scrs_invert))
spp.scrs_invert <- spp.scrs_invert %>%
  rename(MDS1 = NMDS1,
         MDS2 = NMDS2)

# MDS FOR ITRS - Aaron

ggplot(plot_data_batch_1_2_invert, aes(x=MDS1, y=MDS2)) +
  coord_fixed() +
  geom_point(size = 4, aes(color = genusSpecies)) + # set size of points to whatever you want
  scale_color_viridis(discrete = TRUE, end = 0.9, name = "Invert Species") + # my favorite color-blind and b&w friendly palette, look at the viridis package for more details
  # geom_segment(data = spp.scrs_invert,
  #           aes(x = 0, xend = MDS1, y = 0, yend = MDS2),
  #          arrow = arrow(length = unit(0.25, "cm")), color = "grey") +
  # geom_text(data = spp.scrs_invert, aes(label = FA), 
  #      size = 3) +
  # geom_polygon(data = chulls_tax_invert,
  #          aes(x = MDS1, y = MDS2, color = genusSpecies), 
  #         fill = NA) + # optional 'hulls' around points
  theme_classic() + # optional, I just like this theme
  annotate("text", x = -8, y = -9, label = "3D Stress = 0.08")  

# INVERT MDS w/ ICE COVER

site_sample_invert <- all_inverts %>%
  select(SiteID, FAsampleName, `NIC-Klein-Midpoint-Annual`)

invert_cover <- sub_wide_invert %>%
  left_join(site_sample_invert, by = "FAsampleName")

plot_data_batch_1_2_invert <- data.frame(batch_1_2_MDS_points_invert, invert_cover[,c(1,13)])

library(plyr)
# create the list of points that will connect the 'hulls' together from your nMDS point data
chulls_tax_invert <- ddply(plot_data_batch_1_2_invert, .(`NIC.Klein.Midpoint.Annual`), function(df) df[chull(df$MDS1, df$MDS2), ])
# DETACH PLYR so it won't mess with anything!
detach(package:plyr)

# create vectors to plot over MDS
scrs_invert <- as.data.frame(scores(batch_1_2_MDS_invert, display = "sites"))
scrs_invert <- cbind(scrs_invert, IceCover = invert_cover$`NIC-Klein-Midpoint-Annual`)

vf_invert <- envfit(batch_1_2_MDS_invert, sub_wide_trans_invert[3:11], perm = 999)

spp.scrs_invert <- as.data.frame(scores(vf_invert, display = "vectors"))
spp.scrs_invert <- cbind(spp.scrs_invert, FA = rownames(spp.scrs_invert))
spp.scrs_invert <- spp.scrs_invert %>%
  rename(MDS1 = NMDS1,
         MDS2 = NMDS2)

# plot_data_batch_1_2_invert$NIC.Klein.Midpoint.Annual <- as.character(plot_data_batch_1_2_invert$NIC.Klein.Midpoint.Annual)


ggplot(plot_data_batch_1_2_invert, aes(x=MDS1, y=MDS2)) +
  coord_fixed() +
  geom_point(size = 4, aes(color = `NIC.Klein.Midpoint.Annual`)) + # set size of points to whatever you want
  scale_color_viridis(discrete = FALSE, end = 0.9, name = "Ice Cover") + # my favorite color-blind and b&w friendly palette, look at the viridis package for more details
  #geom_segment(data = spp.scrs_invert,
  #           aes(x = 0, xend = MDS1, y = 0, yend = MDS2),
  #          arrow = arrow(length = unit(0.25, "cm")), color = "grey") +
  # geom_text(data = spp.scrs_invert, aes(label = FA), 
  #      size = 3) +
  # geom_polygon(data = chulls_tax_invert,
  #          aes(x = MDS1, y = MDS2, color = genusSpecies), 
  #         fill = NA) + # optional 'hulls' around points
  theme_classic() + # optional, I just like this theme
  annotate("text", x = -8, y = -9, label = "3D Stress = 0.08") 


# PAIRS PANELS OF FA

# all gradients site joining ALGAE

grad_site <- long_algae %>%
  mutate(site = str_sub(sample, 6, 7))

# create site/lat, and add gradient - taken from NIC-Klein-Midpoint-Annual (email from Chuck)
site <- c('02', '03', '04', '05', '07', '08', '09', '10', '11', '12', '13', '14', '15')
lat <- c(-67.5567, -68.1758, -68.6921, -67.5488, -66.0894, -65.5131, -65.1043, -64.9002, -64.7932, -64.7720, -66.0251, -64.7793, -65.2402)
grad <- c('71.5102', '68.61224', '87.67347', '57.87755', '58.36735', '73.95918', '53.5102', '36.12245', '41.10204', '37.26531', '76.77551', '41.06122', '62.85714')

site_lat <- data.frame(site, lat, grad)

all_site <- left_join(grad_site, site_lat, by = "site")
all_site$lat <- as.factor(all_site$lat)

# pivot data wide for pairs
all_site_wide <- all_site %>%
  select(FA, genusSpecies, proportion, sample, site, lat, grad) %>%
  pivot_wider(names_from = FA, values_from = proportion, values_fill = 0)
# remove zero columns
all_site_pairs <- all_site_wide %>%
  select(where(~ any(. != 0)))
# create mean values
all_site_means <- all_site_pairs %>%
  select(-c('site', 'sample')) %>%
  group_by(genusSpecies, lat) %>%
  summarise(across(everything(), mean))
# reduce number of FAs
all_site_reduced <- all_site_pairs %>% # switch between lat/grad for gradients here
  select(genusSpecies, grad, `14:0`, `16:0`, `16:1w7c`, `16:1w5c`, `18:0`, `18:1w7c`, `18:1w9c`, `18:2w6c`,`18:3w3`, `18:4w3c`, `20:0`, `20:4w6`, `20:5w3`, `22:0`) %>%
  #filter(genusSpecies == 'Plocamium cartilagineum')
  #filter(genusSpecies == 'Desmarestia menziesii')
  #filter(genusSpecies == 'Phyllophora antarctica')
  filter(genusSpecies == 'Himantothallus grandifolius')
#filter(genusSpecies == 'Myriogramme manginii')
pairs.panels(all_site_reduced[,2:16])


# FAs of lattitude interest (>|50| Pearson correlation [moderate or stronger]):
# (as of batch 5 for FAs - 
# `14:0`, `16:0...50`, `16:1w7c`, `16:1w5c`, `18:0`, `18:1w7c`, `18:1w9c`, 
# `18:2w6c`,`18:3w3`, `18:4w3c`, `20:0`, `20:4w6`, `20:5w3`, `22:0`)
#
## P. cartilagineum
#     22:0 = 0.63
## D. menziesii
#     20:5w3 = -0.58
## P. antarctica
#     22:0 = -0.58
#     20:5w3 = -0.53
#     20:0 = 0.58
#     18:3w3 = -0.52
#     16:1w5 = 0.87
## H. grandifolius
#     20:5w3 = -0.62
#     20:0 = 0.66
#     18:4w3 = -0.74
#     18:3w3 = -0.50
#     18:2w6 = 0.52
#     18:1w7 = -0.50
#     18:0 = -0.69
#     16:1w7 = 0.76
#     16:0 = 0.51

# FAs of gradient interest (>|50| Pearson correlation [moderate or stronger]):
# (as of batch 5 for FAs - 
# `14:0`, `16:0...50`, `16:1w7c`, `16:1w5c`, `18:0`, `18:1w7c`, `18:1w9c`, 
# `18:2w6c`,`18:3w3`, `18:4w3c`, `20:0`, `20:4w6`, `20:5w3`, `22:0`)
#
## D. menziesii
#     20:5w3 = 0.58
## P. antarctica
#     18:4w3 = 0.56
#     18:2w6 = 0.50
#     18:0 = 0.56
#     16:0 = -0.63
## H. grandifolius
#     20:5w3 = 0.89
#     20:0 = -0.67
#     18:4w3 = 0.96
#     18:3w3 = 0.77
#     18:1w9 = -0.51
#     18:1w7 = 0.58
#     18:0 = 0.64
#     16:1w7 = -0.99
#     16:0 = -0.83
#     14:0 = 0.52


# cluster analysis of FAs that are related
FA_dist <- vegdist(t(sub_wide_trans[3:12]))
FA_clust <- hclust(FA_dist, method="ward.D2") 

plot(FA_clust, las = 1, 
     main="Cluster diagram of algal fatty acids", 
     xlab="Sample", 
     ylab="Euclidean distance")


# cluster analysis of Algae that are related
Alg_dist <- vegdist(sub_wide_trans[3:12])
Alg_clust <- hclust(Alg_dist, method="ward.D2")
Alg_order <- data.frame(sub_wide_trans$species, Alg_clust$order)
Alg_order <- Alg_order[order(Alg_order$Alg_clust.order), ]

plot(Alg_clust, las = 1, 
     main="Cluster diagram of algae", 
     xlab="Sample", 
     ylab="Euclidean distance",
     label = Alg_order$sub_wide_trans.species)





# dotplot of final concentration by species by FA
ggplot(filter(long_algae, FA %in% c("14:0", "16:0", "16:1w7c", "18:0", "18:1w7c", "18:3w3", "18:4w3c", "18:1w9c", "20:4w6", "20:5w3")), 
       aes(x = FA, y = proportion, colour = genusSpecies)) +
  geom_point(size = 4, position = position_jitter(width = .1)) +
  geom_boxplot(data = long_algae %>% 
                 filter(FA %in% c("14:0", "16:0", "16:1w7c", "18:0", "18:1w7c", "18:3w3", "18:4w3c", "18:1w9c", "20:4w6", "20:5w3")), 
               aes(group = FA, y = proportion), color = "black", alpha = 0, show.legend = FALSE) +
  theme_classic() +
  scale_colour_viridis(discrete = TRUE, end = 0.9) +
  labs(x = "Fatty Acid", y = "Proportion of Total Fatty Acids") +
  scale_x_discrete(labels = c("14:0" = "MYR (14:0)",
                              "16:0" = "PAL (16:0)", 
                              "16:1w7c" = "PALO (16:1w7c)", 
                              "18:0" = "STE (18:0)",
                              "18:1w7c" = "VAC (18:1w7c)",
                              "18:3w3" = "LIN (18:3w3)",
                              "18:4w3c" = "SDA (18:4w3c)",
                              "18:1w9c" = "OLE (18:1w9c)", 
                              "20:4w6" = "ARA (20:4w6)", 
                              "20:5w3" = "EPA (20:5w3)")) +
  coord_flip()

# 18:1n7, 18:3n3 < 0.1, rest = 0.1 < other FA

# stacked barplot of percent contribution of each FA to total FA (dominant)
ggplot(filter(long_algae, FA %in% c("14:0", "16:0", "16:1w7c", "18:0", "18:1w7c", "18:3w3", "18:4w3c", "18:1w9c", "20:4w6", "20:5w3")) %>%
         group_by(genusSpecies, FA) %>%
         summarise(proportion = mean(proportion)) %>%
         ungroup(),
       aes(x = proportion, y = genusSpecies, fill = factor(FA, levels = c("14:0", "18:0", "16:1w7c", "18:1w7c", "18:3w3", "18:4w3c", "18:1w9c", "16:0", "20:5w3", "20:4w6")))) +
  geom_col(position = "stack", color = "black") +
  scale_fill_viridis(discrete = TRUE, option = "C") +
  theme_classic() +
  labs(fill = "Fatty Acid") +
  xlab("Mean Proportion of Total Fatty Acids") +
  ylab("Species")


## PAIRS PANEL INVERTS

grad_site <- long_inverts 

grad_site$SiteID <- grad_site$SiteID %>%
  recode('1' = '01',
         '2' = '02',
         '3' = '03',
         '4' = '04',
         '5' = '05',
         '6' = '06',
         '7' = '07',
         '8' = '08',
         '9' = '09',
         '10' = '10',
         '11' = '11',
         '12' = '12',
         '13' = '13',
         '14' = '14',
         '15' = '15')


# create site/lat, and add gradient - taken from NIC-Klein-Midpoint-Annual (email from Chuck)
SiteID <- c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15')
lat <- c(-65.94492, -67.5567, -68.1758, -68.6921, -67.5488, -66.87732, -66.0894, -65.5131, -65.1043, -64.9002, -64.7932, -64.7720, -66.0251, -64.7793, -65.2402)
grad <- c('56.44898', '71.5102', '68.61224', '87.67347', '57.87755', '82.89796', '58.36735', '73.95918', '53.5102', '36.12245', '41.10204', '37.26531', '76.77551', '41.06122', '62.85714')

site_lat <- data.frame(SiteID, lat, grad)

all_site <- left_join(grad_site, site_lat, by = "SiteID")
all_site$lat <- as.factor(all_site$lat)

# pivot data wide for pairs
all_site_wide <- all_site %>%
  select(FA, genusSpecies, proportion, FAsampleName, SiteID, lat, grad) %>%
  pivot_wider(names_from = FA, values_from = proportion, values_fill = 0)
# remove zero columns
all_site_pairs <- all_site_wide %>%
  select(where(~ any(. != 0)))
# create mean values
all_site_means <- all_site_pairs %>%
  select(-c('SiteID', 'FAsampleName')) %>%
  group_by(genusSpecies, lat) %>%
  summarise(across(everything(), mean))


# reduce number of FAs
all_site_reduced <- all_site_pairs %>% # switch between lat/grad for gradients here
  select(genusSpecies, grad,'16:1w7c', '18:1w7c', '18:1w9c', '18:2w6c', '20:0', '20:4w6', '22:0') %>%
  #filter(genusSpecies == 'Odontaster validus')
  #filter(genusSpecies == 'Nacella concinna')
  #filter(genusSpecies == 'Sterechinus neumayeri')
  #filter(genusSpecies == 'Neosmilaster georgianus')
  #filter(genusSpecies == 'Cnemidocarpa sp.')
  #filter(genusSpecies == 'Isotealia antarctica')
  #filter(genusSpecies == 'Perknaster fuscus')
  #filter(genusSpecies == 'Margarella antarctica')
  #filter(genusSpecies == 'Prostebbingia gracilis')
  filter(genusSpecies == 'Dendrilla membranosa')
#filter(genusSpecies == 'Gondogeneia antarctica')
pairs.panels(all_site_reduced[,2:9])


# FAs of gradient interest (>|50| Pearson correlation [moderate or stronger]):
# ('20:5w3', '18:4w3', '18:3w6', '18:3w3', '14:0', '16:0', '18:0',
# '16:1w7c', '18:1w7c', '18:1w9c', '18:2w6c', '20:0', '20:4w6', '22:0')
#
## S. neumayeri
#     20:4w6 = -0.50
#     22:0 = -0.53
## N. georgianus
#     14:0 = 0.50
## P. gracilis
#     16:0 = 0.70
#     22:0 = -0.77
## D. membranosa
#     18:4w3 = 0.77
#     18:3w3 = 0.56
#     14:0 = 0.65
#     18:0 = -0.84
#     16:1w7 = 0.72
#     18:1w9  = 0.55
#     18:2w6 = 0.53
#     20:0 = -0.59
#     20:4w6 = -0.54
#     22:0 = -0.55
## G. antarctica
#     16:1w7 = -0.53
#     18:3w6 = 0.68
#     14:0 = 0.72
## M. antarctica
#     20:4w6 = -0.56
# 


# FIXING TABLE ERROR FOR SUPP 3

# calc mean and sd of each FA for each sp
FA_means <- all_species %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  mutate(across(c(`8:0`:`24:1w9`), function(x) x*100)) 
FA_means <- FA_means %>%
  select(revisedSpecies, phylum, `8:0`:`24:1w9`) %>%
  group_by(phylum, revisedSpecies) %>%
  filter(!is.na(`8:0`)) %>%
  summarise(across(`8:0`:`24:1w9`, list(mean = mean, sd = sd))) 
FA_means <- FA_means %>%
  mutate(across(where(is.numeric), round, 3))
FA_means <- as.data.frame(t(FA_means)) 
colnames(FA_means) <- FA_means[2,]

FA_means <- FA_means %>%
  select(`Curdiea racovitzae`, `Desmarestia anceps`, `Desmarestia antarctica`,
         `Georgiella confluens`, `Meridionella antarctica`, `Pachymenia orbicularis`,
         `Picconiella plumosa`)
FA_means <- rownames_to_column(FA_means)
FA_redMean <- FA_means %>%
  filter(grepl('mean', rowname))
FA_redSD <- FA_means %>%
  filter(grepl('sd', rowname))
FA_redSD <- as.data.frame(lapply(FA_redSD, function(x) paste("±", x, sep=" ")))


write_csv(FA_redMean, "Means.csv")
write_csv(FA_redSD, "SD.csv")
