#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                                                                ##
# Antarctica Algae FA Bsides SI Analyses                                         ##
# Script Created 2023-03-21                                                      ##
# Data source: Antarctic Gradients 2019                                          ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2023-03-21                                                        ##
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
# 1. Combine ggplot PCA plot with vectors  
# 2. Run PERMANOVA and PCA for SI and FA separately
#     - samples in common
#     - FA with additional samples
#     - 'All' samples
# 3. Fill in all values in paper table for FA and SI
# 4. Create PCA labeling phylum instead of species (or shapes)

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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LOAD PACKAGES                                                                ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(tidyverse) # data cleanup
library(vegan) # 'community' analyses
library(viridis) # color palette
library(psych) # pairs panel
library(ggfortify) # PCA visualizations
library(stringi) # order FA's in columns

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# READ IN AND PREPARE DATA                                                     ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Read in all core species from data pipeline and remove duplicated 'all' FA
FASI_QAQC <- read_csv("Data/Biomarkers/FattyAcids/gradients2019_bsides_FASI_QAQC.csv")
all_species <- FASI_QAQC %>%
  filter(targetFA == "standards") %>%
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

  



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# DATA SUMMARY                                                                 ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




###### OVERLAPPING SAMPLES



### PERMANOVA 

# algal FA for adonis
marker_only <- overlap_species %>%
  select(`CN ratio`:`24:1w9`) 

adonis2(abs(marker_only) ~ revisedSpecies, data = overlap_species, method = 'bray', na.rm = TRUE)

# run PCA
PCA_results <-  rda(overlap_species[,c(21:63)], scale = TRUE)
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



###### FA VALUES ONLY



### PERMANOVA 

# algal FA for adonis
FA_only <- all_species %>%
  select(`8:0`:`24:1w9`) 

adonis2(abs(FA_only) ~ revisedSpecies, data = all_species, method = 'bray', na.rm = TRUE)

# PCA


# run PCA
PCA_results <-  rda(all_species[,c(24:63)], scale = TRUE)
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


# Proportional stacked barplot of FA composition for each species

# calc mean of each FA for each sp
FA_means <- all_species %>%
  select(revisedSpecies, phylum, `8:0`:`24:1w9`) %>%
  group_by(phylum, revisedSpecies) %>%
  summarise(across(`8:0`:`24:1w9`, mean)) %>%
  pivot_longer(cols = `8:0`:`24:1w9`, names_to = 'Fatty Acid', values_to = 'value') %>%
  filter(value >= 0.02295) %>%
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
  theme(axis.text.x = element_text(angle = 270, hjust = 0))
  

consumed %>%
  filter(item %in% c('green', 'purple', 'red', 'mussel', 'cucumber')) %>%
  ggplot() +
  geom_col(aes(x = reorder(pycnoID, eaten, function(x){ sum(x)}), y = eaten, fill = item), position = "fill") +
  scale_fill_viridis(discrete = TRUE, option = 5) +
  theme_bw() + 
  scale_x_discrete(label = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K')) +
  geom_text(data = totals, aes(x = pycnoID, y= 1.05, label = total, fill = NULL)) +
  labs(x = "", y = "Proportion Consumed") +
  guides(fill=guide_legend(title="Prey"))




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



