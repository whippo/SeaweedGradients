#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                                                                ##
# Antarctica Algae Descriptive Stats                                             ##
# Data are current as of 2022-02-08                                              ##
# Data source: Antarctic Gradients 2019                                          ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2022-03-07                                                        ##
#                                                                                ##
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# SUMMARY:

# Script to generate descriptive statistical summaries of core fatty acid data 
# from Antarctic Gradients 2019 project. 


# Required Files (check that script is loading latest version):
# Whippo_FA_extraction_log.csv
# gradients2019_corespecies_FA_QAQC.csv


# Associated Scripts:
# FA_data_pipeline.R

# TO DO 
# REMOVE c19 from everything!!!! (and in Quarto too!)

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

# 2021-04-26 Script created
# 2021-05-14 Updated metadata, import and treatment of Insight data
# 2022-02-08 Basically abandoning all previous work. Starting from scratch using
#             Browser instead of Insight. Started w/ batch 1
# 2022-03-07 Finished batch 1 and 2, used code to create basic summary stats
# 2022-03-18 Finished batch 3, added to analyses
# 2022-04-19 batches 4 and 5 now included
# 2022-04-20 cleaned up script to run from new QC data and removed unused parts

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LOAD PACKAGES                                                                ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(tidyverse) # data cleanup
library(vegan) # 'community' analyses
library(viridis) # color palette
library(psych) # pairs panel
library(ggfortify) # PCA visualizations

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# READ IN AND PREPARE DATA                                                     ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Read in batch log for list of all samples in core analysis
# Need to be connected to Dropbox to access 
# LINUX
Whippo_FA_extraction_log <- read_csv("~/Dropbox/OSF/Fatty Acid Extractions/Whippo_FA_extraction_log.csv") #linux
# WINDOWS
Whippo_FA_extraction_log <- read_csv("C:/Users/rossw/Dropbox/OSF/Fatty Acid Extractions/Whippo_FA_extraction_log.csv") #windows

# read in metadata to count all samples
collections <- read_csv("Data/Biomarkers/FattyAcids/B-236_Fatty-Acid_Collections_QAQC.csv")


# Read in all core species from data pipeline
all_species <- read_csv("Data/Biomarkers/FattyAcids/gradients2019_corespecies_FA_QAQC.csv")

# Pull out algae only
all_algae <- all_species %>%
  filter(type == 'algae')

# Create simplified long dataset for analysis of algae
long_algae <- all_algae %>%
  select(FAsampleName, batch, genusSpecies, `8:0`:`22:4w3`) %>%
  pivot_longer(cols = `8:0`:`22:4w3`, names_to = 'FA', values_to = 'proportion')

# subset inverts
all_inverts <- all_species %>%
  filter(type == "invert")

# Create simplified long dataset for analysis of inverts
long_inverts <- all_inverts %>%
  select(FAsampleName, SiteID, genusSpecies, `8:0`:`22:4w3`) %>%
  pivot_longer(cols = `8:0`:`22:4w3`, names_to = 'FA', values_to = 'proportion')
  

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# DATA SUMMARY                                                                 ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

###### Amount of Gradients samples done (from batch log)

# Gradients_all <- Whippo_FA_extraction_log %>%
#   filter(projectID == "Gradients2019") %>%
#   filter(!str_detect(sampleID, "Blank"))

# 106/155 = 0.683871 updated through batch 6



# total reps algae per site table:

# pull out letter and site ID to join with metadata
lettersite <- all_species %>%
  select(SiteID, LetterID)


sample_counts_algae <- collections %>%
  select(genusSpecies, class, LetterID) %>%
  group_by(class, genusSpecies, LetterID) %>%
  count(genusSpecies)

sample_algae_wide <- sample_counts_algae %>%
  pivot_wider(genusSpecies, names_from = 'LetterID', values_from = 'n')


###### DEME only for gradient analysis

 deme <- long_algae %>%
   filter(genusSpecies == "Desmarestia menziesii") %>%
   mutate(site = str_sub(FAsampleName, 6, 7))


 site <- c('02', '03', '04', '05', '07', '08', '09', '10', '11', '12', '13', '14', '15')
 lat <- c(-67.5567, -68.1758, -68.6921, -67.5488, -66.0894, -65.5131, -65.1043, -64.9002, -64.7932, -64.7720, -66.0251, -64.7793, -65.2402)

 site_lat <- data.frame(site, lat)

 deme_site <- left_join(deme, site_lat, by = "site")
 deme_site$lat <- as.factor(deme_site$lat)
 
 # add gradient - taken from visual of sampling sites figure
 grad <- c('80', '80', '90', '70', '60', '80', '60', '40', '40', '30', '80', '40', '70')
 site_grad <- data.frame(site_lat, grad)

 
 # major FA by latitude
 ggplot(filter(deme_site, proportion > 0.1 & FA !="19:0"), aes(x = lat, y = proportion, colour = FA)) +
   geom_point() +
   geom_point(size = 4) +
    theme_classic() +
   scale_colour_viridis(discrete = TRUE, end = 0.9) +
   labs(x = "Latitude", y = "Proportion of Total Fatty Acids") +
   geom_line(data = deme_site %>%
               filter(proportion > 0.1 & FA !="19:0") %>%
               group_by(lat, FA) %>%
               summarise("mean proportion" = mean(proportion)), 
             aes(x = lat, y = `mean proportion`, group = FA, colour = FA), lwd = 1) +
   coord_flip()
 
 deme_grad <- left_join(deme_site, site_grad, by = "site")
 
 # major FA by gradient
 ggplot(filter(deme_grad, proportion > 0.1 & FA !="19:0"), aes(x = grad, y = proportion, colour = FA)) +
   geom_point() +
   geom_point(size = 4) +
   theme_classic() +
   scale_colour_viridis(discrete = TRUE, end = 0.9) +
   labs(x = "Ice Cover (%)", y = "Proportion of Total Fatty Acids") +
   geom_line(data = deme_grad %>%
               filter(proportion > 0.1 & FA !="19:0") %>%
               group_by(grad, FA) %>%
               summarise("mean proportion" = mean(proportion)), 
             aes(x = grad, y = `mean proportion`, group = FA, colour = FA), lwd = 1) +
   coord_flip()
 
 # 20:4w6 and 20:5w3 by latitude
 ggplot(filter(deme_site, FA %in% c("20:4w6", "20:5w3")), aes(x = lat, y = proportion, colour = FA)) +
   geom_point(size = 4) +
   geom_smooth(aes(group = FA), method = lm, formula = y ~ x, fill = "grey") + 
   theme_classic() +
   scale_colour_viridis(discrete = TRUE, end = 0.9) +
   labs(x = "Latitude", y = "Proportion of Total Fatty Acids") 

 # 20:4w6 and 20:5w3 by gradient
 ggplot(filter(deme_grad, FA %in% c("20:4w6", "20:5w3")), aes(x = grad, y = proportion, colour = FA)) +
   geom_point(size = 4) +
   geom_smooth(aes(group = FA), method = lm, formula = y ~ x, fill = "grey") + 
   theme_classic() +
   scale_colour_viridis(discrete = TRUE, end = 0.9) +
   labs(x = "Ice Cover (%)", y = "Proportion of Total Fatty Acids") 


 
 
###### All algal species analyses
 
 # pivot data wide for multivariate
 grad_conc_wide <- long_algae %>%
   select(FA, genusSpecies, proportion, FAsampleName) %>%
   pivot_wider(names_from = FA, values_from = proportion, values_fill = 0)
 # remove zero columns
 grad_conc_PCA <- grad_conc_wide %>%
   select(where(~ any(. != 0)))
 
 
 ### PERMANOVA 
 
 # algal FA for adonis
 FA_only <- grad_conc_PCA %>%
   select(`8:0`:`18:3w1`)

adonis(abs(FA_only) ~ genusSpecies, data = grad_conc_PCA, method = 'bray')

# PCA
 


 # run PCA
PCA_results <-  prcomp(grad_conc_PCA[,c(3:66)], scale = TRUE)
PCA_results$rotation <- -1*PCA_results$rotation
PCA_results$rotation
#reverse the signs of the scores
PCA_results$x <- -1*PCA_results$x
#display the first six scores
head(PCA_results$x)

# plot how much variance
plot(PCA_results)

# calc variance explained
var_explained <- PCA_results$sdev^2/sum(PCA_results$sdev^2)

# biplot of 2 most important axes
PCA_results$x %>% 
  as.data.frame %>%
  ggplot(aes(x=PC1,y=PC2)) + geom_point(aes(color=grad_conc_PCA$genusSpecies),size=4) +
  theme_minimal() +
  scale_color_viridis(discrete = TRUE) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))

autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3)



# run PCA for DEME only
deme_pca <- grad_conc_PCA %>%
  filter(genusSpecies == "Desmarestia menziesii") %>%
  select(where(~ any(. != 0)))
PCA_results <-  prcomp(deme_pca[,c(3:58)], scale = TRUE)
PCA_results$rotation <- -1*PCA_results$rotation
PCA_results$rotation
#reverse the signs of the scores
PCA_results$x <- -1*PCA_results$x
#display the first six scores
head(PCA_results$x)

# plot how much variance
plot(PCA_results)

# calc variance explained
var_explained <- PCA_results$sdev^2/sum(PCA_results$sdev^2)

# biplot of 2 most important axes
PCA_results$x %>% 
  as.data.frame %>%
  ggplot(aes(x=PC1,y=PC2)) + 
  geom_point() +
  theme_minimal() +
  scale_color_viridis(discrete = TRUE) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))

 


# MDS

# pivot data wide for mds
grad_conc_wide <- long_algae %>%
  select(FA, genusSpecies, proportion, sample) %>%
  pivot_wider(names_from = FA, values_from = proportion, values_fill = 0)

# select EPA, ARA, SDA, PAL, OLE, LIN, VAC, and dominant sats (16, 18) 
sub_wide <- grad_conc_wide %>% # fix 16:0 error in inverts
  select(genusSpecies, sample, "14:0", "16:0", "16:1w7c", "18:0", "18:1w7c", "18:3w3", "18:4w3c", "18:1w9c", "20:4w6", "20:5w3")
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

ggplot(plot_data_batch_1_2, aes(x=MDS1, y=MDS2)) +
  coord_fixed() +
  geom_point(size = 4, aes(color = genusSpecies)) + # set size of points to whatever you want
  scale_color_viridis(discrete = TRUE, end = 0.9) + # my favorite color-blind and b&w friendly palette, look at the viridis package for more details
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2),
               arrow = arrow(length = unit(0.25, "cm")), color = "grey") +
  geom_text(data = spp.scrs, aes(label = FA), 
            size = 3) +
  geom_polygon(data = chulls_tax,
               aes(x = MDS1, y = MDS2, color = genusSpecies), 
               fill = NA) + # optional 'hulls' around points
  theme_classic()  # optional, I just like this theme



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


