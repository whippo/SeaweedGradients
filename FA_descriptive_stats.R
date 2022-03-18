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
# 2022-02-08 Basically abandoning all previous work. Starting from scratch using
#             Browser instead of Insight. Started w/ batch 1
# 2022-03-07 Finished batch 1 and 2, used code to create basic summary stats
# 2022-03-18 Finished batch 3, added to analyses

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LOAD PACKAGES                                                                ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(tidyverse)
library(vegan)
library(viridis)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# READ IN AND PREPARE DATA                                                     ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# dummy_data <- read_csv("Data/Biomarkers/FattyAcids/dummy_sample_values.csv")

# OLD Browser data
# Gradients19_FA_Concs <- read_csv("Data/Biomarkers/FattyAcids/Gradients19_FA_Concs.csv", 
#                                                  col_types = cols(Conc = col_double(), 
#                                                  Date.anal = col_character(), Notes = col_character()))

# Insight Data
# Insight_concs <- read_csv("Data/Biomarkers/FattyAcids/Gradients_FA_Concs_Insight.csv", 
#                                col_types = cols(Area = col_number(), 
#                                                Conc = col_number(),
#                                                          Date_Insight_Anal = col_character(), 
#                                                          Number_ID = col_character()))
# Insight_concs <- Insight_concs %>%
#   mutate(across(where(is.character), ~na_if(., "----")))

# Need to be connected to Dropbox to access LINUX
Whippo_FA_extraction_log <- read_csv("~/Dropbox/OSF/Fatty Acid Extractions/Whippo_FA_extraction_log.csv") #linux
# WINDOWS
Whippo_FA_extraction_log <- read_csv("C:/Users/rossw/Dropbox/OSF/Fatty Acid Extractions/Whippo_FA_extraction_log.csv") #windows


# batch 1 import
# batch1 <- read_csv("Data/Biomarkers/FattyAcids/batch1_gradients19_rawquants.csv")
  
# batch1$Conc <- batch1$Conc %>%
#   replace("-----", 0)
#   mutate(Conc_num = case_when(Conc == "-----" ~ "0")) 
#   select(c('Data Filename', 'FA', 'Ret Time', 'Conc'))
  
all_batch <- read_csv("Data/Biomarkers/FattyAcids/core_spp_final_concs.csv")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# MANIPULATE DATA                                                              ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

############### Dummy Data

# convert data to long form
# dummy_trans <- t(dummy_data)
# dummy_trans <- as.data.frame(t(dummy_data[,-1]))
# colnames(dummy_trans) <- dummy_data$X1
# dummy_trans <- as.numeric(dummy_trans)
# dummy_nums <- mutate_all(dummy_trans, function(x) as.numeric(as.character(x)))
# dummy_nums <- dummy_nums %>% 
#  mutate(
#    across(everything(), ~replace_na(.x, 0))
#  )



#long_dummy <- dummy_trans %>%
#  pivot_longer(names_to = "")

########## Amount of Gradients samples done

Gradients_all <- Whippo_FA_extraction_log %>%
  filter(projectID == "Gradients2019") %>%
  filter(!str_detect(sampleID, "Blank"))

53/155

############### Actual data BROWSER NEW 2022!

# create species-genus ID column (BROWSER)
 grad_conc <- all_batch %>%
   mutate(species = case_when(startsWith(sample, "HIGR") ~ "H. grandifolius",
                              startsWith(sample, "DEME") ~ "D. menziesii",
                              startsWith(sample, "IRCO") ~ "I. cordata",
                              startsWith(sample, "PHAN") ~ "P. antarctica",
                              startsWith(sample, "PLCA") ~ "P. cartilagineum",
                              startsWith(sample, "MYMA") ~ "M. manginii"))



# DEME only for gradient analysis

 deme <- grad_conc %>%
   filter(species == "D. menziesii") %>%
   mutate(site = str_sub(sample, 6, 7))


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


### Figure  BROWSER (NEW 2022-03-07)

# MDS

# pivot data wide for mds
grad_conc_wide <- grad_conc %>%
  select(FA, species, proportion, sample) %>%
  pivot_wider(names_from = FA, values_from = proportion, values_fill = 0)

# select EPA, ARA, SDA, PAL, OLE, LIN, VAC, and dominant sats (16, 18)
sub_wide <- grad_conc_wide %>%
  select(species, sample, "14:0", "16:0", "16:1w7c", "18:0", "18:1w7c", "18:3w3", "18:4w3c", "18:1w9c", "20:4w6", "20:5w3")
sub_wide_trans <- (sub_wide[,3:12])
sub_wide_trans <- bind_cols(sub_wide[1:2], sub_wide_trans)

batch_1_2_MDS <- metaMDS(sub_wide_trans[3:12], autotransform = FALSE, distance = "chisq")
batch_1_2_MDS_points <- batch_1_2_MDS$points
batch_1_2_MDS_points <- data.frame(batch_1_2_MDS_points)
plot_data_batch_1_2 <- data.frame(batch_1_2_MDS_points, sub_wide[,1])

library(plyr)
# create the list of points that will connect the 'hulls' together from your nMDS point data
chulls_tax <- ddply(plot_data_batch_1_2, .(species), function(df) df[chull(df$MDS1, df$MDS2), ])
# DETACH PLYR so it won't mess with anything!
detach(package:plyr)

# create vectors to plot over MDS
scrs <- as.data.frame(scores(batch_1_2_MDS, display = "sites"))
scrs <- cbind(scrs, species = sub_wide_trans$species)

vf <- envfit(batch_1_2_MDS, sub_wide_trans[3:12], perm = 999)

spp.scrs <- as.data.frame(scores(vf, display = "vectors"))
spp.scrs <- cbind(spp.scrs, FA = rownames(spp.scrs))
spp.scrs <- spp.scrs %>%
  rename(MDS1 = NMDS1,
         MDS2 = NMDS2)

ggplot(plot_data_batch_1_2, aes(x=MDS1, y=MDS2)) +
  coord_fixed() +
  geom_point(size = 4, aes(color = species)) + # set size of points to whatever you want
  scale_color_viridis(discrete = TRUE, end = 0.9) + # my favorite color-blind and b&w friendly palette, look at the viridis package for more details
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2),
               arrow = arrow(length = unit(0.25, "cm")), color = "grey") +
  geom_text(data = spp.scrs, aes(label = FA), 
            size = 3) +
  geom_polygon(data = chulls_tax,
               aes(x = MDS1, y = MDS2, color = species), 
               fill = NA) + # optional 'hulls' around points
  theme_classic()  # optional, I just like this theme










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
ggplot(filter(grad_conc, FA %in% c("14:0", "16:0", "16:1w7c", "18:0", "18:1w7c", "18:3w3", "18:4w3c", "18:1w9c", "20:4w6", "20:5w3")), aes(x = FA, y = proportion, colour = species)) +
  geom_point(size = 4, position = position_jitter(width = .1)) +
  geom_boxplot(data = grad_conc %>% 
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
ggplot(filter(grad_conc, FA %in% c("14:0", "16:0", "16:1w7c", "18:0", "18:1w7c", "18:3w3", "18:4w3c", "18:1w9c", "20:4w6", "20:5w3")) %>%
         group_by(species, FA) %>%
         summarise(proportion = mean(proportion)) %>%
         ungroup(),
         aes(x = proportion, y = species, fill = factor(FA, levels = c("14:0", "18:0", "16:1w7c", "18:1w7c", "18:3w3", "18:4w3c", "18:1w9c", "16:0", "20:5w3", "20:4w6")))) +
  geom_col(position = "stack", color = "black") +
  scale_fill_viridis(discrete = TRUE, option = "C") +
  theme_classic() +
  labs(fill = "Fatty Acid") +
  xlab("Mean Proportion of Total Fatty Acids") +
  ylab("Species")









# stacked barplot of percent contribution of each FA to total FA (trace)
ggplot(filter(grad_conc, proportion > 0.01 & proportion < 0.1 & FA != '19:0') %>%
         mutate(species = fct_reorder(species, proportion, sum, .desc = TRUE)), aes(x = proportion, y = species, fill = FA)) +
  geom_col(position = "stack", color = "black") +
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