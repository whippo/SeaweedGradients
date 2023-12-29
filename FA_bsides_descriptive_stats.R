#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                                                                ##
# Antarctica Algae FA Bsides SI Analyses                                         ##
# Script Created 2023-03-21                                                      ##
# Data source: Antarctic Gradients 2019                                          ##
# R code prepared by Ross Whippo                                                 ##
# Last updated 2023-12-29                                                        ##
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
library(ggrepel)
library(ggpp)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LOAD FUNCTIONS                                                               ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# fuction for "%notin%
`%notin%` <- Negate(`%in%`)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# READ IN AND PREPARE DATA                                                     ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Read in all core species from data pipeline and remove duplicated 'all' FA
FASI_QAQC <- read_csv("Data/Biomarkers/FattyAcids/gradients2019_bsides_FASI_QAQC_new.csv")
all_species <- FASI_QAQC %>%
  select(!`19:0`) #%>%
  #select_if(~ !is.numeric(.) || sum(., na.rm = TRUE) != 0)

# Create simplified long dataset for analysis, remove non-overlapping samples
long_species <- all_species %>%
  select(ProjID, siteName, revisedSpecies, `Ice cover (NIC-Midpoint-Annual)`,
         targetFA, `CN ratio`:`24:1w9`) %>%
  filter(!is.na(`CN ratio`)) %>%
  filter(!is.na(`8:0`)) %>%
  pivot_longer(cols = `CN ratio`:`24:1w9`, names_to = 'marker', values_to = 'value')

# create wide dataset, remove non-overlapping samples
overlap_species <- all_species %>%
  filter(!is.na(`CN ratio`)) %>%
  filter(!is.na(`8:0`)) #%>%
  #select_if(~ !is.numeric(.) || sum(., na.rm = TRUE) != 0)

# SI only wide dataset
SI_wide <- all_species %>%
  filter(!is.na(`CN ratio`))

# FA only wide dataset
FA_wide <- all_species %>%
  filter(!is.na(`8:0`))





#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FIGURE 2                                                                     ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# Proportional stacked barplot of FA composition for each species

# calc mean of each FA for each sp
FA_quartile <- long_species %>%
  filter(marker %notin% c("d13C", "d15N", "CN ratio")) %>%
  filter(value != 0)
summary(FA_quartile$value)

FA_means <- FA_wide %>%
  select(revisedSpecies, phylum, `8:0`:`24:1w9`) %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  group_by(phylum, revisedSpecies) %>%
  summarise(across(`8:0`:`24:1w9`, mean)) %>%
  pivot_longer(cols = `8:0`:`24:1w9`, names_to = 'FAs', values_to = 'value') %>%
  filter(value >= 0.037) %>% # mean all FA values in dataset
  mutate(phylum = case_when(phylum == "Chlorophyta" ~ "",
                            phylum == "Ochrophyta" ~ "Ochrophyta",
                            phylum == "Rhodophyta" ~ "Rhodophyta")) 
  FA_means <- FA_means %>%
    mutate(`Fatty Acid` = case_when(
            `FAs` %in% "16:1w7c" ~ "16:1ω7",
            `FAs` %in% "16:3w3" ~ "16:3ω3",
            `FAs` %in% "18:1w7c" ~ "18:1ω7",
            `FAs` %in% "18:1w9c" ~ "18:1ω9",
            `FAs` %in% "18:2w6c" ~ "18:2ω6",
            `FAs` %in% "18:3w3" ~ "18:3ω3",
            `FAs` %in% "18:4w3c" ~ "18:4ω3",
            `FAs` %in% "20:3w6" ~ "20:3ω6",
            `FAs` %in% "20:4w6" ~ "20:2ω6",
            `FAs` %in% "20:5w3" ~ "20:5ω3",
            `FAs` %in% "22:5w3" ~ "22:5ω3",
            TRUE ~ as.character(FAs)))

FA_means %>%
  ggplot() +
  geom_col(aes(x = revisedSpecies, y = value, fill = `Fatty Acid`), position = "fill") +
  scale_fill_viridis(discrete = TRUE, option = 6) +
  theme_bw() +
  labs(x = "Species", y = "Mean Proportional Composition") +
  guides(fill = guide_legend(title = "Fatty Acid")) +
  facet_grid(cols = vars(phylum), scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.25, face = "italic")) +
  theme(axis.title.x = element_text(vjust = -1)) +
  theme(axis.title.y = element_text(vjust = 2.5))

# size = 10x6 (use cairo pdf)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FIGURE 3                                                                     ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# cluster analysis of FA only

no_diatoms_meanFA <- FA_wide %>%
  select(revisedSpecies, phylum, `8:0`:`24:1w9`) %>%
  group_by(phylum, revisedSpecies) %>%
  summarise(across(`8:0`:`24:1w9`, mean)) 


Alg_dist <- vegdist(no_diatoms_meanFA[,3:46])
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

# 8 X 8 cairo pdf

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FIGURE 4                                                                     ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


##################
###### FA VALUES ONLY, ORDER LEVEL ANALYSIS
# algal FA for adonis
FA_only <- FA_wide %>%
  select(`8:0`:`24:1w9`) 
FA_tax <- FA_wide %>%
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

# PERMANOVA of order and family ALL FA
adonis2(abs(FA_tax[,24:67]) ~ order, data = FA_tax, method = 'bray', na.rm = TRUE)
adonis2(abs(FA_tax[,24:67]) ~ family, data = FA_tax, method = 'bray', na.rm = TRUE)

# PERMANOVA of order and family REDUCED FA
FA_tax_reduced <- FA_tax %>%
  select(`20:5w3`, `20:4w6`, `16:0`, `18:3w3`, `18:4w3c`, `18:1w9c`, `18:1w7c`,
                         `revisedSpecies`, `order`, `phylum`, `family`)
adonis2(abs(FA_tax_reduced[,1:7]) ~ order, data = FA_tax, method = 'bray', na.rm = TRUE)
adonis2(abs(FA_tax_reduced[,1:7]) ~ family, data = FA_tax, method = 'bray', na.rm = TRUE)

###############
FA_matrix <- FA_only

# run the nMDS
FA_mds <- metaMDS(FA_matrix)
# extract the 'points' from the nMDS that you will plot in ggplot2
FA_mds_points <- FA_mds$points
# turn those plot points into a dataframe that ggplot2 can read
FA_mds_points <- data.frame(FA_mds_points)
# join your plot points with your summed species observations from each habitat type
plot_data_tax <- data.frame(FA_mds_points, FA_tax[,c(7,8,68,69)])
plot_data_tax <- plot_data_tax %>%
  rename("division" = "phylum")


# run the ggplot
phylum_plot <- 
  ggplot(plot_data_tax, aes(x=MDS1, y=MDS2, 
                            fill = division)) +  
  labs(x = "nMDS1", y = "nMDS2") +
  theme_classic() + 
  geom_point(pch = 21, size = 2, color = "black") + 
  scale_fill_viridis(discrete = TRUE, begin = 0.2, end = 0.9, option = "G", name = "division") 
phylum_leg <- as_ggplot(get_legend(phylum_plot))

order_plot <- 
  ggplot(plot_data_tax, aes(x=MDS1, y=MDS2, # note pch and color are linked to *factors*
                            fill = order)) +  
  labs(x = "nMDS1", y = "nMDS2") +
  theme_classic() +
  geom_point(pch = 21, size = 2, color = "black") +  # set size of points to whatever you want
  guides(color=guide_legend(ncol=2)) +
  # geom_polygon(data=chulls_tax, aes(x=MDS1, y=MDS2, group=Habitat), fill=NA) + # optional 'hulls' around points
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "B", name = "order") # my favorite color-blind and b&w friendly palette, look at the viridis package for more details 
order_leg <- as_ggplot(get_legend(order_plot))


family_plot <- 
  ggplot(plot_data_tax, aes(x=MDS1, y=MDS2, # note pch and color are linked to *factors*
                            fill = family)) +  
  labs(x = "nMDS1", y = "nMDS2") +
  theme_classic() + 
  geom_point(pch = 21, size = 2, color = "black") +  # set size of points to whatever you want
  # geom_polygon(data=chulls_tax, aes(x=MDS1, y=MDS2, group=Habitat), fill=NA) + # optional 'hulls' around points
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "H", name = "family") # my favorite color-blind and b&w friendly palette, look at the viridis package for more details 
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

annotate_figure(FigureMDS, top = text_grob("2D stress = 0.11", size = 10))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FIGURE 5                                                                     ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

###### SIMPER TO PULL OUT GEOM_TEXT
simper_FA <- FA_wide %>%
  select(`8:0`:`24:1w9`) 

full_algal_simper <- simper(simper_FA)
simpersum <- summary(full_algal_simper)
simpersum <- data.frame(unclass(simpersum),  # Convert summary to data frame
                        check.names = FALSE)
simpersum <- rownames_to_column(simpersum, "VALUE")
# pull out top 82% contributors to differences
topsimp <- simpersum %>%
  filter(total.cumsum < 0.83)
topsimp$VALUE

### PERMANOVA 

# algal FA for adonis
FA_only <- FA_wide %>%
  select(`8:0`:`24:1w9`) 

# species only
adonis2(abs(FA_only) ~ revisedSpecies, data = FA_wide, method = 'bray', na.rm = TRUE)
# division only
adonis2(abs(FA_only) ~ phylum, data = FA_wide, method = 'bray', na.rm = TRUE)

# PCA



# run PCA
PCA_results <-  rda(FA_only, scale = TRUE)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(FA_wide), 
                       rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
vscores <- rownames_to_column(vscores)
vscores <- vscores 

# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

xvalues <- c(-0.31, 0.18, -0.05, 0.30, 0.26, -0.09, -0.15)
yvalues <- c(-0.14, 0.29, -0.27, -0.13, 0.25, 0.23, -0.15)
# make final ggplot figure
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE, begin = 0.2, end = 0.9, option = "G", guide = guide_legend(title = "division")) +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  geom_point(aes(x = PC1, y = PC2, fill = phylum), pch = 21, color = "black", size = 4) +
 # geom_text(data = subset(vscores, rowname %in% c("16:0", # -0.31, -0.13 
#                                                  "20:4w6", # -0.09, 0.23
#                                                  "18:4w3c", # 0.26, 0.25
#                                                  "18:1w9c", # 0.18, 0.29  
#                                                  "18:3w3", # -0.30, -0.13
#                                                  "18:1w7c", # -0.05, -0.27
#                                                  "20:5w3" # -0.15, -0.15
#                                                  )), aes(x = xvalues, y = yvalues, label = rowname), col = 'red') +
  geom_text_repel(data = subset(vscores, rowname %in% c(topsimp$VALUE)),
                  aes(x = PC1, y = PC2, label = rowname), color = "red",
                  position = position_nudge_center(x = 0.02, y = 0.02, 0, 0),
                  min.segment.length = 1) +
   theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))


# size = 8x6 cairo pdf


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FIGURE 6                                                                     ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


###### FA VALUES ONLY NO DIATOMS REDUCED FA



### PERMANOVA 

# algal FA for adonis
FA_only <- FA_wide %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`20:5w3`, `20:4w6`, `16:0`, `18:3w3`, `18:4w3c`, `18:1w9c`, `18:1w7c`) 

# species only
adonis2(abs(FA_only) ~ revisedSpecies, data = filter(FA_wide, revisedSpecies != "Benthic diatoms"), method = 'bray', na.rm = TRUE)
# division only
adonis2(abs(FA_only) ~ phylum, data = filter(FA_wide, revisedSpecies != "Benthic diatoms"), method = 'bray', na.rm = TRUE)

# PCA


# run PCA
PCA_results <-  rda(FA_only, scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
#ev <- PCA_results$CA$eig
#ev>mean(ev)
# proportion explained
#barplot(ev, main="eigenvalues", col="bisque", las=2)
#abline(h=mean(ev), col="red")
#legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
#biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
#biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
#autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(filter(FA_wide, revisedSpecies != "Benthic diatoms")), rownames_to_column(data.frame(uscores)), by = "rowname")
# add phylum shapes
uscores2 <- uscores1 %>%
  mutate(phylumShape = case_when(phylum == "Ochrophyta" ~ "24",
                   phylum == "Rhodophyta" ~ "22",
                   phylum == "Chlorophyta" ~ "21"))

vscores <- data.frame(PCA_results$CA$v)
vscores <- rownames_to_column(vscores)
vscores <- vscores 
# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

# make final ggplot figure
#xvalues <- c(0.4, -0.10, 0.47, -0.35, -0.50, -0.40, 0.25)
#yvalues <- c(0.41, -0.75, -0.27, 0.29, 0.06, 0.14, 0.45)
ggplot(uscores2) + 
  scale_fill_viridis(discrete = TRUE, guide = guide_legend(title = "species", label.theme = element_text(face = "italic"))) +
  #scale_color_viridis(discrete = TRUE, guide = guide_legend(title = "species", label.theme = element_text(face = "italic")))  +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_color_manual(values = rep("black", 31)) +
  #geom_text(data = vscores, aes(x = xvalues, y = yvalues, label = rownames(vscores)), col = 'red') +
  geom_point(aes(x = PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies, shape = phylum),
                 size = 4) +
  geom_text_repel(data = subset(vscores, rowname %in% c(topsimp$VALUE)),
                  aes(x = PC1, y = PC2, label = rowname), color = "red",
                  position = position_nudge_center(x = 0.02, y = 0.02, 0, 0),
                  min.segment.length = 1) +
  labs(shape = "division", fill = "species") +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  guides(fill = guide_legend(override.aes = list(shape=21), label.theme = element_text(face = "italic", size = 9)), color = "none") +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))

# size = 11x6 cairo pdf


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FIGURE 7 - SI VALUES PANEL                                                   ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# algal SI for panels
SI_values <- SI_wide %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(c("revisedSpecies", "phylum", `CN ratio`:d13C)) %>%
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

# PERMANOVA on order and family
adonis2(abs(SI_values[,3:5]) ~ order, data = SI_values, method = 'bray', na.rm = TRUE)
adonis2(abs(SI_values[,3:5]) ~ family, data = SI_values, method = 'bray', na.rm = TRUE)


# Summarise By Group SPECIES
sumspecies <- SI_values %>% 
  group_by(revisedSpecies, phylum, family, order) %>% 
  summarise(count = n(),
            mC = mean(d13C), 
            sdC = sd(d13C), 
            mN = mean(d15N), 
            sdN = sd(d15N))

# deltas species
deltaspecies <- ggplot(sumspecies) +
  geom_errorbar(data = sumspecies, na.rm = FALSE,
                mapping = aes(x = mC, y = mN,
                              ymin = mN - sdN/sqrt(count), 
                              ymax = mN + sdN/sqrt(count)), 
                width = 0) +
  geom_errorbarh(data = sumspecies, na.rm = FALSE,
                 mapping = aes(x = mC, y = mN,
                               xmin = mC - sdC/sqrt(count),
                               xmax = mC + sdC/sqrt(count)),
                 height = 0) + 
  #geom_point(aes(x = mC, y = mN, color = revisedSpecies),
  #           size = 4) +
  geom_point(aes(x = mC, y = mN, fill = revisedSpecies), pch = 21, color = "black", size = 4) +
  scale_fill_viridis(discrete = TRUE, 
                      guide = guide_legend(title = "species",
                                           label.theme = element_text(face = "italic", size = 9)))  +
  labs(x = "\U03B4\U00B9\u00B3C" , y = "\U03B4\U00B9\U2075N") +
  theme_bw()
species_leg <- ggplot(sumspecies, aes(x = mC, y = mN, color = revisedSpecies)) +
  geom_point(size = 4) +
  theme_bw() +
  scale_color_viridis(discrete = TRUE, 
                      guide = guide_legend(title = "species",
                                           label.theme = element_text(face = "italic", size = 9), ncol = 4))
species_leg <- as_ggplot(get_legend(species_leg)) 

# deltas divisions
deltadivision <- ggplot(sumspecies) +
  scale_fill_viridis(discrete = TRUE, begin = 0.2, end = 0.9, option = "G", name = "division") +
  geom_errorbar(data = sumspecies, 
                mapping = aes(x = mC, y = mN,
                              ymin = mN - sdN/sqrt(count), 
                              ymax = mN + sdN/sqrt(count)), 
                width = 0) +
  geom_errorbarh(data = sumspecies, 
                 mapping = aes(x = mC, y = mN,
                               xmin = mC - sdC/sqrt(count),
                               xmax = mC + sdC/sqrt(count)),
                 height = 0, na.rm = TRUE) + 
 # geom_point(aes(x = mC, y = mN, color = phylum, group = phylum),
#             size = 4) +
  geom_point(aes(x = mC, y = mN, fill = phylum), pch = 21, color = "black", size = 4) +
  labs(x = "\U03B4\U00B9\u00B3C" , y = "\U03B4\U00B9\U2075N") +
  theme_bw()
division_leg <- as_ggplot(get_legend(deltadivision))

# CN ratios

ratios <- SI_values %>%
  mutate(phylum = case_when(phylum == "Chlorophyta" ~ "",
                            phylum == "Ochrophyta" ~ "Ochrophyta",
                            phylum == "Rhodophyta" ~ "Rhodophyta")) %>%
  ggplot(aes(x = revisedSpecies, y = `CN ratio`, fill = phylum), color = "black") + 
  scale_fill_viridis(discrete = TRUE, begin = 0.2, end = 0.9, option = "G") + 
  geom_errorbar(stat="summary", fun.data="mean_se", size = 1) +
  stat_summary(
    geom = "point",
    fun = "mean",
    size = 6,
    shape = 21
  ) +
  facet_grid(cols = vars(phylum), scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "species", y = "C:N ratio") +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.25, face = "italic")) +
  theme(axis.title.x = element_text(vjust = -1)) +
  theme(axis.title.y = element_text(vjust = 2.5))




FigureSI <- ggarrange(
  
  ggarrange(
    
    ggarrange(deltadivision, deltaspecies,
              labels = c("A", "B"),
              ncol = 2, nrow = 1,
              legend = "none"),
    
    ggarrange(
      division_leg, species_leg,
      ncol = 2, nrow = 1, align = "hv",
      widths = c(0.3, 1),
      legend = "none"),
    
    ncol = 1, nrow = 2, heights = c(1, 0.5), legend = "none"),
  
  ratios, labels = "C",
  label.y = 0.1,
  ncol = 1, nrow = 2, legend = "none"
)

FigureSI
# 10x16 portrait PDF








#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Online Resource Fig 1                                                        ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# cluster analysis of SI only

no_diatoms_meanSI <- SI_wide %>%
  select(revisedSpecies, phylum, `CN ratio`:`d13C`) %>%
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

# size = 7x7 cairo pdf

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FIGURE 8                                                                     ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

###### SI VALUES ONLY no diatoms

### PERMANOVA 

# algal SI for adonis
SI_only <- SI_wide %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`CN ratio`:d13C) 

# species
adonis2(abs(SI_only) ~ revisedSpecies, data = filter(SI_wide, revisedSpecies != "Benthic diatoms"), method = 'bray', na.rm = TRUE)
# division
adonis2(abs(SI_only) ~ phylum, data = filter(SI_wide, revisedSpecies != "Benthic diatoms"), method = 'bray', na.rm = TRUE)



# run PCA
PCA_results <-  rda(SI_wide[,c(21:23)], scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
#ev <- PCA_results$CA$eig
#ev>mean(ev)
# proportion explained
#barplot(ev, main="eigenvalues", col="bisque", las=2)
#abline(h=mean(ev), col="red")
#legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
#biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
#biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
#autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(filter(SI_wide, revisedSpecies != "Benthic diatoms")), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]
vscores <- rownames_to_column(vscores)
vscores <- vscores 

xvalues <- c(0.63, 0.35, 0.77)
yvalues <- c(-0.55, 0.90, 0.03)
# make final ggplot figure (Points scaled by 1.5)
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE, guide = guide_legend(title = "species", label.theme = element_text(face = "italic"))) +
  #scale_color_viridis(discrete = TRUE, guide = guide_legend(title = "Species", label.theme = element_text(face = "italic")))  +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_color_manual(values = rep("black", 31)) +
  #geom_text(data = vscores, aes(x = xvalues, y = yvalues, label = rownames(vscores)), col = 'red') +
  #geom_point(aes(x = PC1*1.5, y = PC2*1.5, fill = revisedSpecies, color = revisedSpecies,
   #              shape = phylum), size = 4) +
  geom_point(aes(x = PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies, shape = phylum),
             size = 4) +
  geom_text_repel(data = vscores,
                  aes(x = PC1, y = PC2, label = rowname), color = "red",
                  position = position_nudge_center(x = 0.02, y = 0.02, 0, 0),
                  min.segment.length = 1) +
  labs(shape = "division", fill = "species") +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  guides(fill = guide_legend(override.aes = list(shape=21), label.theme = element_text(face = "italic", size = 9)), color = "none") +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))

# size = 11x6 cairo pdf

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FIGURE 9                                                                     ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


##### WITHOUT DIATOMS

### PERMANOVA 

# algal FA for adonis
marker_only <- overlap_species %>%
  select(`CN ratio`:`24:1w9`) 
#######################


# algal FA/SI for PERMANOVA
overlap_perm <- overlap_species %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`CN ratio`:`24:1w9`, revisedSpecies)  %>%
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

# family
adonis2(abs(overlap_perm[,1:47]) ~ family, data = overlap_perm, method = "bray", na.rm = TRUE)
# order
adonis2(abs(overlap_perm[,1:47]) ~ order, data = overlap_perm, method = "bray", na.rm = TRUE)

# species
adonis2(abs(marker_only) ~ revisedSpecies, 
        data = overlap_species, method = 'bray', na.rm = TRUE)
# division
adonis2(abs(marker_only) ~ phylum, 
        data = overlap_species, method = 'bray', na.rm = TRUE)

####################

  # run PCA
PCA_results <-  rda(marker_only, scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
#ev <- PCA_results$CA$eig
#ev>mean(ev)
# proportion explained
#barplot(ev, main="eigenvalues", col="bisque", las=2)
#abline(h=mean(ev), col="red")
#legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
#biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
#biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
#autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(overlap_species), 
                       rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
vscores <- rownames_to_column(vscores)
vscores <- vscores 

# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]


# pull out top 98% contributors to differences
topsimp <- simpersum %>%
  filter(total.cumsum < 0.99) %>%
  add_row(VALUE = "CN ratio") %>%
  add_row(VALUE = "d15N") %>%
  add_row(VALUE = "d13C")

topsimp$VALUE


# make final ggplot figure
#set.seed(1)
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE, 
                     guide = guide_legend(title = "Species", label.theme = element_text(face = "italic"))) +
  scale_color_viridis(discrete = TRUE, 
                      guide = guide_legend(title = "Species", label.theme = element_text(face = "italic")))  +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_color_manual(values = rep("black", 31)) +
  #geom_point(aes(x = PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies,
        #         shape = phylum), size = 4) +
  geom_point(aes(x = PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies, shape = phylum),
             size = 4) +
  #geom_text(data = vscores, aes(x = PC1, y = PC2, label = rownames(vscores)),
  #          col = 'red',
  #          position = position_jitter(width = 0.05, height = 0.01, seed = 5)) +
  geom_text_repel(data = subset(vscores, rowname %in% c(topsimp$VALUE)),
            aes(x = PC1, y = PC2, label = rowname), color = "red",
            position = position_nudge_center(x = 0.02, y = 0.02, 0, 0),
            min.segment.length = 2) +
  labs(shape = "division", fill = "species") +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  guides(fill = guide_legend(override.aes = list(shape=21), label.theme = element_text(face = "italic", size = 9)), color = "none") +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))

# size 11x6 cairo pdf


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# TABLE 2                                                                      ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


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
FA_means_separate <- FA_means %>%
  select(revisedSpecies, phylum, published, `8:0`:`24:1w9`) %>%
  group_by(phylum, published) %>%
  filter(!is.na(`8:0`)) %>%
  summarise(across(`8:0`:`24:1w9`, list(mean = mean))) 
top_FA_separate <- FA_means_separate %>%
  pivot_longer(`8:0_mean`:`24:1w9_mean`, names_to = "FA", values_to = "Percent") %>%
  arrange(desc(Percent)) %>% 
  group_by(phylum, published) %>%
  slice(1:5) %>%
  mutate(rank = c(1:5)) %>%
  pivot_wider(names_from = "rank", values_from = c("FA", "Percent"))

# grand means

FA_means_together <- FA_means %>%
  select(revisedSpecies, phylum, `8:0`:`24:1w9`) %>%
  group_by(phylum) %>%
  filter(!is.na(`8:0`)) %>%
  summarise(across(`8:0`:`24:1w9`, list(mean = mean))) 
top_FA_together <- FA_means_together %>%
  pivot_longer(`8:0_mean`:`24:1w9_mean`, names_to = "FA", values_to = "Percent") %>%
  arrange(desc(Percent)) %>% 
  group_by(phylum) %>%
  slice(1:10) %>%
  mutate(rank = c(1:10)) %>%
  pivot_wider(names_from = "rank", values_from = c("FA", "Percent"))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# TABLE 4                                                                      ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


simper_FA <- FA_wide %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`8:0`:`24:1w9`) 

full_algal_simper <- simper(simper_FA)
simpersum <- summary(full_algal_simper)
simpersum <- data.frame(unclass(simpersum),  # Convert summary to data frame
                        check.names = FALSE)
simpersum <- rownames_to_column(simpersum, "VALUE")

# write_csv(simpersum, "simper.csv")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# SUPP TABLE 3                                                                 ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

###### BIOMARKER VALUES FOR TABLE SUPP 3 

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

means_only <- marker_means %>%
  filter(grepl('mean', rowname))

 write_csv(means_only, "means_only.csv")

sd_only <- marker_means %>%
  filter(grepl('sd', rowname))

 write_csv(sd_only, "sd_only.csv")

 write_csv(marker_means, "marker_means.csv")



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# VARIOUS SUMMARY STATS                                                        ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# From TABLE 2


# Phylum level means and range of common FA
faprint <- FA_means %>%
  group_by(phylum) %>%
  select(`16:0`, `20:5w3`, `20:4w6`, `16:1w7c`, `18:3w3`, `18:2w6c`) %>%
  summarise(across(
  everything(), 
    .fns = list(Mean = mean, min = min, max = max), na.rm = TRUE, 
    .names = "{col}_{fn}"
  ))
  

# From figure 6

# phylum level values for SI
siprint <- SI_values %>%
  group_by(phylum) %>%
  select(`CN ratio`:`d13C`) %>%
  summarise(across(
    everything(),
    .fns = list(Mean = mean, min = min, max = max, sd = sd), na.rm = TRUE,
    .names = "{col}_{fn}"
  ))





#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Online Resource Figure 2                                                     ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


######## nMDS
# algal SI for adonis
SI_only <- SI_wide %>%
  select(`CN ratio`:`d13C`) 
##############
########### SI VALUES ONLY NO DIATOMS, ORDER LEVEL ANALYSIS
# algal SI for adonis
SI_only <- SI_wide %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`CN ratio`:`d13C`) 
SI_tax <- SI_wide %>%
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

##################
SI_matrix <- SI_only

# run the nMDS
SI_mds <- metaMDS(abs(SI_matrix))
# extract the 'points' from the nMDS that you will plot in ggplot2
SI_mds_points <- SI_mds$points
# turn those plot points into a dataframe that ggplot2 can read
SI_mds_points <- data.frame(SI_mds_points)
# join your plot points with your summed species observations from each habitat type
plot_data_tax <- data.frame(SI_mds_points, SI_tax[,c(7,8,69,70)])
plot_data_tax <- plot_data_tax %>%
  rename("division" = "phylum")


# run the ggplot
phylum_plot <- 
  ggplot(plot_data_tax, aes(x=MDS1, y=MDS2, 
                            fill = division)) +  
  labs(x = "nMDS1", y = "nMDS2") +
  theme_classic() + 
  geom_point(pch = 21, size = 2, color = "black") + 
  scale_fill_viridis(discrete = TRUE, begin = 0.2, end = 0.9, option = "G", name = "division") 
phylum_leg <- as_ggplot(get_legend(phylum_plot))

order_plot <- 
  ggplot(plot_data_tax, aes(x=MDS1, y=MDS2, # note pch and color are linked to *factors*
                            fill = order)) +  
  labs(x = "nMDS1", y = "nMDS2") +
  theme_classic() +
  geom_point(pch = 21, size = 2, color = "black") +  # set size of points to whatever you want
  guides(color=guide_legend(ncol=2)) +
  # geom_polygon(data=chulls_tax, aes(x=MDS1, y=MDS2, group=Habitat), fill=NA) + # optional 'hulls' around points
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "B", name = "order") # my favorite color-blind and b&w friendly palette, look at the viridis package for more details 
order_leg <- as_ggplot(get_legend(order_plot))


family_plot <- 
  ggplot(plot_data_tax, aes(x=MDS1, y=MDS2, # note pch and color are linked to *factors*
                            fill = family)) +  
  labs(x = "nMDS1", y = "nMDS2") +
  theme_classic() + 
  geom_point(pch = 21, size = 2, color = "black") +  # set size of points to whatever you want
  # geom_polygon(data=chulls_tax, aes(x=MDS1, y=MDS2, group=Habitat), fill=NA) + # optional 'hulls' around points
  guides(fill=guide_legend(ncol=2)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "H", name = "family") # my favorite color-blind and b&w friendly palette, look at the viridis package for more details 
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

annotate_figure(FigureMDS, top = text_grob("2D stress = 0.04", size = 10))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Online Resource Figure 3                                                     ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
################

# SI and reduced FA
### PERMANOVA 

# algal FA/SI for PERMANOVA
SIreduced_perm <- overlap_species %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`CN ratio`:`d13C`, `20:5w3`, `20:4w6`, `16:0`, `18:3w3`, `18:4w3c`, `18:1w9c`, `18:1w7c`, revisedSpecies)  %>%
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

# family
adonis2(abs(SIreduced_perm[,1:10]) ~ family, data = SIreduced_perm, method = "bray", na.rm = TRUE)
  # order
adonis2(abs(SIreduced_perm[,1:10]) ~ order, data = SIreduced_perm, method = "bray", na.rm = TRUE)


###################

# algal FA for adonis
marker_only <- overlap_species %>%
  select(`CN ratio`:`d13C`, `20:5w3`, `20:4w6`, `16:0`, `18:3w3`, `18:4w3c`, `18:1w9c`, `18:1w7c`) 

adonis2(abs(marker_only) ~ revisedSpecies, data = overlap_species, method = 'bray', na.rm = TRUE)
adonis2(abs(marker_only) ~ phylum, data = overlap_species, method = 'bray', na.rm = TRUE)



# run PCA
PCA_results <-  rda(marker_only, scale = TRUE)
# check that axes are above the mean (per Numerical Ecology)
#ev <- PCA_results$CA$eig
#ev>mean(ev)
# proportion explained
#barplot(ev, main="eigenvalues", col="bisque", las=2)
#abline(h=mean(ev), col="red")
#legend("topright", "Average eignenvalue", lwd=1, col=2)
# testing different scalings (1 = good for samples/sites, 2 = good for correlations)
#biplot(PCA_results, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
#biplot(PCA_results, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 
#autoplot(PCA_results, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 6, size = 4)


# extract PCA coordinates
uscores <- data.frame(PCA_results$CA$u)
uscores1 <- inner_join(rownames_to_column(filter(overlap_species, revisedSpecies != "Benthic diatoms")), rownames_to_column(data.frame(uscores)), by = "rowname")
vscores <- data.frame(PCA_results$CA$v)
vscores <- rownames_to_column(vscores)
vscores <- vscores 

# extract explanatory percentages
PCA_summary <- summary(PCA_results)
PCA_import <- as.data.frame(PCA_summary[["cont"]][["importance"]])
var_explained <- PCA_import[2, 1:2]

xvalues <- c(0.40, 0.08, 0.25, -0.35, 0.12, -0.39, 0.26, 0.46, 0.46, -0.25)
yvalues <- c(0.2, -0.02, -0.24, 0.24, -0.7, -0.23, 0.41, 0.15, 0.05, 0.49)
# make final ggplot figure
ggplot(uscores1) + 
  scale_fill_viridis(discrete = TRUE, guide = guide_legend(title = "species", label.theme = element_text(face = "italic"))) +
  scale_color_viridis(discrete = TRUE, guide = guide_legend(title = "species", label.theme = element_text(face = "italic")))  +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'grey30') +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_color_manual(values = rep("black", 31)) +
  #geom_point(aes(x = PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies,
  #               shape = phylum), size = 4) +
  geom_point(aes(x = PC1, y = PC2, fill = revisedSpecies, color = revisedSpecies, shape = phylum),
             size = 4) +
  geom_text_repel(data = vscores, aes(x = PC1, y = PC2, label = rowname), color = "red",
                  position = position_nudge_center(x = 0.02, y = 0.02, 0, 0)) +
  labs(shape = "division", fill = "species") +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  guides(fill = guide_legend(override.aes = list(shape=21), label.theme = element_text(face = "italic", size = 9)), color = "none") +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  xlim(-0.45, 0.5)


# size = 11x6 cairo PDF

####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#

# SCRATCH PAD ####

# get means and sd around SI values reported in text

# C at d13
meanvalues <- SI_wide %>%
  filter(revisedSpecies == "Callophyllis atrosanguinea") %>%
  select(d13C) 
mean(meanvalues$d13C)
sd(meanvalues$d13C)

# P ant d13
meanvalues <- SI_wide %>%
  filter(revisedSpecies == "Phyllophora antarctica") %>%
  select(d13C) 
mean(meanvalues$d13C)
sd(meanvalues$d13C)


# CN RHodo
meanvalues <- SI_wide %>%
  filter(phylum == "Rhodophyta") %>%
  select(`CN ratio`) 
mean(meanvalues$`CN ratio`)
sd(meanvalues$`CN ratio`)

# d15 RHodo
meanvalues <- SI_wide %>%
  filter(phylum == "Rhodophyta") %>%
  select(d15N) 
mean(meanvalues$d15N)
sd(meanvalues$d15N)

# d13 RHodo
meanvalues <- SI_wide %>%
  filter(phylum == "Rhodophyta") %>%
  select(d13C) 
mean(meanvalues$d13C)
sd(meanvalues$d13C)

# d15 Och
meanvalues <- SI_wide %>%
  filter(phylum == "Ochrophyta") %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(d15N) 
mean(meanvalues$d15N)
sd(meanvalues$d15N)

# d13 Och
meanvalues <- SI_wide %>%
  filter(phylum == "Ochrophyta") %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(d13C) 
mean(meanvalues$d13C)
sd(meanvalues$d13C)


# CN Och
meanvalues <- SI_wide %>%
  filter(phylum == "Ochrophyta") %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`CN ratio`) 
mean(meanvalues$`CN ratio`)
sd(meanvalues$`CN ratio`)

# d13 Chl
meanvalues <- SI_wide %>%
  filter(phylum == "Chlorophyta") %>%
  select(d13C) 
mean(meanvalues$d13C)
sd(meanvalues$d13C)

# n15 Chl
meanvalues <- SI_wide %>%
  filter(phylum == "Chlorophyta") %>%
  select(d15N) 
mean(meanvalues$d15N)
sd(meanvalues$d15N)



# get means and sd around FA values reported in text

# 16
meanvalues <- FA_wide %>%
  filter(phylum == "Rhodophyta") %>%
  select(`16:0`) 
mean(meanvalues$`16:0`)
sd(meanvalues$`16:0`)

# 20:5w3
meanvalues <- FA_wide %>%
  filter(phylum == "Rhodophyta") %>%
  select(`20:5w3`) 
mean(meanvalues$`20:5w3`)
sd(meanvalues$`20:5w3`)

# 20:4w6
meanvalues <- FA_wide %>%
  filter(phylum == "Rhodophyta") %>%
  select(`20:4w6`) 
mean(meanvalues$`20:4w6`)
sd(meanvalues$`20:4w6`)

# 16:1w7 new Rhodo only
meanvalues <- FA_wide %>%
  filter(phylum == "Rhodophyta") %>%
  filter(revisedSpecies %in% c("Ballia callitricha",
                               "Porphyra plocamiestris",
                               "Paraglossum salicifolium",
                               "Picconiella plumosa",
                               "Meridionella antarctica",
                               "Austropugetia crassa",
                               "Callophyllis atrosanguinea",
                               "Gymnogongrus antarcticus",
                               "Phyllophora antarctica",
                               "Pachymenia orbicularis",
                               "Trematocarpus antarcticus")) %>%
  select(`16:1w7c`) 
mean(meanvalues$`16:1w7c`)
sd(meanvalues$`16:1w7c`)

# OCHRO

# 20:4w6
meanvalues <- FA_wide %>%
  filter(phylum == "Ochrophyta") %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`20:4w6`) 
mean(meanvalues$`20:4w6`)
sd(meanvalues$`20:4w6`)

# 20:5w3
meanvalues <- FA_wide %>%
  filter(phylum == "Ochrophyta") %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`20:5w3`) 
mean(meanvalues$`20:5w3`)
sd(meanvalues$`20:5w3`)

# 16:0
meanvalues <- FA_wide %>%
  filter(phylum == "Ochrophyta") %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`16:0`) 
mean(meanvalues$`16:0`)
sd(meanvalues$`16:0`)

# 18:4w3
meanvalues <- FA_wide %>%
  filter(phylum == "Ochrophyta") %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`18:4w3c`) 
mean(meanvalues$`18:4w3c`)
sd(meanvalues$`18:4w3c`)

# 18:1w9
meanvalues <- FA_wide %>%
  filter(phylum == "Ochrophyta") %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`18:1w9c`) 
mean(meanvalues$`18:1w9c`)
sd(meanvalues$`18:1w9c`)

# CHLOR

# 16:0
meanvalues <- FA_wide %>%
  filter(phylum == "Chlorophyta") %>%
  select(`16:0`) 
mean(meanvalues$`16:0`)
sd(meanvalues$`16:0`)

# 18:3w3
meanvalues <- FA_wide %>%
  filter(phylum == "Chlorophyta") %>%
  select(`18:3w3`) 
mean(meanvalues$`18:3w3`)
sd(meanvalues$`18:3w3`)

# 18:2w6
meanvalues <- FA_wide %>%
  filter(phylum == "Chlorophyta") %>%
  select(`18:2w6c`) 
mean(meanvalues$`18:2w6c`)
sd(meanvalues$`18:2w6c`)



# Rerun of PERMANOVA for FA with site as a factor

### PERMANOVA 

# algal FA for adonis
FA_only <- FA_wide %>%
  filter(revisedSpecies != "Benthic diatoms") %>%
  select(`8:0`:`24:1w9`) 

# species only
adonis2(abs(FA_only) ~ revisedSpecies + siteName, data = filter(FA_wide, revisedSpecies != "Benthic diatoms"), method = 'bray', na.rm = TRUE)



# number of FA detected (mean) for each species
numberFA <- FASI_QAQC %>%
  select(revisedSpecies, `8:0`:`24:1w9`) %>%
  mutate(count=rowSums(.!=0)) %>%
  group_by(revisedSpecies) %>%
  reframe(min(count), max(count), mean(count))

# reviewer 2 16PUFA comment response: 

# Graeve: 18 FA ID
# 
# Georgiella confluens: 14.1%  16PUFA # my GECO likely has trace 16:3 & 16:4
# Neuroglossum ligulatum: 0.6%  16PUFA
# Pantoneura plocamioides: 13.1%  16PUFA # my PAPL likely has trace 16:3 & 16:4
# Gynogongrus turquetii: 0.1%  16PUFA
# Gigartina skottsbergii: 1.4%  16PUFA 
# Rhodymenia subantarctica: 2.5%  16PUFA
# Hymenocladiopsis crustigena: 2.5%  16PUFA
# 
# Desmarestia antarctica: 0.1%  16PUFA
# 
# Pereira: 25 FA ID
# 
# No ANT algae used 
# 
# Dictyota dichotoma: 0.44%  16PUFA
# 
# #########################
# 
# We detected from 21-36 across all samples with a mean ranging from 22-31 per species

# pull out P. plocameoides as example
pplo_us <- FASI_QAQC %>%
  filter(revisedSpecies == "Pantoneura plocamioides" & Transect == "T1") %>%
  select(revisedSpecies, `8:0`:`24:1w9`) %>%
  pivot_longer(`8:0`:`24:1w9`, names_to = "FA")


gcon_us <- FASI_QAQC %>%
  filter(revisedSpecies == "Georgiella confluens") %>%
  select(revisedSpecies, `8:0`:`24:1w9`) %>%
  pivot_longer(`8:0`:`24:1w9`, names_to = "FA") %>%
  group_by(FA) %>%
  summarise(mean(value, na.rm = TRUE))

# What species were collected at what sites?
collectionsites <- FASI_QAQC %>%
  filter(!is.na(siteName)) %>%
  select(revisedSpecies, siteName, `Latitude (dec)`) %>%
  replace_na(list(`Latitude (dec)` = 0)) %>%
  pivot_wider(id_cols = revisedSpecies, names_from = siteName ,values_from = `Latitude (dec)`, values_fn = max) %>%
  relocate("revisedSpecies", "A", "B", "C", "D", "E", "F", "G", "I", "J", "L", "M", "XX")

################# New Figure 5

# PCA


data=subset(mtcars, wt > 4 | mpg > 25)







# SI means per division and species grouping (for mean ranges)
simeans <- SI_values %>%
  group_by(phylum, revisedSpecies) %>%
  summarise(across(
    everything(),
    .fns = list(Mean = mean), na.rm = TRUE,
    .names = "{col}_{fn}"
  ))

SI_values %>%
  group_by(phylum) %>%
  summarise(mean(`CN ratio`))

### PERMANOVA 

# algal FA for adonis
D_only <- FA_wide %>%
  filter(revisedSpecies %in% c("Desmarestia menziesii", "Desmarestia anceps")) %>%
  select(`8:0`:`24:1w9`) 

adonis2(abs(D_only) ~ revisedSpecies, 
        data = filter(FA_wide, revisedSpecies %in% c("Desmarestia menziesii", "Desmarestia anceps")), 
        method = 'bray', na.rm = TRUE, perm = 9999)

# algal SI for adonis

D_only <- SI_wide %>%
  filter(revisedSpecies %in% c("Desmarestia menziesii", "Desmarestia anceps")) %>%
  select(`CN ratio`:`d13C`) 

adonis2(abs(D_only) ~ revisedSpecies, 
        data = filter(FA_wide, revisedSpecies %in% c("Desmarestia menziesii", "Desmarestia anceps")), 
        method = 'bray', na.rm = TRUE, perm = 9999)




# Phyllophora antarctica and Callophyllis atrosanguinea



### PERMANOVA 

# algal FA for adonis
PC_only <- FA_wide %>%
  filter(revisedSpecies %in% c("Phyllophora antarctica", "Callophyllis atrosanguinea")) %>%
  select(`8:0`:`24:1w9`) 

adonis2(abs(PC_only) ~ revisedSpecies, 
        data = filter(FA_wide, revisedSpecies %in% c("Phyllophora antarctica", "Callophyllis atrosanguinea")), 
                      method = 'bray', na.rm = TRUE, perm = 9999)

# algal SI for adonis

PC_only <- SI_wide %>%
  filter(revisedSpecies %in% c("Phyllophora antarctica", "Callophyllis atrosanguinea")) %>%
  select(`CN ratio`:`d13C`) 

adonis2(abs(PC_only) ~ revisedSpecies, 
        data = filter(FA_wide, revisedSpecies %in% c("Phyllophora antarctica", "Callophyllis atrosanguinea")), 
        method = 'bray', na.rm = TRUE, perm = 9999)


# ANOVA of differences between C:N ratios of divisions

library(DHARMa)
library(lme4)
library(rstatix)



# with trial as a random factor
CN_lm <- lm(log10(`CN ratio`) ~ phylum, data = SI_values)

# test if assumptions of model are met
CN_sim <- simulateResiduals(fittedModel = CN_lm, plot = F)
plot(CN_sim)
testDispersion(CN_lm)

CN_ANOVA <- anova(CN_lm)
summary(CN_lm)
CN_ANOVA

pwc <- SI_values %>% tukey_hsd(log10(`CN ratio`) ~ phylum)
pwc






####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#

# SCRATCH PAD ####

listFA <- colnames(FASI_QAQC)
list1 <- listFA[24:64]





############## JULIES REDO OF DATASET

# read in Julie's data and compare samples to existing FA list

























############## JULIES DATASET FOR TOTAL 16 PUFAs

JuliesAlgae <- read_csv("Data/Biomarkers/FattyAcids/JuliesAlgae.csv", 
                        col_types = cols(`C14:1w5` = col_double(), 
                                         `C16:2w6...8` = col_double(), `C16:2w6...9` = col_double(), 
                                         `C17:0` = col_double(), `C16:4w1` = col_double(), 
                                         `C16:3n4` = col_double(), `C18:3w3` = col_double(),
                                         `C18:5w3` = col_double(), `C18:4w3` = col_double(),
                                         `C18:1n` = col_double(), `C18:3w6` = col_double(), 
                                         `C20:4w3` = col_double(), `C22:0` = col_double(), 
                                         `C22:5w3` = col_double(), `C22:6w3` = col_double()))

JuliesAlgae$areaSum <- rowSums(JuliesAlgae[,2:29], na.rm = T)

proportions <- JuliesAlgae[,2:29]/JuliesAlgae$areaSum

proportions$Sample <- JuliesAlgae$Sample
proportions$species <- str_extract(JuliesAlgae$Sample, "^.{4}")

phyla <- proportions %>%
  mutate(phylum = case_when(species == "DEME" ~ "ochrophyta",
         species == "HIGR" ~ "ochrophyta",
         species == "IRCO" ~ "rhodophyta",
         species == "PHAN" ~ "rhodophyta",
         species == "PLCA" ~ "rhodophyta",
         species == "MYMA" ~ "rhodophyta",
         species == "ADUT" ~ "ochrophyta",
         species == "ASCR" ~ "rhodophyta",
         species == "ASMI" ~ "ochrophyta",
         species == "BACA" ~ "rhodophyta",
         species == "CAAT" ~ "rhodophyta",
         species == "CORA" ~ "rhodophyta",
         species == "CYJA" ~ "ochrophyta",
         species == "CYOB" ~ "rhodophyta",
         species == "DEAN" ~ "ochrophyta",
         species == "DEAP" ~ "ochrophyta",
         species == "DEPU" ~ "rhodophyta",
         species == "GECO" ~ "rhodophyta",
         species == "GISK" ~ "rhodophyta",
         species == "GYAN" ~ "rhodophyta",
         species == "HIGR" ~ "ochrophyta",
         species == "HYSP" ~ "rhodophyta",
         species == "LAAN" ~ "chlorophyta",
         species == "MYMA" ~ "rhodophyta",
         species == "MYSM" ~ "rhodophyta",
         species == "PADE" ~ "rhodophyta",
         species == "PAOR" ~ "rhodophyta",
         species == "PAPI" ~ "rhodophyta",
         species == "PASA" ~ "rhodophyta",
         species == "PIPI" ~ "rhodophyta",
         species == "POPL" ~ "rhodophyta",
         species == "RHLA" ~ "rhodophyta",
         species == "SYAU" ~ "rhodophyta",
         species == "PAPL" ~ "rhodophyta",
         species == "PIPL" ~ "rhodophyta",
         species == "TRAN" ~ "rhodophyta"))

PUFA16 <- phyla %>%
  filter(species != "BEDI") %>%
  filter(species != "ANTS") %>%
  select(Sample, species, phylum, `C16:2w6...8`, `C16:2w6...9`, `C16:3n4`, `C16:4w1`) %>%
  mutate(`C16:2w6` = `C16:2w6...8` + `C16:2w6...9`) %>%
  select(Sample, species, phylum, `C16:2w6`, `C16:3n4`, `C16:4w1`)

AllSummary <- PUFA16 %>%
  group_by(species, phylum) %>%
  summarise(across(starts_with("C"),  list(
    mean = ~ mean(.x, na.rm = T), 
    sd = ~ sd(.x, na.rm = T),
    max = ~ max(.x))))
  

speciesSummary <- AllSummary %>%
  replace(is.na(.), 0)

PhylumSummary <- PUFA16 %>%
  group_by(phylum) %>%
  summarise(across(starts_with("C"), list(
    mean = ~ mean(.x, na.rm = T), 
    sd = ~ sd(.x, na.rm = T),
    max = ~ max(.x, na.rm = T))))

phySummary <- PhylumSummary %>%
  replace(is.na(.), 0)

write_csv(phySummary, "C:/Users/Ross.Whippo/Desktop/phylum.csv")
