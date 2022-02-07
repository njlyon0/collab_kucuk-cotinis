## -------------------------------------------------------------- ##
         # Kucuk Cotinis Project - Figure Creation Code
## -------------------------------------------------------------- ##
# Code written by Nicholas J Lyon

# PURPOSE:
## Create figures for the gut microbe communities found in Cotinis nitida

# Clear environment (always better to start with tabula rasa)
rm(list = ls())

# Set working directory
getwd() # should end in ".../Kucuk-CotinisCollab"
myWD <- getwd()

# Necessary libraries
library(tidyverse); library(vegan); library(ape)
library(Rmisc); library(cowplot); library(ggvenn)

## -------------------------------------------------------------- ##
              # Data Retrieval & Housekeeping ####
## -------------------------------------------------------------- ##
# Retrieve the relevant datasets
alpha <- read.csv("./Data/Tidy Data/alpha-diversity-data.csv")
fams <- read.csv("Data/Tidy Data/family-abun.csv")
phyla <- read.csv("Data/Tidy Data/phylum-abun.csv")
beta <- read.csv("./Data/Tidy Data/beta-diversity-data.csv")

# Look at them to be sure nothing obvious is wrong
str(alpha)
str(fams)
str(phyla)

# We also want to customize some plotting aesthetics for our ggplot plots that we can do here
all.cols <- c("Larval paunch" = "#a50026", "Larval ileum" = "#f46d43", "Larval midgut" = "#fee090",
              "Adult midgut" = "#abd9e9", "Adult hindgut" = "#4575b4")
sex.shps <- c("male" = 24, "female" = 25, "larva"= 21)
dodge <- position_dodge(width = 0.5)
pref_theme <- theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'none',
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        legend.background = element_rect(fill = alpha('black', 0)),
        axis.text.x = element_text(angle = 35, hjust = 1))

## -------------------------------------------------------------- ##
          # Figure 1 - Alpha Diversity Superfigure ####
## -------------------------------------------------------------- ##
# Summarize each of the needed response variables
chao <- summarySE(data = alpha, measurevar = "Chao1", groupvars = c("Stage.Gut"))
simp <- summarySE(data = alpha, measurevar = "Simpson", groupvars = c("Stage.Gut"))
ace <- summarySE(data = alpha, measurevar = "ACE", groupvars = c("Stage.Gut"))
shan <- summarySE(data = alpha, measurevar = "Shannon", groupvars = c("Stage.Gut"))

# Re-level the factor column
  ## Pref factor order
stage.gut.lvls <- c("Adult hindgut", "Adult midgut", "Larval midgut",
                    "Larval ileum", "Larval paunch")
  ## Actual re-leveling
chao$Stage.Gut <- factor(chao$Stage.Gut, levels = stage.gut.lvls)
simp$Stage.Gut <- factor(simp$Stage.Gut, levels = stage.gut.lvls)
ace$Stage.Gut <- factor(ace$Stage.Gut, levels = stage.gut.lvls)
shan$Stage.Gut <- factor(shan$Stage.Gut, levels = stage.gut.lvls)

# Chao 1 Diversity
chao.plt <- ggplot(chao, aes(y = Chao1, x = Stage.Gut,
                             fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Chao1 + se, ymin = Chao1 - se), width = 0.2) +
  geom_text(label = "ab", x = 0.7, y = 180, size = 5) +
  geom_text(label = "a", x = 1.7, y =  100, size = 5) +
  geom_text(label = "bc", x = 2.7, y = 740, size = 5) +
  geom_text(label = "b", x = 3.7, y =  590, size = 5) +
  geom_text(label = "c", x = 4.7, y =  975, size = 5) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "Chao1 Index") +
  pref_theme + theme(axis.ticks.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.title.x = element_blank())

# Simpson Diversity
simp.plt <- ggplot(simp, aes(y = Simpson, x = Stage.Gut,
                             fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Simpson + se, ymin = Simpson - se), width = 0.2) +
  geom_text(label = "b", x = 0.7, y =  0.99, size = 5) +
  geom_text(label = "a", x = 1.7, y =  0.82, size = 5) +
  geom_text(label = "b", x = 2.7, y =  1.05, size = 5) +
  geom_text(label = "b", x = 3.7, y =  1, size = 5) +
  geom_text(label = "b", x = 4.7, y =  1.025, size = 5) +
  ylim(0, 1.1) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "Simpson Diversity") +
  pref_theme

# ACE Diversity
ace.plt <- ggplot(ace, aes(y = ACE, x = Stage.Gut,
                           fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = ACE + se, ymin = ACE - se), width = 0.2) +
  geom_text(label = "ab", x = 0.7, y = 180, size = 5) +
  geom_text(label = "a", x = 1.7, y =  100, size = 5) +
  geom_text(label = "bc", x = 2.7, y = 750, size = 5) +
  geom_text(label = "b", x = 3.7, y =  600, size = 5) +
  geom_text(label = "c", x = 4.7, y =  985, size = 5) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "ACE") +
  pref_theme + theme(axis.ticks.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.title.x = element_blank())

shan.plt <- ggplot(shan, aes(y = Shannon, x = Stage.Gut,
                             fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Shannon + se, ymin = Shannon - se), width = 0.2) +
  geom_text(label = "b", x = 0.7, y =  3.8, size = 5) +
  geom_text(label = "a", x = 1.7, y =  2.3, size = 5) +
  geom_text(label = "c", x = 2.7, y =  6.15, size = 5) +
  geom_text(label = "bc", x = 3.7, y = 4.75, size = 5) +
  geom_text(label = "c", x = 4.7, y =  5.45, size = 5) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "Shannon Diversity") +
  pref_theme + theme(legend.position = 'none')

# Examine your handiwork
chao.plt
simp.plt
ace.plt
shan.plt

# Create the superfigure
plot_grid(chao.plt, ace.plt,
          simp.plt, shan.plt,
          align = 'v',
          rel_heights = c(0.7, 1),
          ncol = 2, nrow = 2, 
          labels = c("A", "B", "C", "D"))

# Save it
ggplot2::ggsave("./Figures/Alpha-Diversity-Superfigure.tiff",
                device = 'tiff',
                width = 7, height = 7,
                plot = last_plot())

## -------------------------------------------------------------- ##
            # Figure 2 - Beta Diversity Superfigure ####
## -------------------------------------------------------------- ##
# Made in PowerPoint (sorry!)
  ## Individual graphs within figure found in "cotinis-data-viz.R" script

## -------------------------------------------------------------- ##
            # Figure 3 - Abundance Superfigure ####
## -------------------------------------------------------------- ##
# Look at the phylum- & family-level abundance
head(phyla)
head(fams)

# Calculate truly global abundance for both
phy.total.abun <- sum(phyla$Abundance, na.rm = T)
fam.total.abun <- sum(fams$Abundance, na.rm = T)

# What is the threshold of relative abundance we want?
threshold <- 0.05

# What is the relative abundance of that threshold?
phy.thresh <- phy.total.abun * threshold
fam.thresh <- fam.total.abun * threshold

# Get a global version of phylum-level abundance
  ## (i.e., irrespective of SampleID)
phy.global <- phyla %>%
  # Remove un-needed columns
  dplyr::select(Domain, Phylum, Abundance) %>%
  # Group by phylum
  group_by(Phylum) %>%
  # Sum to get one instance of each
  dplyr::summarise(
    Abundance = sum(Abundance, na.rm = T)
  ) %>%
  # Remove any taxa under the threshold value
  filter(Abundance > phy.thresh) %>%
  # Return a dataframe
  as.data.frame()

# Do the same for family-level abundance
fam.global <- fams %>%
  dplyr::select(Domain, Family, Abundance) %>%
  group_by(Family) %>%
  dplyr::summarise(Abundance = sum(Abundance, na.rm = T)) %>%
  filter(Abundance > fam.thresh) %>%
  as.data.frame()

# Look at what that produced
head(phy.global)
head(fam.global)

# Graph phylum level abundance
phy.plt <- ggplot(phy.global, aes(y = Abundance,
                                  x = reorder(Phylum, -Abundance),
                                  fill = 'x', color = 'x')) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = "#980043") +
  scale_color_manual(values = 'black') +
  labs(x = "Phylum", y = "Abundance") +
  pref_theme + theme(legend.position = 'none',
                     axis.text.x = element_text(size = 8))

# Graph family-level abundance
fam.plt <- ggplot(fam.global, aes(y = Abundance,
                                  x = reorder(Family, -Abundance),
                                  fill = 'x', color = 'x')) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = "#e7298a") +
  scale_color_manual(values = 'black') +
  labs(x = "Family", y = "Abundance") +
  pref_theme + theme(legend.position = 'none',
                     axis.text.x = element_text(size = 8))

# Examine both
phy.plt
fam.plt

# Create combination graph
plot_grid(phy.plt, fam.plt, ncol = 2, nrow = 1, labels = c("A", "B"))

# Save it
ggplot2::ggsave("./Figures/Abundance-Superfigure.tiff",
                device = 'tiff',
                width = 5, height = 4,
                plot = last_plot())

## -------------------------------------------------------------- ##
                # Figure 4 - Relative Abundance ####
## -------------------------------------------------------------- ##
## -------------------------------------------- ##
              # F4 Housekeeping ####
## -------------------------------------------- ##
# Look at the range of relative abundances
  ## Visually
psych::multi.hist(fams$relativeAbun)
psych::multi.hist(phyla$relativeAbun)
  ## Numerically
range(fams$relativeAbun)
range(phyla$relativeAbun)

# Identify the threshold value of "too low"
abun.thresh <- 5

# We have a good alternative to dropping low abundances altogether:
fams.v2 <- fams %>%
  # Remove unneeded columns
  dplyr::select(Sample.ID, Family, relativeAbun) %>%
  dplyr::mutate(
    # Change the family to "other" if it's less than the cutoff
    Family = ifelse(relativeAbun < abun.thresh,
                    yes = paste0("Family < ", abun.thresh, "% Total"),
                    no = Family)
    ) %>%
  # Then, within sample,
  group_by(Sample.ID, Family) %>%
  # Sum through the relative abundances
  dplyr::summarise(
    relativeAbun = sum(relativeAbun, na.rm = T)
  ) %>%
  as.data.frame()

# Do the same for phylum
phyla.v2 <- phyla %>%
  dplyr::select(Sample.ID, Phylum, relativeAbun) %>%
  dplyr::mutate(
    Phylum = ifelse(relativeAbun < abun.thresh,
                    yes = paste0("Phyla < ", abun.thresh, "% Total"),
                    no = Phylum)
  ) %>%
  group_by(Sample.ID, Phylum) %>%
  dplyr::summarise(relativeAbun = sum(relativeAbun, na.rm = T)) %>%
  as.data.frame()

# Check the ranges to see that worked
range(fams.v2$relativeAbun)
range(phyla.v2$relativeAbun)

## -------------------------------------------- ##
               # F4 Graphing ####
## -------------------------------------------- ##
# Family-level *relative* abundance figure
ggplot(fams.v2, aes(x = Sample.ID, y = relativeAbun,
                    fill = Family, color = 'x',
                    order = -relativeAbun)) +
  geom_bar(stat = 'identity', position = 'stack', width = 0.65) +
  scale_color_manual(values = 'black') +
  labs(x = 'Sample ID', y = "Relative Abundance (%)") +
  #scale_fill_manual(values = sp::bpy.colors()) +
  pref_theme + theme(axis.title.x = element_blank())
  #theme(legend.position = 'right', axis.text.x = element_text(size = 8))

# Re-level the Sample.ID factor
phyla.v2$Sample.ID <- factor(phyla.v2$Sample.ID,
                             levels = c(
                               # Adult midgut
                               "Amid7", "Amid12", "Amid39", "Amid40",
                               "Amid41", "Amid46", "Amid48", "Amid49",
                               # Adult hindgut
                               "hind4", "hind7", "hind12", "hind39",
                               "hind40", "hind41", "hind46", "hind48",
                               "hind49",  
                               # Larval midgut
                               "Lmid2", "Lmid3", "Lmid5", "Lmid9",
                               # Larval ileum
                               "Ileum2", "Ileum3", "Ileum5", "Ileum7", "Ileum9",
                               "Ileum12", "Ileum17", "Ileum18", "Ileum19",
                               "Ileum20", "Ileum21",
                               # Larval paunch
                               "Paunch2", "Paunch3", "Paunch5", "Paunch9",
                               "Paunch12", "Paunch17", "Paunch18", "Paunch19",
                               "Paunch20", "Paunch21"))

# Set colors for each phylum
phyla.cols <- c("white", "#ffff99", "#6a3d9a",
                "#cab2d6","#ff7f00","#fdbf6f", 
                "#e31a1c", "#fb9a99", "#33a02c", 
                "#b2df8a", "#1f78b4", "#a6cee3")

# Phylum-level relative abundance
ggplot(phyla.v2, aes(x = Sample.ID, y = relativeAbun,
                     fill = reorder(Phylum, relativeAbun),
                     color = 'x')) +
  geom_bar(stat = 'identity', position = 'stack') +
  scale_color_manual(values = 'black') +
  labs(x = 'Sample ID', y = "Relative Abundance (%)") +
  scale_fill_manual(values = phyla.cols) +
  # Add vertical lines separating gut/life stages
  geom_vline(xintercept = 8.5) +
  geom_vline(xintercept = 17.5) +
  geom_vline(xintercept = 21.5) +
  geom_vline(xintercept = 32.5) +
  # Add text to top of graph defining each chunk (delinated by vertical lines)
  geom_text(label = "Adult Midgut", hjust = 'center', x = 4.25, y = 103) +
  geom_text(label = "Adult Hindgut", hjust = 'center', x = 13, y = 103) +
  geom_text(label = "Larval Mid.", hjust = 'center', x = 19.5, y = 103) +
  geom_text(label = "Ileum", hjust = 'center', x = 27, y = 103) +
  geom_text(label = "Paunch", hjust = 'center', x = 37.5, y = 103) +
  # Remaining aesthetics things
  pref_theme + guides(color = 'none') +
  theme(legend.position = 'right',
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 9),
        legend.key = element_rect(color = "black"))
  
# Save this!
ggplot2::ggsave("./Figures/Relative-Abundance-Superfigure.tiff",
                device = 'tiff',
                width = 11, height = 5,
                plot = last_plot())

## -------------------------------------------------------------- ##
   # Figure 5 - Within Life/Gut Phyla Abundance Superfigure ####
## -------------------------------------------------------------- ##
## -------------------------------------------- ##
             # F5 Housekeeping ####
## -------------------------------------------- ##
# Want to exclude very low abundance phyla to clarify the figures

# Re-identify abundance threshold (%)
new.abun.thresh <- 5

# Create subsets with relative abundance threshold implemented
amid.phyl <- phyla %>%
  dplyr::mutate(
    Phylum = ifelse(relativeAbun < new.abun.thresh,
                    yes = paste0("Phyla < ", new.abun.thresh, "% Total"),
                    no = Phylum)
  ) %>%
  filter(Stage.Gut == "Adult midgut") %>%
  group_by(Stage.Gut, Phylum) %>%
  summarySE(measurevar = "Abundance", groupvars = c("Stage.Gut", "Phylum")) %>%
  as.data.frame()

ahind.phyl <- phyla %>%
  dplyr::mutate(
    Phylum = ifelse(relativeAbun < new.abun.thresh,
                    yes = paste0("Phyla < ", new.abun.thresh, "% Total"),
                    no = Phylum)
  ) %>%
  filter(Stage.Gut == "Adult hindgut") %>%
  group_by(Stage.Gut, Phylum) %>%
  summarySE(measurevar = "Abundance", groupvars = c("Stage.Gut", "Phylum")) %>%
  as.data.frame()

lmid.phyl <- phyla %>%
  dplyr::mutate(
    Phylum = ifelse(relativeAbun < new.abun.thresh,
                    yes = paste0("Phyla < ", new.abun.thresh, "% Total"),
                    no = Phylum)
  ) %>%
  filter(Stage.Gut == "Larval midgut") %>%
  group_by(Stage.Gut, Phylum) %>%
  summarySE(measurevar = "Abundance", groupvars = c("Stage.Gut", "Phylum")) %>%
  as.data.frame()

ileum.phyl <- phyla %>%
  dplyr::mutate(
    Phylum = ifelse(relativeAbun < new.abun.thresh,
                    yes = paste0("Phyla < ", new.abun.thresh, "% Total"),
                    no = Phylum)
  ) %>%
  filter(Stage.Gut == "Larval ileum") %>%
  group_by(Stage.Gut, Phylum) %>%
  summarySE(measurevar = "Abundance", groupvars = c("Stage.Gut", "Phylum")) %>%
  as.data.frame()

paunch.phyl <- phyla %>%
  dplyr::mutate(
    Phylum = ifelse(relativeAbun < new.abun.thresh,
                    yes = paste0("Phyla < ", new.abun.thresh, "% Total"),
                    no = Phylum)
  ) %>%
  filter(Stage.Gut == "Larval paunch") %>%
  group_by(Stage.Gut, Phylum) %>%
  summarySE(measurevar = "Abundance", groupvars = c("Stage.Gut", "Phylum")) %>%
  as.data.frame()

## -------------------------------------------- ##
     # Fig 5 Graphing - Low Abun Excl. ####
## -------------------------------------------- ##
# Make component graphs of phylum-level abundance within each life stage/gut bit
## Adult midgut
amid <- ggplot(amid.phyl, aes(y = Abundance, x = reorder(Phylum, -Abundance),
                              fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  geom_text(label = "Adult Midgut", hjust = 'center',
            x = 3, y = 1475) +
  labs(x = "Phylum", y = "Abundance") +
  pref_theme + theme(legend.position = 'none',
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 9),
                     axis.title.y = element_text(size = 12),
                     axis.title.x = element_blank())

## Adult hindgut
ahind <- ggplot(ahind.phyl, aes(y = Abundance, x = reorder(Phylum, -Abundance),
                                fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  geom_text(label = "Adult Hindgut", hjust = 'center',
            x = 3, y = 8500) +
  labs(x = "Phylum", y = "Abundance") +
  pref_theme + theme(legend.position = 'none',
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 9),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank())

## Larval midgut
lmid <- ggplot(lmid.phyl, aes(y = Abundance, x = reorder(Phylum, -Abundance),
                              fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  geom_text(label = "Larval Midgut", hjust = 'center',
            x = 3.5, y = 1750) +
  labs(x = "Phyla", y = "Abundance") +
  pref_theme + theme(legend.position = 'none',
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 9),
                     axis.title = element_text(size = 12))

## Larval ileum
lile <- ggplot(ileum.phyl, aes(y = Abundance, x = reorder(Phylum, -Abundance),
                               fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  geom_text(label = "Larval Ileum", hjust = 'center',
            x = 5, y = 4150) +
  ylim(c(0, 4250)) +
  labs(x = "Phyla", y = "Abundance") +
  pref_theme + theme(legend.position = 'none',
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 9),
                     axis.title = element_text(size = 12),
                     axis.title.y = element_blank())

## Larval paunch
lpau <- ggplot(paunch.phyl, aes(y = Abundance, x = reorder(Phylum, -Abundance),
                                fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  geom_text(label = "Larval Paunch", hjust = 'center',
            x = 3, y = 13000) +
  labs(x = "Phyla", y = "Abundance") +
  pref_theme + theme(legend.position = 'none',
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 9),
                     axis.title = element_text(size = 12),
                     axis.title.y = element_blank())

# Examine each
amid
ahind
lmid
lile
lpau

# Assemble into larger composite graph
plot_grid(amid, ahind, NA, lmid, lile, lpau,
          ncol = 3, nrow = 2,
          labels = c("A", "B", NA, "C", "D", "E"))

# Option A:
ggplot2::ggsave("./Figures/Relative-Abundance-Subsets-Superfigure.tiff",
                device = 'tiff',
                width = 11, height = 6.5,
                plot = last_plot())


## -------------------------------------------------------------- ##
                   # Figure 6 - Venn Diagram ####
## -------------------------------------------------------------- ##
# Look at beta diversity data
str(beta[1:20])

# Prepare it for use in the venn diagram
venn.df <- beta %>%
  # Pivot to long format
  pivot_longer(
    -Sample.ID:-Sex,
    names_to = 'Taxon',
    values_to = 'Abundance'
  ) %>%
  # Group by Stage.Gut and Taxon
  dplyr::group_by(Stage.Gut, Taxon) %>%
  # Sum across samples within the aforementioned groups
  dplyr::summarise(Abundance = sum(Abundance, na.rm = T)) %>%
  # Do some changes to work with the `ggvenn` function/package
  dplyr::mutate(
    # Change NA abundance to FALSE
    Abundance = ifelse(is.na(Abundance),
                       yes = F, no = Abundance),
    # Change 0 abundance to FALSE
    Abundance = ifelse(Abundance == 0,
                       yes = F, no = Abundance),
    # Change any non-NA and non-zero numbers to TRUE
    Abundance = ifelse(Abundance > 0,
                       yes = T, no = Abundance)
  ) %>%
  # Change 0/1 to F/T
  dplyr::mutate(
    Abundance = ifelse(Abundance == 0,
                       yes = FALSE, no = TRUE)
  ) %>%
  # Pivot to wide format with Stage.Gut as columns
  pivot_wider(
    id_cols = Taxon,
    names_from = Stage.Gut,
    values_from = Abundance
  )

# Check it out again
str(venn.df)

# Make Venn diagram for adults
ggvenn(data = venn.df,
       columns = c("Adult hindgut", "Adult midgut"),
       fill_color = c('blue', 'light blue'),
       text_size = 6, set_name_size = 4.5,
       show_percentage = F)

# Save it
ggplot2::ggsave("./Figures/Venn-Diagram-ADULT.tiff",
                device = 'tiff',
                width = 4, height = 4,
                plot = last_plot())

# Do the same for larvae
  ## Plot
ggvenn(data = venn.df,
       columns = c("Larval paunch", "Larval ileum", "Larval midgut"),
       fill_color = c('red', 'orange', 'yellow'),
       text_size = 5, set_name_size = 4.5,
       show_percentage = F)
  ## Saver
ggplot2::ggsave("./Figures/Venn-Diagram-LARVAE.tiff",
                device = 'tiff',
                width = 4, height = 4,
                plot = last_plot())

## -------------------------------------------------------------- ##
                # Figure 7 - Gilliamella spp. ####
## -------------------------------------------------------------- ##
# Need to process the beta diversity data again (but differently here)
gill.df <- beta %>%
  # Pivot to long format
  pivot_longer(
    -Sample.ID:-Sex,
    names_to = 'Taxon',
    values_to = 'Abundance'
  ) %>%
  # Keep only the genus of interest
  filter(str_detect(Taxon, 'Gilliamella')) %>%
  # Sum across samples within Stage.Gut
  group_by(Stage.Gut) %>%
  dplyr::summarise(Abundance = sum(Abundance, na.rm = T)) %>%
  # Re-level the Stage.Gut column
  dplyr::mutate(
    Stage.Gut = factor(Stage.Gut, levels = stage.gut.lvls)
  ) %>%
  as.data.frame()

# Gilliamella plot
ggplot(gill.df, aes(y = Abundance, x = Stage.Gut,
                             fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region",
       y = expression(paste(italic("Gilliamella"), " spp. Abundance"))) +
  pref_theme

# Save it
ggplot2::ggsave("./Figures/Gilliamella-Figure.tiff",
                device = 'tiff',
                width = 4, height = 4,
                plot = last_plot())


# END ####

