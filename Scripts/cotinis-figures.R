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
library(Rmisc); library(cowplot)

## -------------------------------------------- ##
      # Data Retrieval & Housekeeping ####
## -------------------------------------------- ##
# Retrieve the relevant datasets
alpha <- read.csv("./Data/Tidy Data/alpha-diversity-data.csv")
fams <- read.csv("Data/Tidy Data/family-abun.csv")
phyla <- read.csv("Data/Tidy Data/phylum-abun.csv")

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

## -------------------------------------------- ##
  # Figure 1 - Alpha Diversity Superfigure ####
## -------------------------------------------- ##
# Summarize each of the needed response variables
chao <- summarySE(data = alpha, measurevar = "Chao1", groupvars = c("Stage.Gut"))
simp <- summarySE(data = alpha, measurevar = "Simpson", groupvars = c("Stage.Gut"))
ace <- summarySE(data = alpha, measurevar = "ACE", groupvars = c("Stage.Gut"))

# Re-level the factor column
  ## Pref factor order
stage.gut.lvls <- c("Adult hindgut", "Adult midgut", "Larval midgut",
                    "Larval ileum", "Larval paunch")
  ## Actual re-leveling
chao$Stage.Gut <- factor(chao$Stage.Gut, levels = stage.gut.lvls)
simp$Stage.Gut <- factor(simp$Stage.Gut, levels = stage.gut.lvls)
ace$Stage.Gut <- factor(ace$Stage.Gut, levels = stage.gut.lvls)

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
  pref_theme + theme(axis.ticks.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.title.x = element_blank())

# ACE Diversity
ace.plt <- ggplot(ace, aes(y = ACE, x = Stage.Gut,
                           fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = ACE + se, ymin = ACE - se), width = 0.2) +
  geom_text(label = "ab", x = 0.7, y = 180, size = 5) +
  geom_text(label = "a", x = 1.7, y =  100, size = 5) +
  geom_text(label = "bc", x = 2.7, y = 740, size = 5) +
  geom_text(label = "b", x = 3.7, y =  590, size = 5) +
  geom_text(label = "c", x = 4.7, y =  975, size = 5) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "ACE") +
  pref_theme + theme(legend.position = 'none')

# Examine your handiwork
chao.plt
simp.plt
ace.plt

# Create the superfigure
plot_grid(chao.plt, simp.plt, ace.plt,
          ncol = 1, nrow = 3, 
          labels = c("A", "B", "C"))

# Save it
ggplot2::ggsave("./Figures/Alpha-Diversity-Superfigure.pdf",
                width = 6.5, height = 9,
                plot = last_plot())

## -------------------------------------------- ##
   # Figure 2 - Beta Diversity Superfigure ####
## -------------------------------------------- ##
# Made in PowerPoint (sorry!)

## -------------------------------------------- ##
     # Figure 3 - Abundance Superfigure ####
## -------------------------------------------- ##
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
ggplot2::ggsave("./Figures/Abundance-Superfigure.pdf",
                width = 5, height = 4,
                plot = last_plot())

## -------------------------------------------- ##
      # Figure 4 - Relative Abundance ####
## -------------------------------------------- ##
## ----------------------------- ##
      # Fig 4 Housekeeping ####
## ----------------------------- ##
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
                    yes = paste0("Family <", abun.thresh, "% Total"),
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
                    yes = paste0("Phyla <", abun.thresh, "% Total"),
                    no = Phylum)
  ) %>%
  group_by(Sample.ID, Phylum) %>%
  dplyr::summarise(relativeAbun = sum(relativeAbun, na.rm = T)) %>%
  as.data.frame()

# Check the ranges to see that worked
range(fams.v2$relativeAbun)
range(phyla.v2$relativeAbun)

## ----------------------------- ##
      # Fig 4 Graphing ####
## ----------------------------- ##
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

# Set colors for each phylum
phyla.cols <- c("white", "#ffff99", "#6a3d9a",
                "#cab2d6","#ff7f00","#fdbf6f", 
                "#e31a1c", "#fb9a99", "#33a02c", 
                "#b2df8a", "#1f78b4", "#a6cee3")

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
                               "Paunch20", "Paunch21"
                             ))

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
        axis.text.x = element_text(size = 9))
  
# Save this!
ggplot2::ggsave("./Figures/Relative-Abundance-Superfigure.pdf",
                width = 11, height = 5,
                plot = last_plot())


# END ####

