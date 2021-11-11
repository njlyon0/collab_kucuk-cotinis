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

# END ####

