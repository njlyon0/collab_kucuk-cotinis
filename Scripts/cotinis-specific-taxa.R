## -------------------------------------------------------------- ##
          # Kucuk Cotinis Project - Specific Taxa Figures
## -------------------------------------------------------------- ##
# Code written by Nicholas J Lyon

# PURPOSE:
## Create figures for the specific microbial taxa of interest found in Cotinis nitida

# Clear environment (always better to start with tabula rasa)
rm(list = ls())

# Set working directory
getwd() # should end in ".../Kucuk-Cotinis-Collab"
myWD <- getwd()

# Necessary libraries
# devtools::install_github("NJLyon-Projects/helpR")
library(tidyverse); library(vegan); library(ape); library(helpR)
library(cowplot); library(ggvenn)

## -------------------------------------------------------------- ##
                # Data Retrieval & Housekeeping ####
## -------------------------------------------------------------- ##
# Retrieve the relevant datasets
beta <- read.csv(file.path("Data", "Tidy Data", "beta-diversity-data.csv"))

# Create directory to save these to
dir.create("Special Taxa Graphs", showWarnings = F)

# We also want to customize some plotting aesthetics for our ggplot plots that we can do here
all.cols <- c("Larval paunch" = "#a50026", "Larval ileum" = "#f46d43",
              "Larval midgut" = "#fee090",
              "Adult midgut" = "#abd9e9", "Adult hindgut" = "#4575b4")
sex.shps <- c("male" = 24, "female" = 25, "larva"= 21)
dodge <- position_dodge(width = 0.5)
pref_theme <- theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'none',
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 13, angle = 35, hjust = 1),
        legend.background = element_rect(fill = alpha('black', 0)))

## -------------------------------------------------------------- ##
                  # Specific Taxa Graphing ####
## -------------------------------------------------------------- ##
# We want several specific taxa to be graphed in the same way
## So let's for loop through this

# Create vector of desired taxa
desired_taxa <- c("Gilliamella", "Bacillus", "Turicibacter", "Geobacter",
                  "Alistipes", "Tannerella", "Bradyrhizobium", "Enterococcus",
                  "Diplosphaera", "Pelospora", "Methanobrevibacter", 
                  "Ruminococcaceae", "Dysgonomonadaceae", "Rikenellaceae",
                  "Desulfovibrionaceae", "Lachnospiraceae")

# Loop through that vector to summarize, plot, and export for each
for(k in 1:length(desired_taxa)){
  
  # Identify taxon
  special_taxon <- desired_taxa[k]
  
  # Summarize the dataframe as needed
  sub_df <- beta %>%
    # Pivot to long format
    pivot_longer(-Sample.ID:-Sex, names_to = 'Taxon', values_to = 'Abundance') %>%
    # Keep only the genus of interest
    filter(str_detect(string = Taxon, pattern = special_taxon)) %>%
    # Sum across samples within Stage.Gut
    group_by(Stage.Gut) %>%
    dplyr::summarise(Abundance = sum(Abundance, na.rm = T)) %>%
    # Re-level the Stage.Gut column
    dplyr::mutate(Stage.Gut = factor(Stage.Gut,
                                     levels = c("Larval paunch", "Larval ileum", "Larval midgut",
                                                "Adult midgut", "Adult hindgut"))) %>%
    # Return dataframe
    as.data.frame()
  
  # Assemble y-axis label
  if(str_detect(string = special_taxon, pattern = "aceae") == FALSE){
    y_axis_label <- paste(special_taxon, "spp. Abundance") }
  if(str_detect(string = special_taxon, pattern = "aceae") == TRUE){
    y_axis_label <- paste(special_taxon, "Abundance") }
  
  # Plot that data
  sub_plot <- ggplot(sub_df, aes(y = Abundance, x = Stage.Gut, fill = Stage.Gut, color = 'x')) +
    geom_bar(stat = 'identity') +
    scale_fill_manual(values = all.cols) +
    scale_color_manual(values = 'black') +
    labs(x = "Life Stage & Gut Region", y = y_axis_label) +
    pref_theme
  
  # Save that plot
  ggplot2::ggsave(file.path("Special Taxa Graphs",
                            paste0(special_taxon, "-Figure.tiff")),
                  width = 4, height = 4, plot = sub_plot)
  
  # Alert the user which graph was produced
  message(special_taxon, " graphic produced. ", length(desired_taxa) - k, " remaining." )
  
}

# End ---
