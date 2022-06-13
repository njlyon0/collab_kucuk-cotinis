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
# devtools::install_github("njlyon0/helpR")
library(tidyverse); library(vegan); library(ape); library(helpR)
library(cowplot); library(ggvenn)

## -------------------------------------------------------------- ##
                # Data Retrieval & Housekeeping ####
## -------------------------------------------------------------- ##
# Retrieve the relevant datasets
beta <- read.csv(file.path("Data", "Tidy Data", "beta-diversity-data.csv"))

# Create directory to save these to
dir.create("Special Taxa Graphs", showWarnings = F)

# We also want to make a standard ggplot2 theme
pref_theme <- theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'none',
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 13, angle = 35, hjust = 1),
        legend.background = element_rect(fill = alpha('black', 0)))

# Identify the desired order of gut regions for the factor levels
gut_levels <- c("Larval paunch", "Larval ileum", "Larval midgut",
                "Adult midgut", "Adult hindgut")

# Create vector for coloring as well
## Order is one white then order of 'gut_levels' object above
color_vec <- c("#a50026", "#f46d43", "#fee090", "#abd9e9", "#4575b4")

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

# Suppress how chatty `dplyr::summarise()` wants to be
options(dplyr.summarise.inform = FALSE)

# Create an empty list for storing plots
plot_list <- list()

# Loop through that vector to summarize, plot, and export for each
for(k in 1:length(desired_taxa)){
# for(k in 1){ #retained for ease of modifying loop in future
  
  # Identify taxon
  special_taxon <- desired_taxa[k]
  
  # Summarize the dataframe as needed
  sub_df <- beta %>%
    # Pivot to long format
    tidyr::pivot_longer(-Sample.ID:-Sex, names_to = 'Taxon',
                        values_to = 'Abundance') %>%
    # Coerce all IDs other than the desired one into 'other'
    dplyr::mutate(Taxon_Simp = dplyr::case_when(
      stringr::str_detect(string = Taxon, pattern = special_taxon) ~ special_taxon,
      TRUE ~ "Other") ) %>%
    # Keep only the genus of interest
    # filter() %>%
    # Sum across samples within Stage.Gut
    group_by(Stage.Gut, Taxon_Simp) %>%
    dplyr::summarise(Abundance = sum(Abundance, na.rm = T)) %>%
    # Re-level the Stage.Gut column
    dplyr::mutate(Stage.Gut = factor(Stage.Gut, levels = gut_levels)) %>%
    # Calculate relative abundance
    dplyr::group_by(Stage.Gut) %>%
    dplyr::mutate(
      total = sum(Abundance, na.rm = T),
      percAbun = (Abundance / total) * 100) %>%
    # Drop the 'Other' taxa
    dplyr::filter(Taxon_Simp != "Other") %>%
    # Return dataframe
    as.data.frame()
  
 # Assemble y-axis label
  if(str_detect(string = special_taxon, pattern = "aceae") == FALSE){
    y_axis_label <- paste(special_taxon, "spp. (%)") }
  if(str_detect(string = special_taxon, pattern = "aceae") == TRUE){
    y_axis_label <- paste(special_taxon, "(%)") }

  # Plot that subset of the data
  sub_plot <- ggplot(sub_df, aes(y = percAbun, x = Stage.Gut,
                                 fill = Stage.Gut, color = 'x')) +
    geom_bar(stat = 'identity') +
    scale_fill_manual(values = color_vec,
                      labels = c("Other Taxa", special_taxon)) +
    scale_color_manual(values = 'black') +
    labs(x = "Life Stage & Gut Region", y = y_axis_label) +
    ylim(0, 100) +
    pref_theme
  
  # Remove x-axis labels and ticks for all but what will be the bottom row
  if(k < 13){
    sub_plot <- sub_plot + theme(axis.title.x = element_blank(),
                                 axis.text.x = element_blank())
  }
  
  # Save that plot
  ggplot2::ggsave(file.path("Special Taxa Graphs",
                            paste0(special_taxon, "-Figure.tiff")),
                  width = 4, height = 5, plot = sub_plot)
  
  # Also add each plot to a list before continuing
  plot_list[[k]] <- sub_plot
  
  # Alert the user which graph was produced
  message(special_taxon, " graphic produced. ", length(desired_taxa) - k, " remaining." )
  
}

## -------------------------------------------------------------- ##
                # Combination Figure Creation ####
## -------------------------------------------------------------- ##
# Now we want to make a single grid of all these graphs
length(plot_list)

# Let's make the figure
cowplot::plot_grid(plotlist = plot_list, nrow = 4, ncol = 4,
                   labels = 'AUTO', label_x = 0.2)

# And save that
ggplot2::ggsave(file.path("Figures",
                          "Special-Taxa-Superfigure.tiff"),
                width = 12, height = 12, plot = last_plot())

# End ---
