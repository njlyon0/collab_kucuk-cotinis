## -------------------------------------------------------------- ##
# Kucuk Cotinis Project - Summary Code
## -------------------------------------------------------------- ##
# Code written by Nicholas J Lyon

# PURPOSE:
## Get some useful summaries of the data

# Necessary libraries
# install.packages("librarian")
librarian::shelf(tidyverse, supportR)

# Clear environment (always better to start with tabula rasa)
rm(list = ls())

# Create needed export directory
dir.create(path = file.path("summary"), showWarnings = F)

## -------------------------------------------- ##
      # Data Retrieval & Housekeeping ####
## -------------------------------------------- ##
# Retrieve the relevant datasets
beta <- read.csv(file.path("data", "tidy_data", "beta-diversity-data.csv"))
fams <- read.csv(file.path("data", "tidy_data", "family-abun.csv"))

# Look at them to get familiar with structure
dplyr::glimpse(beta)
dplyr::glimpse(fams)

## -------------------------------------------- ##
          # Summarize Family/Phyla ####
## -------------------------------------------- ##

# Identify total number of families/phyla (across samples)
total_smry <- data.frame(Sample.ID = "All Samples",
                         family_ct = length(unique(fams$Family)),
                         phylum_ct = length(unique(fams$Phylum)))

# Identify number of families & phyla per sample
simp_df <- fams %>%
  dplyr::group_by(Sample.ID, Lifestage, Gut.Region, Stage.Gut, Sex) %>%
  dplyr::summarize(family_ct = length(unique(Family)),
                   phylum_ct = length(unique(Phylum))) %>%
  dplyr::ungroup() %>%
  # Bind on the total information too
  dplyr::bind_rows(y = total_smry)

# Check that out
dplyr::glimpse(simp_df)

## -------------------------------------------- ##
            # Summarize Reads ####
## -------------------------------------------- ##

# Rotate the beta diversity table to long format
beta_long <- beta %>%
  tidyr::pivot_longer(cols = -Sample.ID:-Sex,
                      names_to = "identity",
                      values_to = "count") %>%
  # Drop any 0 count rows that result
  dplyr::filter(count > 0)

# Check structure of new dataframe
dplyr::glimpse(beta_long)

# Identify the total unique read IDs
total_reads <- data.frame(Sample.ID = "All Samples",
                          read_ct = length(unique(beta_long$identity)))

# Identify unique reads per sample
simp_reads <- beta_long %>%
  dplyr::group_by(Sample.ID, Lifestage, Gut.Region, Stage.Gut, Sex) %>%
  dplyr::summarize(read_ct = length(unique(identity))) %>%
  dplyr::ungroup() %>%
  # Attach the total read count dataframe (across samples)
  dplyr::bind_rows(y = total_reads)

# Check structure
dplyr::glimpse(simp_reads)

## -------------------------------------------- ##
# Combine Summaries ####
## -------------------------------------------- ##

# Need to combine these into one object
total_summary <- simp_reads %>%
  # Left join on the read count
  dplyr::left_join(y = simp_df,
                   by = c("Sample.ID", "Lifestage", "Gut.Region", "Stage.Gut", "Sex")) %>%
  # Rename the read count column
  dplyr::rename(asv_ct = read_ct) %>%
  # Re-arrange so 'all samples' is on top
  dplyr::arrange(Sample.ID)

# Check structure
dplyr::glimpse(total_summary)

# Save this out
write.csv(x = total_summary, file.path("summary", "cotinis_summary_table.csv"),
          row.names = F, na = "")

# End ####
