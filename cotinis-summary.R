## -------------------------------------------------------------- ##
# Kucuk Cotinis Project - Summary Code
## -------------------------------------------------------------- ##
# Code written by Nicholas J Lyon

# PURPOSE:
## Get some useful summaries of the data

# Necessary libraries
# install.packages("librarian")
librarian::shelf(tidyverse, supportR, jbisanz/qiime2R)

# Clear environment (always better to start with tabula rasa)
rm(list = ls())

# Create needed export directory
dir.create(path = file.path("summary"), showWarnings = F)

## -------------------------------------------- ##
# Summarize Reads ####
## -------------------------------------------- ##

# Read in the denoised reads QZA file
# reads_qza <- qiime2R::read_qza(file = file.path("data", "raw_data", "rep-seqs.qza"))
reads_qza <- qiime2R::read_qza(file = file.path("data", "raw_data", "reads-joined.qza"))

# Strip out the dataframe portion
reads <- as.data.frame(reads_qza$data)

# Check structure
dplyr::glimpse(reads)

# Identify number of reads globally
total_reads <- data.frame(Sample.ID = "All Samples",
                         read_ct = length(unique(reads$x)))

# Check that out
total_reads

## -------------------------------------------- ##
          # Summarize Family/Phyla ####
## -------------------------------------------- ##
# Read in family/phyla information
fams <- read.csv(file.path("data", "tidy_data", "family-abun.csv"))

# Check structure
dplyr::glimpse(fams)

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
            # Summarize ASVs ####
## -------------------------------------------- ##
# Retrieve the ASVs info (basically 'species' for our purposes)
beta <- read.csv(file.path("data", "tidy_data", "beta-diversity-data.csv"))

# Look at structure
dplyr::glimpse(beta)

# Rotate the beta diversity table to long format
beta_long <- beta %>%
  tidyr::pivot_longer(cols = -Sample.ID:-Sex,
                      names_to = "identity",
                      values_to = "count") %>%
  # Drop any 0 count rows that result
  dplyr::filter(count > 0)

# Check structure of new dataframe
dplyr::glimpse(beta_long)

# Identify the total unique ASV IDs (essentially 'species')
total_asvs <- data.frame(Sample.ID = "All Samples",
                          asv_ct = length(unique(beta_long$identity)))

# Identify unique ASVs per sample
simp_asvs <- beta_long %>%
  dplyr::group_by(Sample.ID, Lifestage, Gut.Region, Stage.Gut, Sex) %>%
  dplyr::summarize(asv_ct = length(unique(identity))) %>%
  dplyr::ungroup() %>%
  # Attach the total ASV count dataframe (across samples)
  dplyr::bind_rows(y = total_asvs)

# Check structure
dplyr::glimpse(simp_asvs)

## -------------------------------------------- ##
# Combine Summaries ####
## -------------------------------------------- ##

# Need to combine these into one object
total_summary <- simp_asvs %>%
  # Left join on the ASV count
  dplyr::left_join(y = simp_df,
                   by = c("Sample.ID", "Lifestage", "Gut.Region", "Stage.Gut", "Sex")) %>%
  # Rename the ASV count column
  dplyr::rename(species_ct = asv_ct) %>%
  # Re-arrange so 'all samples' is on top
  dplyr::arrange(Sample.ID)

# Check structure
dplyr::glimpse(total_summary)

# Save this out
write.csv(x = total_summary, file.path("summary", "cotinis_summary_table.csv"),
          row.names = F, na = "")

# End ####
