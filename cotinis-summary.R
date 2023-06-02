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
# Read in the joined, joined + filtered, and denoised reads QZA files
reads_join_qza <- qiime2R::read_qza(file = file.path("data", "raw_data", "reads-joined.qza"))
reads_jofi_qza <- qiime2R::read_qza(file = file.path("data", "raw_data", "reads-joined-filtered.qza"))
reads_deno_qza <- qiime2R::read_qza(file = file.path("data", "raw_data", "rep-seqs.qza"))

# Strip out the dataframe portion of each
reads_join <- as.data.frame(reads_join_qza$data) %>%
  # Rename a column
  dplyr::rename(join_size = size)
reads_jofi <- as.data.frame(reads_jofi_qza$data) %>%
  dplyr::rename(jofi_size = size)
reads_deno <- as.data.frame(reads_deno_qza$data)

# Check structure
dplyr::glimpse(reads_join)
dplyr::glimpse(reads_jofi)
dplyr::glimpse(reads_deno)

# Combine the join and join + filter dataframes
reads_by_sample <- reads_join %>%
  dplyr::left_join(y = reads_jofi, by = "files") %>%
  # Make both columns truly numeric
  dplyr::mutate(join_read_sum = as.numeric(join_size),
                joinfilter_read_sum = as.numeric(jofi_size)) %>%
  # Drop data from unfed individuals
  dplyr::filter(stringr::str_sub(string = files, start = 1, end = 5) != "Unfed") %>%
  # Pare down the files column into just sample ID
  dplyr::mutate(Sample.ID = gsub(pattern = "_[[:digit:]]{1,2}_L001_R1_001.fastq.gz",
                                 replacement = "", x = files),
                .before = dplyr::everything()) %>%
  # Drop unwanted non-sample information too
  dplyr::filter(!Sample.ID %in% c("MANIFEST", "metadata.yml")) %>%
  # Pare down to only desired columns
  dplyr::select(Sample.ID, join_read_sum, joinfilter_read_sum)

# Re-check structure
dplyr::glimpse(reads_by_sample)

# Identify number of reads globally
total_reads <- data.frame(Sample.ID = "All Samples",
                          join_read_sum = sum(reads_by_sample$join_read_sum, na.rm = T),
                          joinfilter_read_sum = sum(reads_by_sample$joinfilter_read_sum, na.rm = T),
                          denoised_read_ct = length(unique(reads_deno$x)))

# Check that out
total_reads

# Combine these into one data object
reads_actual <- total_reads %>%
  dplyr::bind_rows(y = reads_by_sample)

# Check structure
dplyr::glimpse(reads_actual)

## -------------------------------------------- ##
          # Summarize Family/Phyla ####
## -------------------------------------------- ##
# Read in family/phyla information
fams <- read.csv(file.path("data", "tidy_data", "family-abun.csv"))

# Check structure
dplyr::glimpse(fams)

# Identify number of families & phyla per sample
fams_by_sample <- fams %>%
  dplyr::group_by(Sample.ID, Lifestage, Gut.Region, Stage.Gut, Sex) %>%
  dplyr::summarize(family_ct = length(unique(Family)),
                   phylum_ct = length(unique(Phylum))) %>%
  dplyr::ungroup()

# Check that out
dplyr::glimpse(fams_by_sample)

# Identify total number of families/phyla (across samples)
total_fams <- data.frame(Sample.ID = "All Samples",
                         family_ct = length(unique(fams$Family)),
                         phylum_ct = length(unique(fams$Phylum)))

# Check it out
total_fams

# Combine these
fams_actual <- total_fams %>%
  dplyr::bind_rows(y = fams_by_sample)

# Re-check structure
dplyr::glimpse(fams_actual)

## -------------------------------------------- ##
            # Summarize ASVs ####
## -------------------------------------------- ##
# Retrieve the ASVs info (basically 'species' for our purposes)
spp <- read.csv(file.path("data", "tidy_data", "beta-diversity-data.csv"))

# Look at structure
dplyr::glimpse(spp)

# Rotate the beta diversity table to long format
spp_long <- spp %>%
  tidyr::pivot_longer(cols = -Sample.ID:-Sex,
                      names_to = "identity",
                      values_to = "count") %>%
  # Drop any 0 count rows that result
  dplyr::filter(count > 0)

# Check structure of new dataframe
dplyr::glimpse(spp_long)

# Identify unique ASVs per sample
asvs_by_sample <- spp_long %>%
  dplyr::group_by(Sample.ID, Lifestage, Gut.Region, Stage.Gut, Sex) %>%
  dplyr::summarize(asv_ct = length(unique(identity))) %>%
  dplyr::ungroup()

# Glimpse
dplyr::glimpse(asvs_by_sample)

# Identify the total unique ASV IDs (essentially 'species')
total_asvs <- data.frame(Sample.ID = "All Samples",
                          asv_ct = length(unique(spp_long$identity)))

# Look at it
total_asvs

# Combine the two
asvs_actual <- total_asvs %>%
  dplyr::bind_rows(asvs_by_sample) %>%
  # Rename ASV column
  dplyr::rename(species_ct = asv_ct)

# Check structure
dplyr::glimpse(asvs_actual)

## -------------------------------------------- ##
          # Combine Summaries ####
## -------------------------------------------- ##
# Need to combine these into one object
total_summary <- reads_actual %>%
  # Left join on family information
  dplyr::left_join(y = fams_actual, by = c("Sample.ID")) %>%
  # Left join on species (i.e., ASVs)
  dplyr::left_join(y = asvs_actual, c("Sample.ID", "Lifestage", "Gut.Region", "Stage.Gut", "Sex")) %>%
  # Move sample info to the left
  dplyr::relocate(Lifestage:Sex, .after = Sample.ID) %>%
  # Drop Lmid19 (included in reads but not other data)
  dplyr::filter(Sample.ID != "Lmid19")

# Check structure
dplyr::glimpse(total_summary)

# Save this out
write.csv(x = total_summary, file.path("summary", "cotinis_summary_table.csv"),
          row.names = F, na = "")

# End ####
