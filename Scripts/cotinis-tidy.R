## ---------------------------------------------------------------------------------- ##
                  # Kucuk Cotinis Project - Visualization Code
## ---------------------------------------------------------------------------------- ##
# Code written by Nicholas J Lyon

# PURPOSE:
  ## Transform the data from .qza files into dataframes that R can understand
  ## Those "regular" dataframes will be used by subsequent codes

# Clear environment (always better to start with tabula rasa)
rm(list = ls())

# Set working directory
getwd() # should end in ".../Kucuk-CotinisCollab"
myWD <- getwd()

# Get the latest version of the package that will allow for the import of .qza files
  ## .qza is the format of the data from QIIME2
#remotes::install_github("jbisanz/qiime2R")
    ### This line will prompt you to update any needed packages
    ### Responding with "3" (no quotes) will skip the update step

# Necessary libraries
library(remotes); library(qiime2R); library(tidyverse)
library(stringr); library(vegan); library(fossil)

## ------------------------------------------- ##
                # Load Data ####
## ------------------------------------------- ##
# Read in the data
  ## Community dataset
comm.qza <- read_qza("./Data/cotinisdeblur/table.qza")
comm.data <- as.data.frame(comm.qza$data)

  ## FeatureID to taxonomy file
tax.qza <- read_qza("./Data/diversity6/taxonomy.qza")
tax.data <- as.data.frame(tax.qza$data)

  ## Metadata file connecting sample ID with sex/life stage/etc.
meta <- read.csv("./Data/cotinis-metadata.csv")

# Look at what we've got
str(comm.data)
str(tax.data)
str(meta)

# We can only start with a transposed dataframe
comm.data.v2 <- as.data.frame(t(comm.data))

## ------------------------------------------- ##
        # Integrate & Tidy Metadata ####
## ------------------------------------------- ##
# Process data
comm.data.v3 <- comm.data.v2 %>%
  ## Create a column from the rownames
  dplyr::mutate(X.SampleID = rownames(comm.data.v2)) %>%
  ## Bring in the metadata
  left_join(meta, by = "X.SampleID") %>%
  ## Get a differently-named version of Sample ID
  dplyr::mutate(Sample.ID = X.SampleID) %>%
  ## Remove unnecessary columns
  select(-X.SampleID:-Fraction., -Generalgutregion) %>%
  ## Rename the gut region column as two new columns
  dplyr::mutate(
    Gut.Region = Specificgutregion,
    Stage.Gut = Specificgutregion
      ) %>%
  ## Remove the specific gut region column too
  select(-Specificgutregion) %>%
  ## Move the columns we kept to the front
  relocate(Sample.ID, Lifestage, Gut.Region, Stage.Gut, Sex) %>%
  ## Filter out the unfed adults
  filter(Sample.ID != "Unfedfemale1hind" & Sample.ID != "Unfedfemale1mid" &
           Sample.ID != "UnfedFemale1Mid" & Sample.ID != "UnfedFemale2hind" &
           Sample.ID != "UnfedMale1hind" & Sample.ID != "UnfedMale1Mid") %>%
  ## Return dataframe
  as.data.frame()

# Look at what we brought in
names(comm.data.v3[1:10])

# Tidy those metadata columns
sort(unique(comm.data.v3$Sample.ID))
sort(unique(comm.data.v3$Lifestage))
sort(unique(comm.data.v3$Stage.Gut))
sort(unique(comm.data.v3$Sex))
  ## These look good

# The gut region needs simplification
sort(unique(comm.data.v3$Gut.Region))
comm.data.v3$Gut.Region <- gsub("Adult |Larval ", "", comm.data.v3$Gut.Region)
sort(unique(comm.data.v3$Gut.Region))

# Examine the starts and ends of the dataframe
names(comm.data.v3[c(1:10, (ncol(comm.data.v3)-10):ncol(comm.data.v3))])

# Save this out for beta diversity calculation in another script
write.csv(comm.data.v3, "./Data/beta-diveristy-data.csv", row.names = F)

## ------------------------------------------- ##
    # Calculate Alpha Diversity Indices ####
## ------------------------------------------- ##
# Calculate each index
dive.indexes <- comm.data.v3 %>%
  ## The following calculations should be done for every row
  rowwise() %>%
  ## Calculate the indexes
  mutate(
    Shannon = vegan::diversity(across(-Sample.ID:-Sex),
                               index = 'shannon'),
    Simpson = vegan::diversity(across(-Sample.ID:-Sex),
                               index = 'simpson'),
    Richness = vegan::specnumber(across(-Sample.ID:-Sex)),
    Pielous = Shannon / log(Simpson),
    ACE = fossil::ACE(across(-Sample.ID:-Sex),
                      taxa.row = T),
    Chao1 = fossil::chao1(across(-Sample.ID:-Sex),
                          taxa.row = T)
    ) %>%
  ## Remove the actual community (for ease of plotting)
  select(-`9a5171a5b50ffc0b6b48abc366e0076b`:-`05584c5da9aef5857819f6cedd30adcf`,
         -Richness) %>%
  as.data.frame()

head(dive.indexes)

# This dataframe is ready for use in stats/visualization!
write.csv(dive.indexes, "./Data/alpha-diveristy-data.csv", row.names = F)

# END ####

