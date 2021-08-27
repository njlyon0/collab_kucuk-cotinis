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

## ------------------------------------------- ##
     # Get a More Meaningful Taxon ID ####
## ------------------------------------------- ##
# Split the Taxon by the semicolons (it is too dense as is)
tax.data.v2 <- tax.data %>%
  separate(col = Taxon,
           into = c("Domain", "Phylum", "Subphylum",
                    "Order", "Family", "Genus",
                    "Species"),
           sep = ";D_", remove = F)

# Check the structure now
str(tax.data.v2)

# Examine each of the new columns
sort(unique(tax.data.v2$Domain)) # starts with "D_0__"
sort(unique(tax.data.v2$Phylum)) # starts with "1__"
sort(unique(tax.data.v2$Subphylum)) # starts with "2__"
sort(unique(tax.data.v2$Order)) # starts with "3__"
sort(unique(tax.data.v2$Family)) # starts with "4__"
sort(unique(tax.data.v2$Genus)) # starts with "5__"
sort(unique(tax.data.v2$Species)) # starts with "6__"

# Remove the preceeding numbers for each column
  ## starts with "D_0__"
tax.data.v2$Domain <- gsub("D_0__", "", tax.data.v2$Domain)
sort(unique(tax.data.v2$Domain))

  ## starts with "1__"
tax.data.v2$Phylum <- gsub("1__", "", tax.data.v2$Phylum)
sort(unique(tax.data.v2$Phylum)) 

  ## starts with "2__"
tax.data.v2$Subphylum <- gsub("2__", "", tax.data.v2$Subphylum)
sort(unique(tax.data.v2$Subphylum))

  ## starts with "3__"
tax.data.v2$Order <- gsub("3__", "", tax.data.v2$Order)
sort(unique(tax.data.v2$Order)) 

  ## starts with "4__"
tax.data.v2$Family <- gsub("4__", "", tax.data.v2$Family)
sort(unique(tax.data.v2$Family)) 

  ## starts with "5__"
tax.data.v2$Genus <- gsub("5__", "", tax.data.v2$Genus)
sort(unique(tax.data.v2$Genus)) 

  ## starts with "6__"
tax.data.v2$Species <- gsub("6__", "", tax.data.v2$Species)
sort(unique(tax.data.v2$Species)) 

# How many NAs in the lower classifications?
plyr::count(is.na(tax.data.v2$Species))
plyr::count(is.na(tax.data.v2$Genus))
plyr::count(is.na(tax.data.v2$Family))
plyr::count(is.na(tax.data.v2$Order))
plyr::count(is.na(tax.data.v2$Subphylum))
plyr::count(is.na(tax.data.v2$Phylum))
plyr::count(is.na(tax.data.v2$Domain))

# What columns do we have access to?
names(tax.data.v2)

# Some columns are duplicated
  # (i.e., family, genus, and species have identical entries)
# Need to fix that
tax.data.v3 <- tax.data.v2 %>%
  # Rename the first version of the Taxon column while we're here
  dplyr::rename(Old.Taxon = Taxon) %>%
  # Replace each cell with NA if it's an exact match with its broader version
  dplyr::mutate(
    Species = ifelse(test = Species == Genus,
                     yes = NA,
                     no = Species),
    Genus = ifelse(test = Genus == Family,
                     yes = NA,
                     no = Genus),
    Family = ifelse(test = Family == Order,
                     yes = NA,
                     no = Family),
    Order = ifelse(test = Order == Subphylum,
                    yes = NA,
                    no = Order),
    Subphylum = ifelse(test = Subphylum == Phylum,
                    yes = NA,
                    no = Subphylum),
    Phylum = ifelse(test = Phylum == Domain,
                       yes = NA,
                       no = Phylum)
  ) %>%
  
  # Return a dataframe
  as.data.frame()

# Re-check NA frequency
plyr::count(is.na(tax.data.v2$Species))
plyr::count(is.na(tax.data.v3$Species)) # 460 more

plyr::count(is.na(tax.data.v2$Genus))
plyr::count(is.na(tax.data.v3$Genus)) # 197 more

plyr::count(is.na(tax.data.v2$Family))
plyr::count(is.na(tax.data.v3$Family)) # 101 more

plyr::count(is.na(tax.data.v2$Order))
plyr::count(is.na(tax.data.v3$Order)) # 1 more

plyr::count(is.na(tax.data.v2$Subphylum))
plyr::count(is.na(tax.data.v3$Subphylum)) # 274 more

plyr::count(is.na(tax.data.v2$Phylum))
plyr::count(is.na(tax.data.v3$Phylum)) # same

# We're going to make a diagnostic column
  ## to figure out how we can make the new taxon column
tax.data.v4 <- tax.data.v3 %>%
  dplyr::mutate(
    Species.Bin = ifelse(test = is.na(Species) == T,
                     yes = 0,
                     no = 1),
    Genus.Bin = ifelse(test = is.na(Genus) == T,
                       yes = 0,
                       no = 1),
    Family.Bin = ifelse(test = is.na(Family) == T,
                    yes = 0,
                    no = 1),
    Order.Bin = ifelse(test = is.na(Order) == T,
                       yes = 0,
                       no = 1),
    Subphylum.Bin = ifelse(test = is.na(Subphylum) == T,
                           yes = 0,
                           no = 1),
    Phylum.Bin = ifelse(test = is.na(Phylum) == T,
                        yes = 0,
                        no = 1),
    Domain.Bin = ifelse(test = is.na(Domain) == T,
                        yes = 0,
                        no = 1),
    Barcode = paste(Domain.Bin, Phylum.Bin, Subphylum.Bin,
                           Order.Bin, Family.Bin, Genus.Bin,
                           Species.Bin, sep = '-')
  ) %>%
  # return a df
  as.data.frame()

# Check out our new dataframe
str(tax.data.v4)

# Look at our barcode for specificity
sort(unique(tax.data.v4$Barcode))

# We can use that to assign correctly-formatted Taxa IDs!
tax.data.v5 <- tax.data.v4 %>%
  # This will be a nested ifelse but it'll be okay
  dplyr::mutate(
    BetterID = ifelse(Barcode == "1-0-0-0-0-0-0",
                   yes = Domain,
                   no = ifelse(
                    Barcode == "1-1-0-0-0-0-0" |
                    Barcode == "1-1-1-0-0-0-0",
                    yes = paste(Domain, Phylum,
                                sep = '_'),
                    no = ifelse(
                      Barcode == "1-1-0-1-0-0-0" | 
                      Barcode == "1-1-1-1-0-0-0",
                      yes = paste(Domain, Phylum,
                                  Order, sep = '_'),
                      no = ifelse(
                        Barcode == "1-1-0-1-1-0-0" |
                        Barcode == "1-1-1-1-1-0-0",
                        yes = paste(Domain, Phylum,
                                    Order, Family, sep = '_'),
                        no = ifelse(
                          Barcode == "1-1-0-1-1-1-0" |
                          Barcode == "1-1-1-1-1-1-0",
                          yes = paste(Domain, Phylum,
                                      Order, Family, Genus,
                                      sep = '_'),
                          no = ifelse(
                            Barcode == "1-1-0-1-1-1-1" |
                              Barcode == "1-1-1-1-1-1-1",
                            yes = paste(Domain, Phylum,
                                        Order, Family,
                                        Genus, Species, sep = '_'),
                            no = NA))))))
  ) %>%
  # Group by this ID
  group_by(BetterID) %>%
  # Get an integer counter for duplicate IDs
  dplyr::mutate(
        Dup.Counter = seq_along(BetterID)
  ) %>%
  # Ungroup
  ungroup(BetterID) %>%
  # Make a finalized Taxon column
  dplyr::mutate(
    Taxon = paste(BetterID, Dup.Counter, sep = '_')
                ) %>%
  # Remove some now-unneeded columns
  select(-Species.Bin:-Dup.Counter) %>%
  # Return a dataframe
  as.data.frame()

# Check out the structure of the data
str(tax.data.v5)

# And look at our new ID
sort(unique(tax.data.v5$Taxon))
plyr::count(is.na(tax.data.v5$Taxon))
  ## No NAs!

# Save this out for later reference
write.csv(tax.data.v5,
          "./Data/taxonomy.csv",
          row.names = F)

## ------------------------------------------- ##
   # Exchange FeatureID for new Taxon Col ####
## ------------------------------------------- ##
# Get featureID out of rownames and into a column
comm.data$FeatureID <- rownames(comm.data)
str(comm.data)

# Use that to bring in the TaxonID
comm.data$Taxon <- tax.data.v5$Taxon[match(comm.data$FeatureID, tax.data.v5$Feature.ID)]
str(comm.data)

# Make taxon the new rowname
rownames(comm.data) <- comm.data$Taxon

# Remove both the FeatureID and Taxon columns
comm.data.v1 <- comm.data %>%
  select(-FeatureID, -Taxon) %>%
  as.data.frame()

# We can only start with a transposed dataframe
comm.data.v2 <- as.data.frame(t(comm.data.v1))

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
write.csv(comm.data.v3, "./Data/beta-diversity-data.csv", row.names = F)

## ------------------------------------------- ##
    # Calculate Alpha Diversity Indices ####
## ------------------------------------------- ##
# Calculate each index
dive.indexes <- comm.data.v3 %>%
  ## The following calculations should be done for every row
  rowwise() %>%
  ## Calculate the indexes
  dplyr::mutate(
    Shannon = vegan::diversity(across(-Sample.ID:-Sex),
                               index = 'shannon'),
    Simpson = vegan::diversity(across(-Sample.ID:-Sex),
                               index = 'simpson'),
    Richness = vegan::specnumber(across(-Sample.ID:-Sex)),
    Pielous = (Shannon / log(Richness)),
    ACE = fossil::ACE(across(-Sample.ID:-Sex),
                      taxa.row = T),
    Chao1 = fossil::chao1(across(-Sample.ID:-Sex),
                          taxa.row = T)
    ) %>%
  ## Remove the actual community (for ease of plotting)
  select(Sample.ID:Sex, Shannon:Chao1) %>%
  as.data.frame()

head(dive.indexes)

# This dataframe is ready for use in stats/visualization!
write.csv(dive.indexes, "./Data/alpha-diversity-data.csv", row.names = F)

# END ####

