## ---------------------------------------------------------------------------------- ##
                  # Kucuk Cotinis Project - Wrangling Code
## ---------------------------------------------------------------------------------- ##
# Code written by Nicholas J Lyon

# PURPOSE:
  ## Transform the data from .qza files into dataframes that R can understand
  ## Those "regular" dataframes will be used by subsequent codes

# Clear environment (always better to start with tabula rasa)
rm(list = ls())

# Load needed libraries
# install.packages("librarian")
librarian::shelf(tidyverse, vegan, fossil, ape, phyloseq, jbisanz/qiime2R)

# Create folder to export tidy data to
dir.create("Tidy Data", showWarnings = F)

## ------------------------------------------- ##
                # Load Data ####
## ------------------------------------------- ##
# Read in the data
  ## Community dataset
comm.qza <- qiime2R::read_qza(file = file.path("Data", "Raw Data", "table.qza"))
comm.data <- as.data.frame(comm.qza$data)

  ## FeatureID to taxonomy information
tax.qza <- qiime2R::read_qza(file = file.path("Data", "Raw Data", "taxonomy.qza"))
tax.data <- as.data.frame(tax.qza$data)

  ## Unrooted tree
tree.qza <- qiime2R::read_qza(file = file.path("Data", "Raw Data", "unrooted-tree.qza"))
tree.data <- tree.qza$data

  ## Metadata file connecting sample ID with sex/life stage/etc.
meta <- read.csv(file = file.path("Data", "Raw Data", "cotinis-metadata.csv"))

# Look at what we've got
str(comm.data)
str(tax.data)
str(tree.data)
str(meta)

## -------------------------------------------------------------- ##
               # Part 1: Taxonomy Metadata Tidying ####
## -------------------------------------------------------------- ##
## ------------------------------------------- ##
    # P1-1: Fix Concatenated Taxon ID ####
## ------------------------------------------- ##
# Split the Taxon by the semicolons (it is too dense as is)
tax.data.v2 <- tax.data %>%
  tidyr::separate(col = Taxon, sep = ";D_", remove = F,
           into = c("Domain", "Phylum", "Subphylum", "Order", "Family", "Genus", "Species"))

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

## ------------------------------------------- ##
  # P1-2: Fix Duplicate IDs (in same row) ####
## ------------------------------------------- ##
# Some columns are duplicated
  # (i.e., family, genus, and species have identical entries)
# Need to fix that
tax.data.v3 <- tax.data.v2 %>%
  # Rename the first version of the Taxon column while we're here
  dplyr::rename(Old.Taxon = Taxon) %>%
  # Replace each cell with NA if it's an exact match with its broader version
  dplyr::mutate(Species = ifelse(test = Species == Genus,
                                 yes = NA, no = Species),
                Genus = ifelse(test = Genus == Family,
                               yes = NA, no = Genus),
                Family = ifelse(test = Family == Order,
                                yes = NA, no = Family),
                Order = ifelse(test = Order == Subphylum,
                               yes = NA, no = Order),
                Subphylum = ifelse(test = Subphylum == Phylum,
                                   yes = NA, no = Subphylum),
                Phylum = ifelse(test = Phylum == Domain,
                                yes = NA, no = Phylum)) %>%
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
  dplyr::mutate(Species.Bin = ifelse(test = is.na(Species) == T,
                                     yes = 0,  no = 1),
                Genus.Bin = ifelse(test = is.na(Genus) == T,
                                   yes = 0, no = 1),
                Family.Bin = ifelse(test = is.na(Family) == T,
                                    yes = 0, no = 1),
                Order.Bin = ifelse(test = is.na(Order) == T,
                                   yes = 0, no = 1),
                Subphylum.Bin = ifelse(test = is.na(Subphylum) == T,
                                       yes = 0, no = 1),
                Phylum.Bin = ifelse(test = is.na(Phylum) == T,
                                    yes = 0, no = 1),
                Domain.Bin = ifelse(test = is.na(Domain) == T,
                                    yes = 0, no = 1),
                Barcode = paste(Domain.Bin, Phylum.Bin, Subphylum.Bin,
                                Order.Bin, Family.Bin, Genus.Bin,
                                Species.Bin, sep = '-')) %>%
  # return a df
  as.data.frame()

# Check out our new dataframe
str(tax.data.v4)

# Look at our barcode for specificity
sort(unique(tax.data.v4$Barcode))

# We can use that to assign correctly-formatted Taxa IDs!
tax.data.v5 <- tax.data.v4 %>%
  dplyr::mutate(BetterID = dplyr::case_when(
    Barcode == "1-0-0-0-0-0-0" ~ Domain,
    Barcode %in% c("1-1-0-0-0-0-0", "1-1-1-0-0-0-0") ~ paste(Domain, Phylum, sep = '_'),
    Barcode %in% c("1-1-0-1-0-0-0", "1-1-1-1-0-0-0") ~ paste(Domain, Phylum, Order,
                                                             sep = '_'),
    Barcode %in% c("1-1-0-1-1-0-0", "1-1-1-1-1-0-0") ~ paste(Domain, Phylum, Order,
                                                             Family, sep = '_'),
    Barcode %in% c("1-1-0-1-1-1-0", "1-1-1-1-1-1-0") ~ paste(Domain, Phylum, Order,
                                                             Family, Genus, sep = '_'),
    Barcode %in% c("1-1-0-1-1-1-1", "1-1-1-1-1-1-1") ~ paste(Domain, Phylum, Order,
                                                             Family, Genus, Species,
                                                             sep = '_'),
    # If the barcode *isn't* one of those, fill with NA
    TRUE ~ NA)) %>%
  # Get an integer counter for duplicate IDs
  dplyr::group_by(BetterID) %>%
  dplyr::mutate(Dup.Counter = seq_along(BetterID)) %>%
  dplyr::ungroup() %>%
  # Make a finalized Taxon column
  dplyr::mutate(Taxon = paste(BetterID, Dup.Counter, sep = '_')) %>%
  # Remove some now-unneeded columns
  dplyr::select(-Species.Bin:-Dup.Counter) %>%
  # Return a dataframe
  as.data.frame()

# Check out the structure of the data
str(tax.data.v5)

# And look at our new ID
sort(unique(tax.data.v5$Taxon))
plyr::count(is.na(tax.data.v5$Taxon))
  ## No NAs!

# Save this out for later reference
write.csv(x = tax.data.v5, row.names = F, na = '',
          file = file.path("Data", "Tidy Data", "taxonomy.csv"))

## -------------------------------------------------------------- ##
                # Part 2: Community Data Tidying ####
## -------------------------------------------------------------- ##
## ------------------------------------------- ##
    # P2-1: Exchange FeatureID for Taxon ####
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
     # P2-2: Integrate & Tidy Metadata ####
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
write.csv(comm.data.v3,
          "./Data/Tidy Data/beta-diversity-data.csv",
          row.names = F)

## ------------------------------------------- ##
     # P2-3: Calculate Alpha Diversity ####
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
write.csv(x = dive.indexes, row.names = F, na = '',
          file = file.path("Data", "Tidy Data", "alpha-diversity-data.csv"))

## -------------------------------------------------------------- ##
            # Part 3: Phylogenetic Tree Tidying ####
## -------------------------------------------------------------- ##
# Look at the tree file
str(tree.data)

# Make a duplicate
tree.data.v2 <- tree.data

## ------------------------------------------- ##
          # P3-1: Fix Tip Labels ####
## ------------------------------------------- ##
# Tip labels are named with FeatureID which is not informative
  ## Want to replace with our tidied "Taxon" column
tree.data.v2$tip.label.actual <- tax.data.v5$Taxon[match(tree.data$tip.label, tax.data.v5$Feature.ID)]

# Check it out
str(tree.data.v2)

# Remove the old tip label
tree.data.v2$tip.label <- NULL

# Re-name the new tip label so later functions don't get suprised and frightened
tree.data.v2$tip.label <- tree.data.v2$tip.label.actual
tree.data.v2$tip.label.actual <- NULL

# Double check what we're left with
str(tree.data.v2)

# Save it
ape::write.tree(phy = tree.data.v2, file = file.path("Data", "Tidy Data", "unrooted-tree.tre"))

## -------------------------------------------------------------- ##
              # Part 4: Taxon Abundance Tidying ####
## -------------------------------------------------------------- ##
## ------------------------------------------- ##
    # P4-1: Migrate Phylum and Family ####
## ------------------------------------------- ##
# Check structure of taxon metadata
str(tax.data.v5)

# And the community data
str(comm.data.v3)

# Create the necessary community dataset
comm.data.v4 <- comm.data.v3 %>%
  # Pivot to long format
  pivot_longer(
    cols = -Sample.ID:-Sex,
    names_to = "Taxon",
    values_to = "Abundance"
  ) %>%
  # Bring over all of the columns from the taxonomic data
  left_join(tax.data.v5, by = "Taxon") %>%
  # Keep only desired columns (non-specified columns are dropped)
  dplyr::select(Sample.ID, Lifestage, Gut.Region, Stage.Gut, Sex,
                Domain, Phylum, Family, Abundance) %>%
  # Return df
  as.data.frame()

# Check what we've created
str(comm.data.v4)

## ------------------------------------------- ##
   # P4-2: Check Phylum & Family Columns ####
## ------------------------------------------- ##
# Check domain (shouldn't have issues but better safe than sorry)
sort(unique(comm.data.v4$Domain))

# Check phylum
sort(unique(comm.data.v4$Phylum))
  ## "RsaHF231" is a candidate phylum so we'll allow it (for now)

# Check family
sort(unique(comm.data.v4$Family))
  ## Some worrisome so let's get more detail on those

# Make a vector of potentially incorrect family IDs
check.fam <- c("0319-6G20", "67-14", "A0839", "A4b", "AKYG1722", "bacteriap25",
             "BIrii41", "Brachyspirales Incertae Sedis", "Clade I",
             "Clostridiaceae 1", "Clostridiales vadinBB60 group",
             "Coriobacteriales Incertae Sedis", "D05-2", "Family XIII",
             "JG30-KF-CM45", "KF-JG30-B3", "metagenome",
             "Rhizobiales Incertae Sedis", "SC-I-84",
             "Solibacteraceae (Subgroup 3)", "uncultured", "uncultured bacterium",
             "uncultured Chloroflexi bacterium", "uncultured delta proteobacterium",
             "uncultured Firmicutes bacterium", "uncultured Mollicutes bacterium",
             "uncultured rumen bacterium", "uncultured soil bacterium",
             "Unknown Family", "wastewater metagenome", "WD2101 soil group")

# Examine the taxonomic data we have on these dubious "families"
dubious.fams <- tax.data.v5 %>%
  dplyr::filter(Family %in% check.fam) %>%
  as.data.frame()

# Look at this
# View(dubious.fams)

# Because it's abundance, we can leave this alone (though may need to return)

## ------------------------------------------- ##
      # P4-3: Family-Level Abundance ####
## ------------------------------------------- ##
# Want family-level abundance
fam.abun <- comm.data.v4 %>%
  # Group by needed columns
  group_by(Sample.ID, Lifestage, Gut.Region, Stage.Gut, Sex,
           Domain, Phylum, Family) %>%
  # Sum abundance within these groups
  summarise(Abundance = sum(Abundance, na.rm = T)) %>%
  # Ungroup
  ungroup() %>%
  # Remove any NAs
  filter(!is.na(Family)) %>%
  # Remove any 0 abundances
  filter(Abundance != 0) %>%
  # Also want to calculate relative abundance
    ## "relative" means within sample calculations are needed
  group_by(Sample.ID) %>%
  dplyr::mutate(totalAbun = sum(Abundance, na.rm = T)) %>%
  ungroup() %>%
  # Calculate relative abundance
  dplyr::mutate(relativeAbun = (Abundance / totalAbun) * 100) %>%
  # Return df
  as.data.frame()

# Should have lost *a lot* of rows
nrow(comm.data.v4) - nrow(fam.abun)
dim(comm.data.v4); dim(fam.abun)

# Save this file for later use
write.csv(x = fam.abun, row.names = F, na = '',
          file = file.path("Data", "Tidy Data", "family-abun.csv"))

## ------------------------------------------- ##
      # P4-4: Phylum-Level Abundance ####
## ------------------------------------------- ##
# Want phylum-level abundance
phyla.abun <- comm.data.v4 %>%
  # Remove family column
    ## This would be done implicitly by group_by()
  dplyr::select(-Family) %>%
  # Group by needed columns
  group_by(Sample.ID, Lifestage, Gut.Region, Stage.Gut, Sex, Domain, Phylum) %>%
  # Sum abundance within these groups
  summarise(Abundance = sum(Abundance, na.rm = T)) %>%
  # Ungroup
  ungroup() %>%
  # Remove any NAs
  filter(!is.na(Phylum)) %>%
  # Remove any 0 abundances
  filter(Abundance != 0) %>%
  # Also want to calculate relative abundance
    ## "relative" means within sample calculations are needed
  group_by(Sample.ID) %>%
  # Calculate total abundance for each sample
  dplyr::mutate(totalAbun = sum(Abundance)) %>%
  # Ungroup
  ungroup() %>%
  # Calculate relative abundance
  dplyr::mutate(relativeAbun = (Abundance / totalAbun) * 100) %>%
  # Return df
  as.data.frame()

# Should have lost *a lot* of rows
nrow(comm.data.v4) - nrow(phyla.abun)
dim(comm.data.v4); dim(phyla.abun)

# Save this file for later use
write.csv(x = phyla.abun, row.names = F, na = '',
          file = file.path("Data", "Tidy Data", "phylum-abun.csv"))

## ------------------------------------------- ##
       # P4-5: Genus-Level Abundance ####
## ------------------------------------------- ##
# Examine genera
sort(unique(tax.data.v5$Genus))

# Create the necessary community dataset
gen.abun <- comm.data.v3 %>%
  # Pivot to long format
  pivot_longer(
    cols = -Sample.ID:-Sex,
    names_to = "Taxon",
    values_to = "Abundance"
  ) %>%
  # Bring over all of the columns from the taxonomic data
  left_join(tax.data.v5, by = "Taxon") %>%
  # Keep only desired columns (non-specified columns are dropped)
  dplyr::select(Sample.ID, Lifestage, Gut.Region, Stage.Gut, Sex,
                Genus, Abundance) %>%
  # Group by needed columns
  group_by(Sample.ID, Lifestage, Gut.Region, Stage.Gut, Sex, Genus) %>%
  # Sum abundance within these groups
  summarise(Abundance = sum(Abundance, na.rm = T)) %>%
  # Ungroup
  ungroup() %>%
  # Remove rows without a genus-level ID
  filter(!is.na(Genus)) %>%
  # Return df
  as.data.frame()

# Check what we've created
str(gen.abun)

# Save it
write.csv(x = gen.abun, row.names = F, na = '',
          file = file.path("Data", "Tidy Data", "genus-abun.csv"))

## -------------------------------------------------------------- ##
             # Part 5: Calculate Unifrac Distances ####
## -------------------------------------------------------------- ##
# Re-check what you may need for unifrac distance calculation
  ## Abundance table
str(comm.data.v3)
class(comm.data.v3)

  ## Phylogenetic tree
str(tree.data.v2)
class(tree.data.v2)

## ------------------------------------------- ##
         # P5-1: Weighted Unifrac ####
## ------------------------------------------- ##
# https://rdrr.io/github/MadsAlbertsen/ampvis2/man/unifrac.html
# unifrac::unifrac(tree.data.v2)



# END ####

