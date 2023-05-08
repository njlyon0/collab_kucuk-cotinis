## ---------------------------------------------------------------------------------- ##
                  # Kucuk Cotinis Project - Wrangling Code
## ---------------------------------------------------------------------------------- ##
# Code written by Nicholas J Lyon

# PURPOSE:
  ## Transform the data from .qza files into dataframes that R can understand
  ## Those "regular" dataframes will be used by subsequent codes

# Load needed libraries
# install.packages("librarian")
librarian::shelf(tidyverse, vegan, fossil, ape, phyloseq, jbisanz/qiime2R)

# Clear environment
rm(list = ls())

# Create folder to export tidy_data to
dir.create(path = file.path("data", "tidy_data"), showWarnings = F)

## ------------------------------------------- ##
                # Load Data ####
## ------------------------------------------- ##
# Read in the data
## Community dataset
comm_qza <- qiime2R::read_qza(file = file.path("data", "raw_data", "table.qza"))
comm_data <- as.data.frame(comm_qza$data)

## FeatureID to taxonomy information
tax_qza <- qiime2R::read_qza(file = file.path("data", "raw_data", "taxonomy.qza"))
tax_data <- as.data.frame(tax_qza$data)

## Unrooted tree
tree_qza <- qiime2R::read_qza(file = file.path("data", "raw_data", "unrooted-tree.qza"))
tree_data <- tree_qza$data

## Metadata file connecting sample ID with sex/life stage/etc.
meta <- read.csv(file = file.path("data", "raw_data", "cotinis-metadata.csv"))

# Look at what we've got
dplyr::glimpse(comm_data)
dplyr::glimpse(tax_data)
str(tree_data)
dplyr::glimpse(meta)

## -------------------------------------------------------------- ##
               # Part 1: Taxonomy Metadata Tidying ####
## -------------------------------------------------------------- ##
## ------------------------------------------- ##
    # P1-1: Fix Concatenated Taxon ID ####
## ------------------------------------------- ##
# Split the Taxon by the semicolons (it is too dense as is)
tax_v2 <- tax_data %>%
  tidyr::separate_wider_delim(cols = Taxon, delim = ";D_",
                              cols_remove = F, too_few = "align_start",
                              names = c("Domain", "Phylum", "Subphylum", "Order",
                                        "Family", "Genus", "Species"))

# Check the structure now
dplyr::glimpse(tax_v2)

# Examine each of the new columns
sort(unique(tax_v2$Domain)) # starts with "D_0__"
sort(unique(tax_v2$Phylum)) # starts with "1__"
sort(unique(tax_v2$Subphylum)) # starts with "2__"
sort(unique(tax_v2$Order)) # starts with "3__"
sort(unique(tax_v2$Family)) # starts with "4__"
sort(unique(tax_v2$Genus)) # starts with "5__"
sort(unique(tax_v2$Species)) # starts with "6__"

# Remove the preceeding numbers for each column
tax_v3 <- tax_v2 %>%
  dplyr::mutate(Domain = gsub(pattern = "D_0__", replacement = "", x = Domain),
                Phylum = gsub(pattern = "1__", replacement = "", x = Phylum),
                Subphylum = gsub(pattern = "2__", replacement = "", x = Subphylum),
                Order = gsub(pattern = "3__", replacement = "", x = Order),
                Family = gsub(pattern = "4__", replacement = "", x = Family),
                Genus = gsub(pattern = "5__", replacement = "", x = Genus),
                Species = gsub(pattern = "6__", replacement = "", x = Species))

# Check that worked
sort(unique(tax_v3$Domain))
sort(unique(tax_v3$Phylum))
sort(unique(tax_v3$Subphylum))
sort(unique(tax_v3$Order))
sort(unique(tax_v3$Family))
sort(unique(tax_v3$Genus))
sort(unique(tax_v3$Species))

# How many NAs in the lower classifications?
(na_freq <- tax_v3 %>%
  # Pivot longer
  tidyr::pivot_longer(cols = Domain:Species) %>%
  # Filter to only NAs
  dplyr::filter(is.na(value)) %>%
  # Count number of NAs per classification
  dplyr::group_by(name) %>%
  dplyr::summarize(raw_na_ct = dplyr::n()))

# What columns do we have access to?
dplyr::glimpse(tax_v3)

## ------------------------------------------- ##
  # P1-2: Fix Duplicate IDs (in same row) ####
## ------------------------------------------- ##
# Some columns are duplicated
  # (i.e., family, genus, and species have identical entries)
# Need to fix that
tax_v4 <- tax_v3 %>%
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
                                yes = NA, no = Phylum))

# Re-check NA frequency
tax_v4 %>%
  # Pivot longer
  tidyr::pivot_longer(cols = Domain:Species) %>%
  # Filter to only NAs
  dplyr::filter(is.na(value)) %>%
  # Count number of NAs per classification
  dplyr::group_by(name) %>%
  dplyr::summarize(actual_na_ct = dplyr::n()) %>%
  # Join on previous frequency assessment
  dplyr::full_join(y = na_freq, by = "name")

# More holistic structure assessment
dplyr::glimpse(tax_v4)

## ------------------------------------------- ##
        # P1-3: Refine Taxon Column ####
## ------------------------------------------- ##

# Get a cleaned up complete feature taxonomy
better_taxon <- tax_v4 %>%
  # Pivot long
  tidyr::pivot_longer(cols = Domain:Species,
                      names_to = "group",
                      values_to = "group_id") %>%
  # Drop NAs
  dplyr::filter(!is.na(group_id)) %>%
  # Group by feature ID and collapse re-formatted classifications
  dplyr::group_by(Feature.ID, Old.Taxon, Confidence) %>%
  dplyr::summarize(improved_taxon = paste(group_id, collapse = "_")) %>%
  dplyr::ungroup() %>%
  # Group by this revized ID Count number of recurring IDs
  dplyr::group_by(improved_taxon) %>%
  dplyr::mutate(rep_ct = seq_along(improved_taxon)) %>%
  dplyr::ungroup() %>%
  # Attach the replicate counter to the simplified taxon names
  dplyr::mutate(Taxon = ifelse(test = rep_ct != 1,
                               yes = paste0(improved_taxon, "_", rep_ct),
                               no = improved_taxon)) %>%
  # Remove unwanted columns
  dplyr::select(-improved_taxon, -rep_ct)

# Check out contents / NA frequency
sort(unique(better_taxon$Taxon))
better_taxon %>%
  dplyr::filter(is.na(Taxon)) %>%
  dplyr::summarize(na_ct = dplyr::n())

# Attach this to the actual dataframe
tax_v5 <- tax_v4 %>%
  dplyr::left_join(y = better_taxon, by = c("Feature.ID", "Old.Taxon", "Confidence")) %>%
  # Reorder columns a little as well
  dplyr::relocate(Old.Taxon, Taxon, .after = Feature.ID)

# Check out what we're left with
dplyr::glimpse(tax_v5)

# Save this out for later reference
write.csv(x = tax_v5, row.names = F, na = '',
          file = file.path("data", "tidy_data", "taxonomy.csv"))

## -------------------------------------------------------------- ##
                # Part 2: Community Data Tidying ####
## -------------------------------------------------------------- ##
## ------------------------------------------- ##
    # P2-1: Exchange FeatureID for Taxon ####
## ------------------------------------------- ##
# Do preliminary wrangling
comm_v2 <- comm_data %>%
  # Get featureID out of rownames and into a column
  dplyr::mutate(Feature.ID = rownames(.),
                .before = dplyr::everything()) %>%
  # Attach the refined taxonomy column we made above
  dplyr::mutate(Taxon = tax_v5$Taxon[match(.$Feature.ID, tax_v5$Feature.ID)],
                .after = Feature.ID)

# Re-check structure
dplyr::glimpse(comm_v2)

# Make the revised Taxon column into the new rownames
rownames(comm_v2) <- comm_v2$Taxon

# Drop the now-unwanted columns
comm_v3 <- comm_v2 %>%
  dplyr::select(-Feature.ID, -Taxon)

# Transpose dataframe to get features as columns and samples as rows
comm_v4 <- as.data.frame(t(comm_v3)) %>%
  # Add the rownames back in as a column
  dplyr::mutate(X.SampleID = rownames(.),
                .before = dplyr::everything())

# Check out what this leaves us with
dplyr::glimpse(comm_v4[1:5])

## ------------------------------------------- ##
     # P2-2: Integrate & Tidy Metadata ####
## ------------------------------------------- ##
# Integrate metadata
comm_v5 <- comm_v4 %>%
  # Attach all metadata
  dplyr::left_join(y = meta, by = "X.SampleID") %>%
  # Rename sample ID & specific gut region
  dplyr::rename(Sample.ID = X.SampleID,
                Stage.Gut = Specificgutregion) %>%
  # Get a 'gut region only' column
  dplyr::mutate(Gut.Region = gsub(pattern = "Adult |Larval ", replacement = "",
                                  x = Stage.Gut)) %>%
  # Remove unwanted columns
  dplyr::select(-BarcodeSequence, -LinkerPrimer, -Nucleotide,
                -Fraction., -Generalgutregion) %>%
  # Move metadata columns to left
  dplyr::relocate(Lifestage, Gut.Region, Stage.Gut, Sex,
                  .after = Sample.ID) %>%
  # Drop unwanted samples (i.e., "unfed" samples)
  dplyr::filter(!Sample.ID %in% c("Unfedfemale1hind", "Unfedfemale1mid", "UnfedFemale1Mid",
                                  "UnfedFemale2hind", "UnfedMale1hind", "UnfedMale1Mid"))

# Check structure of first bit of dataframe
dplyr::glimpse(comm_v5[1:5])

# What samples remain?
sort(unique(comm_v5$Sample.ID))

# Explore other metadata columns
sort(unique(comm_v5$Lifestage))
sort(unique(comm_v5$Stage.Gut))
sort(unique(comm_v5$Gut.Region))
sort(unique(comm_v5$Sex))

# Save this out for beta diversity calculation in another script
write.csv(x = comm_v5, row.names = F, na = '',
          file = file.path("data", "tidy_data", "beta-diversity-data.csv"))

## ------------------------------------------- ##
     # P2-3: Calculate Alpha Diversity ####
## ------------------------------------------- ##
# Calculate each index
## NOTE: takes a few seconds to complete
dive_indexes <- comm_v5 %>%
  # The following calculations should be done for every row
  rowwise() %>%
  # Calculate the indexes
  dplyr::mutate(Shannon = vegan::diversity(dplyr::across(-Sample.ID:-Sex),
                                           index = 'shannon'),
                Simpson = vegan::diversity(dplyr::across(-Sample.ID:-Sex),
                                           index = 'simpson'),
                Richness = vegan::specnumber(dplyr::across(-Sample.ID:-Sex)),
                Pielous = (Shannon / log(Richness)),
                ACE = fossil::ACE(dplyr::across(-Sample.ID:-Sex),
                                  taxa.row = T),
                Chao1 = fossil::chao1(dplyr::across(-Sample.ID:-Sex),
                                      taxa.row = T) ) %>%
  # Remove the actual community (for ease of plotting)
  select(Sample.ID:Sex, Shannon:Chao1)

# Check structure
dplyr::glimpse(dive_indexes)

# This dataframe is ready for use in stats/visualization!
write.csv(x = dive_indexes, row.names = F, na = '',
          file = file.path("data", "tidy_data", "alpha-diversity-data.csv"))

## -------------------------------------------------------------- ##
            # Part 3: Phylogenetic Tree Tidying ####
## -------------------------------------------------------------- ##
# Look at the tree file
str(tree_data)

# Make a duplicate
tree_v2 <- tree_data

## ------------------------------------------- ##
          # P3-1: Fix Tip Labels ####
## ------------------------------------------- ##
# Tip labels are named with FeatureID which is not informative
  ## Want to replace with our tidied "Taxon" column
tree_v2$tip.label.actual <- tax_v5$Taxon[match(tree_data$tip.label, tax_v5$Feature.ID)]

# Check it out
str(tree_v2)

# Remove the old tip label
tree_v2$tip.label <- NULL

# Re-name the new tip label so later functions don't get suprised and frightened
tree_v2$tip.label <- tree_v2$tip.label.actual
tree_v2$tip.label.actual <- NULL

# Double check what we're left with
str(tree_v2)

# Save it
ape::write.tree(phy = tree_v2, file = file.path("data", "tidy_data", "unrooted-tree.tre"))

## -------------------------------------------------------------- ##
              # Part 4: Taxon Abundance Tidying ####
## -------------------------------------------------------------- ##
## ------------------------------------------- ##
    # P4-1: Migrate Phylum and Family ####
## ------------------------------------------- ##
# Check structure of taxon metadata
dplyr::glimpse(tax_v5)

# And the community data
dplyr::glimpse(comm_v5)

# Create the necessary community dataset
comm_v6 <- comm_v5 %>%
  # Pivot to long format
  tidyr::pivot_longer(cols = -Sample.ID:-Sex,
                      names_to = "Taxon",
                      values_to = "Abundance") %>%
  # Bring over all of the columns from the taxonomic data
  left_join(tax_v5, by = "Taxon") %>%
  # Keep only desired columns (non-specified columns are dropped)
  dplyr::select(Sample.ID, Lifestage, Gut.Region, Stage.Gut, Sex,
                Domain, Phylum, Family, Abundance)

# Check what we've created
dplyr::glimpse(comm_v6)

## ------------------------------------------- ##
   # P4-2: Check Phylum & Family Columns ####
## ------------------------------------------- ##
# Check domain (shouldn't have issues but better safe than sorry)
sort(unique(comm_v6$Domain))

# Check phylum
sort(unique(comm_v6$Phylum))
  ## "RsaHF231" is a candidate phylum so we'll allow it (for now)

# Check family
sort(unique(comm_v6$Family))
  ## Some worrisome so let's get more detail on those

# Make a vector of potentially incorrect family IDs
check_fam <- c("0319-6G20", "67-14", "A0839", "A4b", "AKYG1722", "bacteriap25",
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
(dubious_fams <- tax_v5 %>%
  dplyr::filter(Family %in% check_fam))

# Because it's abundance, we can leave this alone (though may need to return)

## ------------------------------------------- ##
      # P4-3: Family-Level Abundance ####
## ------------------------------------------- ##
# Want family-level abundance
fam_abun <- comm_v6 %>%
  # Sum abundance within family (retaining other, coarser groupings)
  dplyr::group_by(Sample.ID, Lifestage, Gut.Region, Stage.Gut, Sex,
                  Domain, Phylum, Family) %>%
  dplyr::summarise(Abundance = sum(Abundance, na.rm = T)) %>%
  dplyr::ungroup() %>%
  # Remove any NAs and 0 abundances
  dplyr::filter(!is.na(Family) & Abundance > 0) %>%
  # Also want to calculate relative abundance
    ## "relative" means within sample calculations are needed
  dplyr::group_by(Sample.ID) %>%
  dplyr::mutate(totalAbun = sum(Abundance, na.rm = T)) %>%
  dplyr::ungroup() %>%
  # Calculate relative abundance
  dplyr::mutate(relativeAbun = (Abundance / totalAbun) * 100)

# Glimpse this
dplyr::glimpse(fam_abun)

# Should have lost *a lot* of rows
nrow(comm_v6) - nrow(fam_abun)
dim(comm_v6); dim(fam_abun)

# Save this file for later use
write.csv(x = fam_abun, row.names = F, na = '',
          file = file.path("data", "tidy_data", "family-abun.csv"))

## ------------------------------------------- ##
      # P4-4: Phylum-Level Abundance ####
## ------------------------------------------- ##
# Want phylum-level abundance
phyla_abun <- comm_v6 %>%
  # Sum abundance within phylum (retaining other, coarser groupings)
  dplyr::group_by(Sample.ID, Lifestage, Gut.Region, Stage.Gut, Sex,
                  Domain, Phylum) %>%
  dplyr::summarise(Abundance = sum(Abundance, na.rm = T)) %>%
  dplyr::ungroup() %>%
  # Remove any NAs and 0 abundances
  dplyr::filter(!is.na(Phylum) & Abundance > 0) %>%
  # Also want to calculate relative abundance
  ## "relative" means within sample calculations are needed
  dplyr::group_by(Sample.ID) %>%
  dplyr::mutate(totalAbun = sum(Abundance, na.rm = T)) %>%
  dplyr::ungroup() %>%
  # Calculate relative abundance
  dplyr::mutate(relativeAbun = (Abundance / totalAbun) * 100)

# Glimpse this
dplyr::glimpse(phyla_abun)

# Should have lost *a lot* of rows
nrow(comm_v6) - nrow(phyla_abun)
dim(comm_v6); dim(phyla_abun)

# Save this file for later use
write.csv(x = phyla_abun, row.names = F, na = '',
          file = file.path("data", "tidy_data", "phylum-abun.csv"))

## ------------------------------------------- ##
       # P4-5: Genus-Level Abundance ####
## ------------------------------------------- ##
# Need to create the necessary community dataset
gen_abun <- comm_v5 %>%
  # Pivot to long format
  tidyr::pivot_longer(cols = -Sample.ID:-Sex,
                      names_to = "Taxon",
                      values_to = "Abundance") %>%
  # Bring over all of the columns from the taxonomic data
  dplyr::left_join(tax_v5, by = "Taxon") %>%
  # Get generic abundance totals
  dplyr::group_by(Sample.ID, Lifestage, Gut.Region, Stage.Gut, Sex, Genus) %>%
  dplyr::summarise(Abundance = sum(Abundance, na.rm = T)) %>%
  dplyr::ungroup() %>%
  # Remove rows without a genus-level ID / those with abundance of 0
  filter(!is.na(Genus) & Abundance > 0)

# Check that out
dplyr::glimpse(gen_abun)

# Save it
write.csv(x = gen_abun, row.names = F, na = '',
          file = file.path("data", "tidy_data", "genus-abun.csv"))

## -------------------------------------------------------------- ##
             # Part 5: Calculate Unifrac Distances ####
## -------------------------------------------------------------- ##
# Re-check what you may need for unifrac distance calculation
## Abundance table
str(comm_v6)
class(comm_v6)

## Phylogenetic tree
str(tree_v2)
class(tree_v2)

## ------------------------------------------- ##
         # P5-1: Weighted Unifrac ####
## ------------------------------------------- ##
# https://rdrr.io/github/MadsAlbertsen/ampvis2/man/unifrac.html
# unifrac::unifrac(tree_v2)



# END ####

