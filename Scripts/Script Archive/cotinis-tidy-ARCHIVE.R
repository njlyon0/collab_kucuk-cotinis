## ------------------------------------------- ##
          # Get Meaningful Taxon ID ####
## ------------------------------------------- ##

# This approach was abandoned because there are more FeatureIDs than taxa
  ## This make sense because FeatureIDs are genetic
  ## So I moved this out in case I want it later
          
          
#names(comm.data)
comm.data.v2 <- comm.data %>%
  ## Create a column from the rownames
  dplyr::mutate(Feature.ID = rownames(comm.data)) %>%
  ## Left join the taxonomy info by that column
  left_join(tax.data, by = "Feature.ID") %>%
  ## Remove unfed individuals, now-unuecessary feature ID, and confidence cols
  select(-UnfedFemale1Mid:-Unfedfemale1mid, -Feature.ID, -Confidence) %>%
  ## Go to long format
  pivot_longer(!Taxon, names_to = "Sample.ID", values_to = "Count") %>%
  ## Some of those samples are not unique
  ### Group by sample ID
  group_by(Sample.ID, Taxon) %>%
  ### Summarise through those non-unique values
  dplyr::summarise(Count = sum(Count)) %>%
  ## Return a dataframe
  as.data.frame()

# Split the Taxon by the semicolons (it is too dense as is)
comm.data.v3 <- comm.data.v2 %>%
  separate(col = Taxon,
           into = c("Domain", "Phylum", "Subphylum",
                    "Order", "Family", "Genus",
                    "Species"),
           sep = ";D_", remove = T)

# Examine each of the new columns
sort(unique(comm.data.v3$Domain)) # starts with "D_0__"
sort(unique(comm.data.v3$Phylum)) # starts with "1__"
sort(unique(comm.data.v3$Subphylum)) # starts with "2__"
sort(unique(comm.data.v3$Order)) # starts with "3__"
sort(unique(comm.data.v3$Family)) # starts with "4__"
sort(unique(comm.data.v3$Genus)) # starts with "5__"
sort(unique(comm.data.v3$Species)) # starts with "6__"

# Remove the preceeding numbers
comm.data.v3$Domain <- gsub("D_0__", "", comm.data.v3$Domain)
sort(unique(comm.data.v3$Domain)) # starts with "D_0__"
comm.data.v3$Phylum <- gsub("1__", "", comm.data.v3$Phylum)
sort(unique(comm.data.v3$Phylum)) # starts with "1__"
comm.data.v3$Subphylum <- gsub("2__", "", comm.data.v3$Subphylum)
sort(unique(comm.data.v3$Subphylum)) # starts with "2__"
comm.data.v3$Order <- gsub("3__", "", comm.data.v3$Order)
sort(unique(comm.data.v3$Order)) # starts with "3__"
comm.data.v3$Family <- gsub("4__", "", comm.data.v3$Family)
sort(unique(comm.data.v3$Family)) # starts with "4__"
comm.data.v3$Genus <- gsub("5__", "", comm.data.v3$Genus)
sort(unique(comm.data.v3$Genus)) # starts with "5__"
comm.data.v3$Species <- gsub("6__", "", comm.data.v3$Species)
sort(unique(comm.data.v3$Species)) # starts with "6__"

# How many NAs in the lower classifications?
plyr::count(is.na(comm.data.v3$Species))
plyr::count(is.na(comm.data.v3$Genus))
plyr::count(is.na(comm.data.v3$Family))
plyr::count(is.na(comm.data.v3$Order))
plyr::count(is.na(comm.data.v3$Subphylum))
plyr::count(is.na(comm.data.v3$Phylum))
plyr::count(is.na(comm.data.v3$Domain))

# Re-create a better Taxon ID
names(comm.data.v3)
comm.data.v3$Taxon <- with(comm.data.v3, paste0(
  Domain, Phylum, Subphylum, Order,
  Family, " (", Genus, " ", Species, ")"
))
sort(unique(comm.data.v3$Taxon))

# Enough unique IDs?
## Need as many unique IDs as the original taxon column had
length(unique(comm.data.v2$Taxon))
length(unique(comm.data.v3$Taxon))

# Need to get back to wide format
comm.data.v4 <- comm.data.v3 %>%
  ## Remove now-unnecessary columns
  select(-Domain:-Species) %>%
  ## Pivot back to wide format, this time with Taxon as the columns
  pivot_wider(id_cols = Sample.ID,
              names_from = Taxon,
              values_from = Count) %>%
  ## Return dataframe
  as.data.frame()


## ------------------------------------------- ##
# Ensure No Data is Lost ####
## ------------------------------------------- ##
# Those steps could have lost data which would be bad
## We need to check this possibility

# Interestingly, there are far more FeatureIDs than there are taxonomic IDs
## That's why our summarise step (line 62) is needed
length(unique(tax.data$Feature.ID))
length(unique(tax.data$Taxon))

# Which ones are non-unique?
DupCheck <- plyr::count(as.factor(tax.data$Taxon))
DupCheck <- filter(DupCheck, freq > 1)
#View(DupCheck)

# Compare dimensions
dim(DimCheck)
dim(comm.data.v4)  

names(comm.data.v2)
names(comm.data.v2[1:5])

comm.data.v3 <- t((comm.data.v2))
View(comm.data.v3)
