## -------------------------------------------------------------- ##
              # Kucuk Cotinis Project - Analysis Code
## -------------------------------------------------------------- ##
# Code written by Nicholas J Lyon

# PURPOSE:
  ## Test the three hypotheses posed by Roy (Kucuk)
  ### 1) bacterial communities differ between adults and larvae
  ### 2) bacterial communities differ among gut stages
  ### 3) bacterial communities differ between sexes
  #### Note that 1 & 2 are tested together (see paper methods for explanation)

# Clear environment (always better to start with tabula rasa)
rm(list = ls())

# Set working directory
getwd() # should end in ".../Kucuk-CotinisCollab"
myWD <- getwd()

# Necessary libraries
library(RRPP); library(tidyverse); library(vegan); library(ape)

## -------------------------------------------------------------- ##
                # Data Retrieval & Housekeeping ####
## -------------------------------------------------------------- ##
# Retrieve the relevant datasets
alpha <- read.csv("./Data/Tidy Data/alpha-diversity-data.csv")
beta <- read.csv("./Data/Tidy Data/beta-diversity-data.csv")
fams <- read.csv("Data/Tidy Data/family-abun.csv")
phyla <- read.csv("Data/Tidy Data/phylum-abun.csv")
wtd.frc.dist <- read.csv("./Data/wtd-unfrc-dist.csv")[-1]
uwt.frc.dist <- read.csv("./Data/unwtd-unfrc-dist.csv")[-1]

# Look at them to be sure nothing obvious is wrong
str(alpha)
str(beta)
str(fams)
str(phyla)
str(wtd.frc.dist)

# Calculate the Bray Curtis and Jaccard distances
bc.dist <- vegdist(beta[-c(1:5)], method = 'bray')
jac.dist <- vegdist(beta[-c(1:5)], method = 'jaccard')

# We also want to order the factor levels of what we're using as groups for PCoAs
unique(beta$Stage.Gut)
beta$Stage.Gut <- factor(beta$Stage.Gut, 
                         levels = c("Larval paunch", "Larval ileum", "Larval midgut",
                                    "Adult midgut", "Adult hindgut"))
unique(beta$Stage.Gut)

# Strip out the ID missing from the unifrac data
frc.data <- beta %>%
  filter(Sample.ID != "Amid48") %>%
  as.data.frame()

# Create an adults-only subset of alpha diversity data
alpha.sex <- alpha %>%
  filter(Sex == "female" | Sex == "male") %>%
  as.data.frame()

# Do the same for beta diversity
beta.sex <- beta %>%
  filter(Sex == "female" | Sex == "male") %>%
  as.data.frame()

## -------------------------------------------------------------- ##
                  # Life Stage / Gut Region Tests ####
## -------------------------------------------------------------- ##
## ------------------------------------------- ##
          # Alpha Div. ~ Life/Gut ####
## ------------------------------------------- ##
# Incl. Chao, ACE, Simpson, and Shannon diversity
str(alpha)

# Chao1 Diversity
  ## Fit the model
chao <- lm.rrpp(Chao1 ~ Stage.Gut, data = alpha, iter = 9999)
  ## Results?
anova(chao)
  ## If significant, run pairwise comparisons
summary(pairwise(chao, groups = alpha$Stage.Gut))

# Shannon Diversity
shan <- lm.rrpp(Shannon ~ Stage.Gut, data = alpha, iter = 9999)
anova(shan)
summary(pairwise(shan, groups = alpha$Stage.Gut))

# Simpson Diversity
simp <- lm.rrpp(Simpson ~ Stage.Gut, data = alpha, iter = 9999)
anova(simp)
summary(pairwise(simp, groups = alpha$Stage.Gut))

# Pielou's Evenness
piel <- lm.rrpp(Pielous ~ Stage.Gut, data = alpha, iter = 9999)
anova(piel)
summary(pairwise(piel, groups = alpha$Stage.Gut))

# ACE
ace <- lm.rrpp(ACE ~ Stage.Gut, data = alpha, iter = 9999)
anova(ace)
summary(pairwise(ace, groups = alpha$Stage.Gut))

## ------------------------------------------- ##
            # Beta Div. ~ Life/Gut ####
## ------------------------------------------- ##
# Get a version of the beta diversity community that lacks all other values
beta.resp <- as.matrix(select(beta, -Sample.ID:-Sex))

# Now plug that whole dataframe in as the response variable
beta.test <- lm.rrpp(beta.resp ~ Stage.Gut, data = beta)

# Look at the test
anova(beta.test)

# Get pairwise results
summary(pairwise(beta.test, groups = alpha$Stage.Gut))

## -------------------------------------------------------------- ##
                          # Sex Tests ####
## -------------------------------------------------------------- ##
## ------------------------------------------- ##
            # Alpha Div. ~ Sex ####
## ------------------------------------------- ##
# Incl. Chao, ACE, Simpson, and Shannon diversity
str(alpha.sex)

# Chao1 Diversity
anova(lm.rrpp(Chao1 ~ Sex, data = alpha.sex, iter = 9999))

# Shannon Diversity
anova(lm.rrpp(Shannon ~ Sex, data = alpha.sex, iter = 9999))

# Simpson Diversity
anova(lm.rrpp(Simpson ~ Sex, data = alpha.sex, iter = 9999))

# Pielou's Evenness
anova(lm.rrpp(Pielous ~ Sex, data = alpha.sex, iter = 9999))

# ACE
anova(lm.rrpp(ACE ~ Sex, data = alpha.sex, iter = 9999))

## ------------------------------------------- ##
              # Beta Div. ~ Sex ####
## ------------------------------------------- ##
# Get a version of the beta diversity community that lacks all other values
beta.sex.resp <- as.matrix(select(beta.sex, -Sample.ID:-Sex))

# Now plug that whole dataframe in as the response variable
beta.sex.test <- lm.rrpp(beta.sex.resp ~ Sex, data = beta.sex)

# Look at the test
anova(beta.sex.test)

## -------------------------------------------------------------- ##
                    # Family Abundance Tests ####
## -------------------------------------------------------------- ##
# Create relevant subsets
amid.fam <- fams %>%
  filter(Stage.Gut == "Adult midgut") %>%
  as.data.frame()
ahind.fam <- fams %>%
  filter(Stage.Gut == "Adult hindgut") %>%
  as.data.frame()
lmid.fam <- fams %>%
  filter(Stage.Gut == "Larval midgut") %>%
  as.data.frame()
paunch.fam <- fams %>%
  filter(Stage.Gut == "Larval paunch") %>%
  as.data.frame()
ileum.fam <- fams %>%
  filter(Stage.Gut == "Larval ileum") %>%
  as.data.frame()

# Among family abundance differences for:

# Adult midgut
anova(lm.rrpp(Abundance ~ Family, data = amid.fam, iter = 9999))
  ## NS

# Adult hindgut
anova(lm.rrpp(Abundance ~ Family, data = ahind.fam, iter = 9999))
  ## sig

# Larval midgut
anova(lm.rrpp(Abundance ~ Family, data = lmid.fam, iter = 9999))
  ## sig

# Larval Paunch
anova(lm.rrpp(Abundance ~ Family, data = paunch.fam, iter = 9999))
  ## sig

# Larval ileum
anova(lm.rrpp(Abundance ~ Family, data = ileum.fam, iter = 9999))
  ## sig

## -------------------------------------------------------------- ##
                    # Phylum Abundance Tests ####
## -------------------------------------------------------------- ##
# Create relevant subsets
amid.phyl <- phyla %>%
  filter(Stage.Gut == "Adult midgut") %>%
  as.data.frame()
ahind.phyl <- phyla %>%
  filter(Stage.Gut == "Adult hindgut") %>%
  as.data.frame()
lmid.phyl <- phyla %>%
  filter(Stage.Gut == "Larval midgut") %>%
  as.data.frame()
paunch.phyl <- phyla %>%
  filter(Stage.Gut == "Larval paunch") %>%
  as.data.frame()
ileum.phyl <- phyla %>%
  filter(Stage.Gut == "Larval ileum") %>%
  as.data.frame()

# Among phylum abundance differences for:

# Adult midgut
anova(lm.rrpp(Abundance ~ Phylum, data = amid.phyl, iter = 9999))
  ## sig

# Adult hindgut
anova(lm.rrpp(Abundance ~ Phylum, data = ahind.phyl, iter = 9999))
  ## sig

# Larval midgut
anova(lm.rrpp(Abundance ~ Phylum, data = lmid.phyl, iter = 9999))
  ## sig

# Larval Paunch
anova(lm.rrpp(Abundance ~ Phylum, data = paunch.phyl, iter = 9999))
  ## sig

# Larval ileum
anova(lm.rrpp(Abundance ~ Phylum, data = ileum.phyl, iter = 9999))
  ## sig

  

# END ####

