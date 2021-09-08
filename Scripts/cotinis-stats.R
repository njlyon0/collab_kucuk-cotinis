## ---------------------------------------------------------------------------------- ##
                  # Kucuk Cotinis Project - Visualization Code
## ---------------------------------------------------------------------------------- ##
# Code written by Nicholas J Lyon

# PURPOSE:
  ## Do the analysis that Roy (Kucuk) wants

# Clear environment (always better to start with tabula rasa)
rm(list = ls())

# Set working directory
getwd() # should end in ".../Kucuk-CotinisCollab"
myWD <- getwd()

# Necessary libraries
library(RRPP); library(tidyverse); library(vegan); library(ape)

## ------------------------------------------- ##
      # Data Retrieval & Housekeeping ####
## ------------------------------------------- ##
# Retrieve the relevant datasets
alpha <- read.csv("./Data/alpha-diversity-data.csv")
beta <- read.csv("./Data/beta-diversity-data.csv")
wtd.frc.dist <- read.csv("./Data/wtd-unfrc-dist.csv")[-1]
uwt.frc.dist <- read.csv("./Data/unwtd-unfrc-dist.csv")[-1]

# Look at them to be sure nothing obvious is wrong
str(alpha)
str(beta)
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

## ------------------------------------------- ##
        # Alpha Diversity Tests ####
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
          # Beta Diversity Tests ####
## ------------------------------------------- ##
# Get a version of the beta diversity community that lacks all other values
beta.resp <- as.matrix(select(beta, -Sample.ID:-Sex))

# Now plug that whole dataframe in as the response variable
beta.test <- lm.rrpp(beta.resp ~ Stage.Gut, data = beta)

# Look at the test
anova(beta.test)

# Get pairwise results
summary(pairwise(beta.test, groups = alpha$Stage.Gut))

#END ####

