## ---------------------------------------------------------------------------------- ##
                  # Kucuk Cotinis Project - Visualization Code
## ---------------------------------------------------------------------------------- ##
# Code written by Nicholas J Lyon

# PURPOSE:
  ## Do the analysis that Roy (Kucuk) wants

# Clear environment (always better to start with tabula rasa)
rm(list = ls())

# Set working directory
setwd("~/Documents/_Publications/2022_Kucuk_Cotinis/Kucuk-CotinisCollab")

# Necessary libraries
library(RRPP); library(tidyverse); library(vegan); library(ape)

## ------------------------------------------- ##
      # Data Retrieval & Housekeeping ####
## ------------------------------------------- ##


## ------------------------------------------- ##
          # Beta Diversity Tests ####
## ------------------------------------------- ##
# Incl. Bray Curtis dissimilarity and Jaccard distance


# Diversity ~ Life Stage + Gut Section, data = all

# Initial perMANOVA
anova(lm.rrpp(resp ~ factor, data = ref), effect.type = "F")
## Interpretation: at least one group is significantly different from the others
## Yes these groups don't mean anything, but it's still illustrative.

# In reality would need to track down which groups are different from which (pairwise comparisons)
## but we'll leave that alone here



## ------------------------------------------- ##
          # Alpha Diversity Tests ####
## ------------------------------------------- ##
# Incl. Chao, ACE, Simpson, and Shannon diversity


# Diversity ~ Life Stage + Gut Section, data = all




## ------------------------------------------- ##
            # Larval Sex Tests ####
## ------------------------------------------- ##
# Same as above beta/alpha tests but subsetted for only larvae and incl. sex as explanatory variable


# Diversity ~ Life Stage + Gut Section + Sex, data = larvae


#END ####

