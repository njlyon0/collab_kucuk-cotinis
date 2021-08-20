## ---------------------------------------------------------------------------------- ##
                  # Kucuk Cotinis Project - Visualization Code
## ---------------------------------------------------------------------------------- ##
# Code written by Nicholas J Lyon

# PURPOSE:
  ## Perform data visualization of the gut microbe communities found in Cotinis nitida

# Clear environment (always better to start with tabula rasa)
rm(list = ls())

# Set working directory
setwd("~/Documents/_Publications/2022_Kucuk_Cotinis/Kucuk-CotinisCollab")

# Necessary libraries
library(tidyverse); library(vegan); library(ape); library(Rmisc)

## -------------------------------------------- ##
      # Data Retrieval & Housekeeping ####
## -------------------------------------------- ##


vegdist(bz.patch.rsp, method = "jaccard")

# Retrieve the relevant datasets
bc.data <- read.csv("./Data/Beta Diversity/bc-data.csv")
jacc.data <- read.csv("./Data/Beta Diversity/jacc-data.csv")
uwt.frc.data <- read.csv("./Data/Beta Diversity/unwtd-unfrc-data.csv")
wtd.frc.data <- read.csv("./Data/Beta Diversity/wtd-unfrc-data.csv")
alpha <- read.csv("./Data/Alpha Diversity/dive-data.csv")

# Look at them to be sure nothing obvious is wrong
str(bc.data)
str(jacc.data)
str(uwt.frc.data)
str(wtd.frc.data)
str(alpha)

# The nature of distance *matrices* makes saving them then re-loading them a pain
  ## So, we will quickly load and fix the raw distance matrices for both metrics (BC & Jaccard)

# Bray Curtis dissimilarity matrix
  ## Load qza file
bc.qza <- read_qza("./Data/diversity6/bray_curtis_distance_matrix.qza")
  ## Get just the part of that we want
bc.dist.v0 <- as.matrix(bc.qza$data)
  ## Remove comparisons to unfed individuals (pilot study that is irrelevant here)
bc.dist <- bc.dist.v0[-c(33:36), -c(33:36)]

# Jaccard distance matrix
jacc.qza <- read_qza("./Data/diversity6/jaccard_distance_matrix.qza")
jacc.dist.v0 <- as.matrix(jacc.qza$data)
jacc.dist <- jacc.dist.v0[-c(33:36), -c(33:36)]

# Weighted & Unweighted Unifrac distance matrices
wtd.frc.qza <- read_qza("./Data/diversity6/weighted_unifrac_distance_matrix.qza")
uwt.frc.qza <- read_qza("./Data/diversity6/unweighted_unifrac_distance_matrix.qza")
wtd.frc.v0 <- as.matrix(wtd.frc.qza$data)
uwt.frc.v0 <- as.matrix(uwt.frc.qza$data)
wtd.frc.dist <- wtd.frc.v0[-c(33:36), -c(33:36)]
uwt.frc.dist <- uwt.frc.v0[-c(33:36), -c(33:36)]

# We also want to order the factor levels of what we're using as groups for PCoAs
  ## Bray Curtis
unique(bc.data$Stage.n.Gut)
bc.data$Stage.n.Gut <- factor(bc.data$Stage.n.Gut, 
                              levels = c("Larval paunch", "Larval ileum", "Larval midgut",
                                         "Adult midgut", "Adult hindgut"))
unique(bc.data$Stage.n.Gut)

  ## Jaccard
unique(jacc.data$Stage.n.Gut)
jacc.data$Stage.n.Gut <- factor(jacc.data$Stage.n.Gut, 
                              levels = c("Larval paunch", "Larval ileum", "Larval midgut",
                                         "Adult midgut", "Adult hindgut"))
unique(jacc.data$Stage.n.Gut)

  ## Unweighted Unifrac
unique(uwt.frc.data$Stage.n.Gut)
uwt.frc.data$Stage.n.Gut <- factor(uwt.frc.data$Stage.n.Gut, 
                                levels = c("Larval paunch", "Larval ileum", "Larval midgut",
                                           "Adult midgut", "Adult hindgut"))
unique(uwt.frc.data$Stage.n.Gut)

  ## Weighted Unifrac
unique(wtd.frc.data$Stage.n.Gut)
wtd.frc.data$Stage.n.Gut <- factor(wtd.frc.data$Stage.n.Gut, 
                               levels = c("Larval paunch", "Larval ileum", "Larval midgut",
                                          "Adult midgut", "Adult hindgut"))
unique(wtd.frc.data$Stage.n.Gut)


# Roy (Kucuk) also wants some plots/ordinations that include larvae only
alpha.larv <- filter(alpha, Life.Stage == "larva")

# Retrieve my custom Principal Coordinates Analysis (PCoA) ordination function
pcoa.5.ord <- function(mod, groupcol, g1, g2, g3, g4, g5,
                       lntp1 = 1, lntp2 = 1, lntp3 = 1, lntp4 = 1, lntp5 = 1,
                       legcont, legpos = "topleft") {
  ## mod = object returned by ape::pcoa
  ## groupcol = group column in the dataframe that contains those (not the matrix used in vegdist)
  ## g1 - g5 = how each group appears in your dataframe (in quotes)
  ## lntp1 - 5 = what sort of line each ellipse will be made of (accepts integers between 1 and 6 for diff lines)
  ## legcont = single object for what you want the content of the legend to be
  ## legpos = legend position, either numeric vector of x/y coords or shorthand accepted by "legend" function
  
  # Create plot
  plot(mod$vectors, display = 'sites', choice = c(1, 2), type = 'none',
       xlab = paste0("PC1 (", round(mod$values$Relative_eig[1] * 100, digits = 2), "%)"),
       ylab = paste0("PC2 (", round(mod$values$Relative_eig[2] * 100, digits = 2), "%)"))
  ## Probably want the relative eigenvalues (% variation explained per axis) on the plot in an obvious way
  
  # Set colors (easier for you to modify if we set this now and call these objects later)
  col1 <- "#a50026" # dark red
  col2 <- "#f46d43" # red
  col3 <- "#fee090" # yellow-orange
  col4 <- "#abd9e9" # light blue
  col5 <- "#4575b4" # blue
  
  # Add points for each group with a different color per group
  points(mod$vectors[groupcol == g1, 1], mod$vectors[groupcol == g1, 2], pch = 21, bg = col1)
  points(mod$vectors[groupcol == g2, 1], mod$vectors[groupcol == g2, 2], pch = 22, bg = col2)
  points(mod$vectors[groupcol == g3, 1], mod$vectors[groupcol == g3, 2], pch = 23, bg = col3)
  points(mod$vectors[groupcol == g4, 1], mod$vectors[groupcol == g4, 2], pch = 24, bg = col4)
  points(mod$vectors[groupcol == g5, 1], mod$vectors[groupcol == g5, 2], pch = 25, bg = col5)
  ## As of right now the colors are colorblind safe and each group is also given its own shape
  
  # Get a single vector of your manually set line types for the ellipses
  lntps <- c(lntp1, lntp2, lntp3, lntp4, lntp5)
  
  # Ordinate SD ellipses around the centroid
  vegan::ordiellipse(mod$vectors, groupcol, 
                     col = c(g1 = col1, g2 = col2, g3 = col3, g4 = col4, g5 = col5),
                     display = "sites", kind = "sd", lwd = 2, lty = lntps, label = F)
  
  # Add legend
  legend(legpos, legend = legcont, bty = "n", 
         title = NULL,  cex = 1.15, 
         pch = c(21, 22, 23, 24, 25),
         pt.bg = c(col1, col2, col3, col4, col5))
  
}

# We also want to customize some plotting aesthetics for our ggplot plots that we can do here
all.cols <- c("Larval paunch" = "#a50026", "Larval ileum" = "#f46d43", "Larval midgut" = "#fee090",
              "Adult midgut" = "#abd9e9", "Adult hindgut" = "#4575b4")
sex.shps <- c("male" = 24, "female" = 25, "larva"= 21)
dodge <- position_dodge(width = 0.5)
pref_theme <- theme_classic() +
  theme(legend.title = element_blank(),
        #legend.position = 'none',
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        legend.background = element_rect(fill = alpha('black', 0)),
        axis.text.x = element_text(angle = 35, hjust = 1))

## -------------------------------------------- ##
       # Beta Diversity Ordinations ####
## -------------------------------------------- ##
# Incl. Bray Curtis dissimilarity, Jaccard distance, etc.
  ## (Diversity ~ Life Stage + Gut Section, data = all)

# Run a PCoA on all distance matrices
bc.pnts <- ape::pcoa(bc.dist)
jacc.pnts <- ape::pcoa(jacc.dist)
wtd.frc.pnts <- ape::pcoa(wtd.frc.dist)
uwt.frc.pnts <- ape::pcoa(uwt.frc.dist)

# Now make the ordination so that you can look at it (i.e., so that it prints to the viewer in R)
pcoa.5.ord(mod = bc.pnts, groupcol = bc.data$Stage.n.Gut,
           g1 = "Larval paunch", g2 = "Larval ileum", g3 = "Larval midgut",
           g4 = "Adult midgut", g5 = "Adult hindgut",
           legcont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
           legpos = "bottomleft")

# Once that looks good save it out using this dev.off bit
jpeg(file = "./Graphs/PCoA-bray-curtis.jpg")
pcoa.5.ord(mod = bc.pnts, groupcol = bc.data$Stage.n.Gut,
           g1 = "Larval paunch", g2 = "Larval ileum", g3 = "Larval midgut",
           g4 = "Adult midgut", g5 = "Adult hindgut",
           legcont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
           legpos = "bottomleft")
dev.off()

# Now make the Jaccard PCoA
pcoa.5.ord(mod = jacc.pnts, groupcol = jacc.data$Stage.n.Gut,
           g1 = "Larval paunch", g2 = "Larval ileum", g3 = "Larval midgut",
           g4 = "Adult midgut", g5 = "Adult hindgut",
           legcont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
           legpos = "bottomleft")

# Once that looks good save it out using this dev.off bit
jpeg(file = "./Graphs/PCoA-jaccard.jpg")
pcoa.5.ord(mod = jacc.pnts, groupcol = jacc.data$Stage.n.Gut,
           g1 = "Larval paunch", g2 = "Larval ileum", g3 = "Larval midgut",
           g4 = "Adult midgut", g5 = "Adult hindgut",
           legcont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
           legpos = "bottomleft")
dev.off()

# Now make the Weighted Unifrac PCoA
pcoa.5.ord(mod = wtd.frc.pnts, groupcol = wtd.frc.data$Stage.n.Gut,
           g1 = "Larval paunch", g2 = "Larval ileum", g3 = "Larval midgut",
           g4 = "Adult midgut", g5 = "Adult hindgut",
           legcont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
           legpos = "topright")

# Save the ordination
jpeg(file = "./Graphs/PCoA-weighted-unifrac.jpg")
pcoa.5.ord(mod = wtd.frc.pnts, groupcol = wtd.frc.data$Stage.n.Gut,
           g1 = "Larval paunch", g2 = "Larval ileum", g3 = "Larval midgut",
           g4 = "Adult midgut", g5 = "Adult hindgut",
           legcont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
           legpos = "bottomright")

dev.off()

# Finally, do the Unweighted Unifrac PCoA
pcoa.5.ord(mod = uwt.frc.pnts, groupcol = uwt.frc.data$Stage.n.Gut,
           g1 = "Larval paunch", g2 = "Larval ileum", g3 = "Larval midgut",
           g4 = "Adult midgut", g5 = "Adult hindgut",
           legcont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
           legpos = "topleft")

# Save the ordination
jpeg(file = "./Graphs/PCoA-unweighted-unifrac.jpg")
pcoa.5.ord(mod = uwt.frc.pnts, groupcol = uwt.frc.data$Stage.n.Gut,
           g1 = "Larval paunch", g2 = "Larval ileum", g3 = "Larval midgut",
           g4 = "Adult midgut", g5 = "Adult hindgut",
           legcont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
           legpos = "topleft")

dev.off()

## -------------------------------------------- ##
       # Alpha Diversity Scatterplots ####
## -------------------------------------------- ##
# Incl. Chao, ACE, Simpson, and Shannon diversity & Pielou's Evenness
  ## (Diversity ~ Life Stage + Gut Section, data = all)

# Shannon diversity
  ## Plot
ggplot(alpha, aes(y = shannon, x = Stage.n.Gut, fill = Stage.n.Gut)) +
  geom_jitter(pch = 24, size = 2.5, width = 0.2) +
  scale_fill_manual(values = all.cols) +
  labs(x = "Life Stage & Gut Region", y = "Shannon Diversity") +
  pref_theme +
  theme(legend.position = c(0.85, 0.2))

  ## Save plot
ggsave(plot = last_plot(), "./Graphs/scatter-shannon.pdf", width = 5, height = 5, units = 'in')

# Pielou's Evenness
  ## Plot
ggplot(alpha, aes(y = pielou_e, x = Stage.n.Gut, fill = Stage.n.Gut)) +
  geom_jitter(pch = 24, size = 2.5, width = 0.2) +
  scale_fill_manual(values = all.cols) +
  labs(x = "Life Stage & Gut Region", y = "Pielou's Evenness") +
  pref_theme +
  theme(legend.position = c(0.85, 0.2))

  ## Save plot
ggsave(plot = last_plot(), "./Graphs/scatter-pielou.pdf", width = 5, height = 5, units = 'in')

## -------------------------------------------- ##
        # Alpha Diversity Bar Graphs ####
## -------------------------------------------- ##
# Incl. Chao, ACE, Simpson, and Shannon diversity & Pielou's Evenness
  ## (Diversity ~ Life Stage + Gut Section, data = all)

# Summarize each of the response variables
shan <- summarySE(data = alpha, measurevar = "shannon", groupvars = c("Stage.n.Gut"))
piel <- summarySE(data = alpha, measurevar = "pielou_e", groupvars = c("Stage.n.Gut"))

# Shannon diversity bar graph
ggplot(shan, aes(y = shannon, x = Stage.n.Gut,
                fill = Stage.n.Gut, color = rep('dummy', nrow(shan)))) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = shannon + se, ymin = shannon - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "Shannon Diversity") +
  pref_theme + theme(legend.position = 'none')
ggsave(plot = last_plot(), "./Graphs/bar-shannon.pdf", width = 5, height = 5, units = 'in')

# Pielou's evenness bar graph
ggplot(piel, aes(y = pielou_e, x = Stage.n.Gut,
                 fill = Stage.n.Gut, color = rep('dummy', nrow(piel)))) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = pielou_e + se, ymin = pielou_e - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "Pielou's Evenness") +
  pref_theme + theme(legend.position = 'none')
ggsave(plot = last_plot(), "./Graphs/bar-pielou.pdf", width = 5, height = 5, units = 'in')
  

## -------------------------------------------- ##
 # Larval Sex Ordinations & Bar/Scatterplots ####
## -------------------------------------------- ##
# Same as above beta/alpha tests but subsetted for only larvae and incl. sex as explanatory variable
  ## Diversity ~ Life Stage + Gut Section + Sex, data = larvae

# Issue here:
  ## Larvae don't have an identified sex in the metadata
  ## So I can't make these plots...

#END ####

