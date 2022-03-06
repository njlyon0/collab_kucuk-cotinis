## -------------------------------------------------------------- ##
           # Kucuk Cotinis Project - Visualization Code
## -------------------------------------------------------------- ##
# Code written by Nicholas J Lyon

# PURPOSE:
  ## Perform data visualization of the gut microbe communities found in Cotinis nitida

# Clear environment (always better to start with tabula rasa)
rm(list = ls())

# Set working directory
getwd() # should end in ".../Kucuk-CotinisCollab"
myWD <- getwd()

# Necessary libraries
library(tidyverse); library(vegan); library(ape); library(Rmisc)

## -------------------------------------------- ##
      # Data Retrieval & Housekeeping ####
## -------------------------------------------- ##
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

# Roy (Kucuk) also wants some plots/ordinations that include larvae only
alpha.larv <- filter(alpha, Lifestage == "larva")

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
jacc.pnts <- ape::pcoa(jac.dist)
wtd.frc.pnts <- ape::pcoa(wtd.frc.dist)
uwt.frc.pnts <- ape::pcoa(uwt.frc.dist)

# Now make the ordination so that you can look at it (i.e., so that it prints to the viewer in R)
pcoa.5.ord(mod = bc.pnts, groupcol = beta$Stage.Gut,
           g1 = "Larval paunch", g2 = "Larval ileum", g3 = "Larval midgut",
           g4 = "Adult midgut", g5 = "Adult hindgut",
           legcont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
           legpos = "bottomleft")

# Once that looks good save it out using this dev.off bit
tiff(file = "./Graphs/PCoA-bray-curtis.tiff")
pcoa.5.ord(mod = bc.pnts, groupcol = beta$Stage.Gut,
           g1 = "Larval paunch", g2 = "Larval ileum", g3 = "Larval midgut",
           g4 = "Adult midgut", g5 = "Adult hindgut",
           legcont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
           legpos = "bottomleft")
dev.off()

# Now make the Jaccard PCoA
pcoa.5.ord(mod = jacc.pnts, groupcol = beta$Stage.Gut,
           g1 = "Larval paunch", g2 = "Larval ileum", g3 = "Larval midgut",
           g4 = "Adult midgut", g5 = "Adult hindgut",
           legcont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
           legpos = "bottomleft")

# Once that looks good save it out using this dev.off bit
tiff(file = "./Graphs/PCoA-jaccard.tiff")
pcoa.5.ord(mod = jacc.pnts, groupcol = beta$Stage.Gut,
           g1 = "Larval paunch", g2 = "Larval ileum", g3 = "Larval midgut",
           g4 = "Adult midgut", g5 = "Adult hindgut",
           legcont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
           legpos = "bottomleft")
dev.off()

# Now make the Weighted Unifrac PCoA
pcoa.5.ord(mod = wtd.frc.pnts, groupcol = frc.data$Stage.Gut,
           g1 = "Larval paunch", g2 = "Larval ileum", g3 = "Larval midgut",
           g4 = "Adult midgut", g5 = "Adult hindgut",
           legcont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
           legpos = "topright")

# Save the ordination
tiff(file = "./Graphs/PCoA-weighted-unifrac.tiff")
pcoa.5.ord(mod = wtd.frc.pnts, groupcol = frc.data$Stage.Gut,
           g1 = "Larval paunch", g2 = "Larval ileum", g3 = "Larval midgut",
           g4 = "Adult midgut", g5 = "Adult hindgut",
           legcont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
           legpos = "bottomright")
dev.off()

# Finally, do the Unweighted Unifrac PCoA
pcoa.5.ord(mod = uwt.frc.pnts, groupcol = frc.data$Stage.Gut,
           g1 = "Larval paunch", g2 = "Larval ileum", g3 = "Larval midgut",
           g4 = "Adult midgut", g5 = "Adult hindgut",
           legcont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
           legpos = "topleft")

# Save the ordination
tiff(file = "./Graphs/PCoA-unweighted-unifrac.tiff")
pcoa.5.ord(mod = uwt.frc.pnts, groupcol = frc.data$Stage.Gut,
           g1 = "Larval paunch", g2 = "Larval ileum", g3 = "Larval midgut",
           g4 = "Adult midgut", g5 = "Adult hindgut",
           legcont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
           legpos = "topleft")
dev.off()

## -------------------------------------------- ##
        # Alpha Diversity Bar Graphs ####
## -------------------------------------------- ##
# Incl. Chao, ACE, Simpson, and Shannon diversity & Pielou's Evenness
  ## (Diversity ~ Life Stage + Gut Section, data = all)

# Summarize each of the response variables
shan <- summarySE(data = alpha, measurevar = "Shannon", groupvars = c("Stage.Gut"))
piel <- summarySE(data = alpha, measurevar = "Pielous", groupvars = c("Stage.Gut"))
simp <- summarySE(data = alpha, measurevar = "Simpson", groupvars = c("Stage.Gut"))
ace <- summarySE(data = alpha, measurevar = "ACE", groupvars = c("Stage.Gut"))
chao <- summarySE(data = alpha, measurevar = "Chao1", groupvars = c("Stage.Gut"))

# Re-level the factor column
  ## Pref factor order
stage.gut.lvls <- c("Adult hindgut", "Adult midgut", "Larval midgut",
                    "Larval ileum", "Larval paunch")
  ## Actual re-leveling
shan$Stage.Gut <- factor(shan$Stage.Gut, levels = stage.gut.lvls)
piel$Stage.Gut <- factor(piel$Stage.Gut, levels = stage.gut.lvls)
simp$Stage.Gut <- factor(simp$Stage.Gut, levels = stage.gut.lvls)
chao$Stage.Gut <- factor(chao$Stage.Gut, levels = stage.gut.lvls)
ace$Stage.Gut <- factor(ace$Stage.Gut, levels = stage.gut.lvls)

# Shannon diversity bar graph
ggplot(shan, aes(y = Shannon, x = Stage.Gut,
                fill = Stage.Gut, color = rep('dummy', 5))) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Shannon + se, ymin = Shannon - se), width = 0.2) +
  geom_text(label = "b", x = 0.7, y =  3.8, size = 5) +
  geom_text(label = "a", x = 1.7, y =  2.3, size = 5) +
  geom_text(label = "c", x = 2.7, y =  6.15, size = 5) +
  geom_text(label = "bc", x = 3.7, y = 4.75, size = 5) +
  geom_text(label = "c", x = 4.7, y =  5.45, size = 5) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "Shannon Diversity") +
  pref_theme + theme(legend.position = 'none')
ggsave(plot = last_plot(), "./Graphs/bar-shannon.pdf",
       width = 5, height = 5, units = 'in')

# Pielou's evenness bar graph
ggplot(piel, aes(y = Pielous, x = Stage.Gut,
                 fill = Stage.Gut, color = rep('dummy', 5))) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Pielous + se, ymin = Pielous - se), width = 0.2) +
  geom_text(label = "b", x = 0.7, y =  0.775, size = 5) +
  geom_text(label = "a", x = 1.7, y =  0.65, size = 5) +
  geom_text(label = "c", x = 2.7, y =  0.97, size = 5) +
  geom_text(label = "b", x = 3.7, y =  0.795, size = 5) +
  geom_text(label = "bc", x = 4.7, y = 0.815, size = 5) +
  ylim(0, 1.05) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "Pielou's Evenness") +
  pref_theme + theme(legend.position = 'none')
ggsave(plot = last_plot(), "./Graphs/bar-pielou.pdf",
       width = 5, height = 5, units = 'in')

# Simpson Diversity
ggplot(simp, aes(y = Simpson, x = Stage.Gut,
                 fill = Stage.Gut, color = rep('dummy', 5))) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Simpson + se, ymin = Simpson - se), width = 0.2) +
  geom_text(label = "b", x = 0.7, y =  0.99, size = 5) +
  geom_text(label = "a", x = 1.7, y =  0.82, size = 5) +
  geom_text(label = "b", x = 2.7, y =  1.05, size = 5) +
  geom_text(label = "b", x = 3.7, y =  1, size = 5) +
  geom_text(label = "b", x = 4.7, y =  1.025, size = 5) +
  ylim(0, 1.1) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "Simpson Diversity") +
  pref_theme + theme(legend.position = 'none')
ggsave(plot = last_plot(), "./Graphs/bar-simpson.pdf",
       width = 5, height = 5, units = 'in')

# ACE Diversity
ggplot(ace, aes(y = ACE, x = Stage.Gut, fill = Stage.Gut,
                color = rep('dummy', 5))) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = ACE + se, ymin = ACE - se), width = 0.2) +
  geom_text(label = "ab", x = 0.7, y = 180, size = 5) +
  geom_text(label = "a", x = 1.7, y =  100, size = 5) +
  geom_text(label = "bc", x = 2.7, y = 740, size = 5) +
  geom_text(label = "b", x = 3.7, y =  590, size = 5) +
  geom_text(label = "c", x = 4.7, y =  975, size = 5) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "ACE") +
  pref_theme + theme(legend.position = 'none')
ggsave(plot = last_plot(), "./Graphs/bar-ace.pdf",
       width = 5, height = 5, units = 'in')

# Chao 1 Diversity
ggplot(chao, aes(y = Chao1, x = Stage.Gut, fill = Stage.Gut,
                color = rep('dummy', 5))) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Chao1 + se, ymin = Chao1 - se), width = 0.2) +
  geom_text(label = "ab", x = 0.7, y = 180, size = 5) +
  geom_text(label = "a", x = 1.7, y =  100, size = 5) +
  geom_text(label = "bc", x = 2.7, y = 740, size = 5) +
  geom_text(label = "b", x = 3.7, y =  590, size = 5) +
  geom_text(label = "c", x = 4.7, y =  975, size = 5) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "Chao1 Index") +
  pref_theme + theme(legend.position = 'none')
ggsave(plot = last_plot(), "./Graphs/bar-chao.pdf",
       width = 5, height = 5, units = 'in')

## -------------------------------------------- ##
        # Family Abundance Bar Graphs ####
## -------------------------------------------- ##
# Summarize relevant subsets
amid.fam <- fams %>%
  filter(Stage.Gut == "Adult midgut") %>%
  summarySE(measurevar = "Abundance", groupvars = c("Stage.Gut", "Family")) %>%
  as.data.frame()
ahind.fam <- fams %>%
  filter(Stage.Gut == "Adult hindgut") %>%
  summarySE(measurevar = "Abundance", groupvars = c("Stage.Gut", "Family")) %>%
  as.data.frame()
lmid.fam <- fams %>%
  filter(Stage.Gut == "Larval midgut") %>%
  summarySE(measurevar = "Abundance", groupvars = c("Stage.Gut", "Family")) %>%
  as.data.frame()
paunch.fam <- fams %>%
  filter(Stage.Gut == "Larval paunch") %>%
  summarySE(measurevar = "Abundance", groupvars = c("Stage.Gut", "Family")) %>%
  as.data.frame()
ileum.fam <- fams %>%
  filter(Stage.Gut == "Larval ileum") %>%
  summarySE(measurevar = "Abundance", groupvars = c("Stage.Gut", "Family")) %>%
  as.data.frame()

# Make and save each graph
  ## Adult midgut
ggplot(amid.fam, aes(y = Abundance, x = reorder(Family, -Abundance),
                     fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Family", y = "Adult Midgut Abundance") +
  pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))
ggsave(plot = last_plot(), "./Graphs/fam-abun-adult-midgut.pdf",
       width = 6, height = 6, units = 'in')

  ## Adult hindgut
ggplot(ahind.fam, aes(y = Abundance, x = reorder(Family, -Abundance),
                      fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Family", y = "Adult Hindgut Abundance") +
  pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))
ggsave(plot = last_plot(), "./Graphs/fam-abun-adult-hindgut.pdf",
       width = 6, height = 6, units = 'in')

  ## Larval midgut
ggplot(lmid.fam, aes(y = Abundance, x = reorder(Family, -Abundance),
                     fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Family", y = "Larval Midgut Abundance") +
  pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))
ggsave(plot = last_plot(), "./Graphs/fam-abun-larval-midgut.pdf",
       width = 6, height = 6, units = 'in')

  ## Larval paunch
ggplot(paunch.fam, aes(y = Abundance, x = reorder(Family, -Abundance),
                       fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Family", y = "Larval Paunch Abundance") +
  pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))
ggsave(plot = last_plot(), "./Graphs/fam-abun-larval-paunch.pdf",
       width = 6, height = 6, units = 'in')

  ## Larval ileum
ggplot(ileum.fam, aes(y = Abundance, x = reorder(Family, -Abundance),
                      fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Family", y = "Larval Ileum Abundance") +
  pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))
ggsave(plot = last_plot(), "./Graphs/fam-abun-larval-ileum.pdf",
       width = 6, height = 6, units = 'in')

## -------------------------------------------- ##
        # Phylum Abundance Bar Graphs ####
## -------------------------------------------- ##
# Summarize relevant subsets
amid.phyl <- phyla %>%
  filter(Stage.Gut == "Adult midgut") %>%
  summarySE(measurevar = "Abundance", groupvars = c("Stage.Gut", "Phylum")) %>%
  as.data.frame()
ahind.phyl <- phyla %>%
  filter(Stage.Gut == "Adult hindgut") %>%
  summarySE(measurevar = "Abundance", groupvars = c("Stage.Gut", "Phylum")) %>%
  as.data.frame()
lmid.phyl <- phyla %>%
  filter(Stage.Gut == "Larval midgut") %>%
  summarySE(measurevar = "Abundance", groupvars = c("Stage.Gut", "Phylum")) %>%
  as.data.frame()
paunch.phyl <- phyla %>%
  filter(Stage.Gut == "Larval paunch") %>%
  summarySE(measurevar = "Abundance", groupvars = c("Stage.Gut", "Phylum")) %>%
  as.data.frame()
ileum.phyl <- phyla %>%
  filter(Stage.Gut == "Larval ileum") %>%
  summarySE(measurevar = "Abundance", groupvars = c("Stage.Gut", "Phylum")) %>%
  as.data.frame()

# Make and save each graph
## Adult midgut
ggplot(amid.phyl, aes(y = Abundance, x = reorder(Phylum, -Abundance),
                      fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Phylum", y = "Adult Midgut Abundance") +
  pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))
ggsave(plot = last_plot(), "./Graphs/phyla-abun-adult-midgut.pdf",
       width = 6, height = 6, units = 'in')

## Adult hindgut
ggplot(ahind.phyl, aes(y = Abundance, x = reorder(Phylum, -Abundance),
                       fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Phylum", y = "Adult Hindgut Abundance") +
  pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))
ggsave(plot = last_plot(), "./Graphs/phyla-abun-adult-hindgut.pdf",
       width = 6, height = 6, units = 'in')

## Larval midgut
ggplot(lmid.phyl, aes(y = Abundance, x = reorder(Phylum, -Abundance),
                      fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Phylum", y = "Larval Midgut Abundance") +
  pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))
ggsave(plot = last_plot(), "./Graphs/phyla-abun-larval-midgut.pdf",
       width = 6, height = 6, units = 'in')

## Larval paunch
ggplot(paunch.phyl, aes(y = Abundance, x = reorder(Phylum, -Abundance),
                        fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Phylum", y = "Larval Paunch Abundance") +
  pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))
ggsave(plot = last_plot(), "./Graphs/phyla-abun-larval-paunch.pdf",
       width = 6, height = 6, units = 'in')

## Larval ileum
ggplot(ileum.phyl, aes(y = Abundance, x = reorder(Phylum, -Abundance),
                       fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Phylum", y = "Larval Ileum Abundance") +
  pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))
ggsave(plot = last_plot(), "./Graphs/phyla-abun-larval-ileum.pdf",
       width = 6, height = 6, units = 'in')

## -------------------------------------------- ##
       # Global Phylum Abundance Plot ####
## -------------------------------------------- ##
# Look at the phylum-level abundance
head(phyla)

# Get a global version of that (i.e., one that ditches the sample information)
phyl.global <- phyla %>%
  # Remove un-needed columns
  dplyr::select(Domain, Phylum, Abundance) %>%
  # Group by phylum
  group_by(Phylum) %>%
  # Sum to get one instance of each
  dplyr::summarise(
    Abundance = sum(Abundance, na.rm = T)
  ) %>%
  as.data.frame()

# Look at what that produced
head(phyl.global)

# Graph that
ggplot(phyl.global, aes(y = Abundance, x = reorder(Phylum, -Abundance),
                       fill = 'x', color = 'x')) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = "#41ab5d") +
  scale_color_manual(values = 'black') +
  labs(x = "Phylum", y = "Abundance") +
  pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))
ggsave(plot = last_plot(), "./Graphs/phyla-abun-all.pdf",
       width = 6, height = 6, units = 'in')

# END ####

