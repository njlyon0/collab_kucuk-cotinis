## -------------------------------------------------------------- ##
         # Kucuk Cotinis Project - Figure Creation Code
## -------------------------------------------------------------- ##
# Code written by Nicholas J Lyon

# PURPOSE:
## Create figures for the gut microbe communities found in Cotinis nitida

# Clear environment (always better to start with tabula rasa)
rm(list = ls())

# Necessary libraries
# install.packages("librarian")
librarian::shelf(tidyverse, vegan, psych, ape, cowplot, supportR, ggvenn)

## -------------------------------------------------------------- ##
              # Data Retrieval & Housekeeping ####
## -------------------------------------------------------------- ##
# Retrieve the relevant datasets
alpha <- read.csv(file.path("Data", "Tidy Data", "alpha-diversity-data.csv"))
beta <- read.csv(file.path("Data", "Tidy Data", "beta-diversity-data.csv"))
fams <- read.csv(file.path("Data", "Tidy Data", "family-abun.csv"))
phyla <- read.csv(file.path("Data", "Tidy Data", "phylum-abun.csv"))

# Look at them to be sure nothing obvious is wrong
str(alpha)
str(fams)
str(phyla)

# Create needed export directories
dir.create("Graphs", showWarnings = F)
dir.create("Figures", showWarnings = F)

# We also want to customize some plotting aesthetics for our ggplot plots that we can do here
all.cols <- c("Larval paunch" = "#a50026", "Larval ileum" = "#f46d43",
              "Larval midgut" = "#fee090",
              "Adult midgut" = "#abd9e9", "Adult hindgut" = "#4575b4")
sex.shps <- c("male" = 24, "female" = 25, "larva"= 21)
dodge <- position_dodge(width = 0.5)
pref_theme <- theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'none',
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        legend.background = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1))

## -------------------------------------------------------------- ##
          # Figure 1 - Alpha Diversity Superfigure ####
## -------------------------------------------------------------- ##
# Summarize the response variables
alpha_plotdf <- alpha %>%
  group_by(Stage.Gut) %>%
  dplyr::summarise(
    shan = mean(Shannon, na.rm = T),
    shan_se = ( sd(Shannon, na.rm = T) / sqrt(dplyr::n()) ),
    simp = mean(Simpson, na.rm = T),
    simp_se = ( sd(Simpson, na.rm = T) / sqrt(dplyr::n()) ),
    ace = mean(ACE, na.rm = T),
    ace_se = ( sd(ACE, na.rm = T) / sqrt(dplyr::n()) ),
    chao = mean(Chao1, na.rm = T),
    chao_se = ( sd(Chao1, na.rm = T) / sqrt(dplyr::n()) )) %>%
  # Re-level stage.gut factor order
  dplyr::mutate(Stage.Gut = factor(Stage.Gut,
                                   levels = c("Adult hindgut", "Adult midgut", "Larval midgut",
                                              "Larval ileum", "Larval paunch")))

# Shannon diversity bar graph
(shan.plt <- ggplot(alpha_plotdf, aes(y = shan, x = Stage.Gut,
                         fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = shan + shan_se, ymin = shan - shan_se), width = 0.2) +
  geom_text(label = "b", x = 0.7, y =  3.8, size = 5) +
  geom_text(label = "a", x = 1.7, y =  2.3, size = 5) +
  geom_text(label = "c", x = 2.7, y =  6.15, size = 5) +
  geom_text(label = "bc", x = 3.7, y = 4.75, size = 5) +
  geom_text(label = "c", x = 4.7, y =  5.45, size = 5) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "Shannon Diversity") +
  pref_theme + theme(legend.position = 'none'))

# Simpson Diversity
(simp.plt <- ggplot(alpha_plotdf, aes(y = simp, x = Stage.Gut,
                         fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = simp + simp_se, ymin = simp - simp_se), width = 0.2) +
  geom_text(label = "b", x = 0.7, y =  0.99, size = 5) +
  geom_text(label = "a", x = 1.7, y =  0.82, size = 5) +
  geom_text(label = "b", x = 2.7, y =  1.05, size = 5) +
  geom_text(label = "b", x = 3.7, y =  1, size = 5) +
  geom_text(label = "b", x = 4.7, y =  1.025, size = 5) +
  ylim(0, 1.1) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "Simpson Diversity") +
  pref_theme + theme(legend.position = 'none'))

# ACE Diversity
(ace.plt <- ggplot(alpha_plotdf, aes(y = ace, x = Stage.Gut,
                         fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = ace + ace_se, ymin = ace - ace_se), width = 0.2) +
  geom_text(label = "ab", x = 0.7, y = 180, size = 5) +
  geom_text(label = "a", x = 1.7, y =  100, size = 5) +
  geom_text(label = "bc", x = 2.7, y = 750, size = 5) +
  geom_text(label = "b", x = 3.7, y =  600, size = 5) +
  geom_text(label = "c", x = 4.7, y =  985, size = 5) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "ACE") +
  pref_theme + theme(legend.position = 'none'))

# Chao 1 Diversity
(chao.plt <- ggplot(alpha_plotdf, aes(y = chao, x = Stage.Gut,
                         fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = chao + chao_se, ymin = chao - chao_se), width = 0.2) +
  geom_text(label = "ab", x = 0.7, y = 180, size = 5) +
  geom_text(label = "a", x = 1.7, y =  100, size = 5) +
  geom_text(label = "bc", x = 2.7, y = 750, size = 5) +
  geom_text(label = "b", x = 3.7, y =  600, size = 5) +
  geom_text(label = "c", x = 4.7, y =  985, size = 5) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "Chao1 Index") +
  pref_theme + theme(legend.position = 'none'))

# Create the superfigure
cowplot::plot_grid(chao.plt, ace.plt,  simp.plt, shan.plt,
                   align = 'v', rel_heights = c(0.7, 1),
                   ncol = 2, nrow = 2, labels = "AUTO")

# Save it
ggplot2::ggsave(file.path("Figures", "Alpha-Diversity-Superfigure.tiff"),
                device = 'tiff',
                width = 7, height = 7,
                plot = last_plot())

## -------------------------------------------------------------- ##
            # Figure 2 - Beta Diversity Superfigure ####
## -------------------------------------------------------------- ##
# Re-order factor levels for PCoA
unique(beta$Stage.Gut)
beta$Stage.Gut <- factor(beta$Stage.Gut,
                         levels = c("Larval paunch", "Larval ileum", "Larval midgut",
                                    "Adult midgut", "Adult hindgut"))
unique(beta$Stage.Gut)

# Strip out the ID missing from the unifrac data
frc.data <- beta %>%
  filter(Sample.ID != "Amid48") %>%
  as.data.frame()

# Calculate/read in distances
bc.dist <- vegdist(beta[-c(1:5)], method = 'bray')
jac.dist <- vegdist(beta[-c(1:5)], method = 'jaccard')
wtd.frc.dist <- read.csv(file.path("Data", "wtd-unfrc-dist.csv"))[-1]
uwt.frc.dist <- read.csv(file.path("Data", "unwtd-unfrc-dist.csv"))[-1]

# Run a PCoA on all distance matrices
bc.pnts <- ape::pcoa(bc.dist)
jacc.pnts <- ape::pcoa(jac.dist)
wtd.frc.pnts <- ape::pcoa(wtd.frc.dist)
uwt.frc.pnts <- ape::pcoa(uwt.frc.dist)

# No unique ordination combinations are used as figures
## Only a single ordination gets to be a figure in the main text of the paper

## Bray Curtis
supportR::pcoa_ord(mod = bc.pnts, groupcol = beta$Stage.Gut,
                   leg_cont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
                   colors = all.cols, leg_pos = "bottomleft")

## Jaccard
supportR::pcoa_ord(mod = jacc.pnts, groupcol = beta$Stage.Gut,
                   leg_cont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
                   colors = all.cols, leg_pos = "bottomleft")

## Weighted Unifrac
supportR::pcoa_ord(mod = wtd.frc.pnts, groupcol = frc.data$Stage.Gut,
                   leg_cont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
                   colors = all.cols, leg_pos = "topright")

## Unweighted Unifrac
supportR::pcoa_ord(mod = uwt.frc.pnts, groupcol = frc.data$Stage.Gut,
                   leg_cont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
                   colors = all.cols, leg_pos = "topleft")

## -------------------------------------------------------------- ##
                # Relative Abundance - Global ####
## -------------------------------------------------------------- ##
# Look at the phylum- & family-level abundance
head(phyla)
head(fams)

# What is the threshold of relative abundance we want?
threshold <- 0.05

# Identify global abundance totals
phy.global.abun <- sum(phyla$Abundance, na.rm = T)
fam.global.abun <- sum(fams$Abundance, na.rm = T)

# Now identify abundance of over/under threshold percent
phy.global.thresh <- phy.global.abun * threshold
fam.global.thresh <- fam.global.abun * threshold

# Slim down the phyla dataframe
phy.global.df <- phyla %>%
  # Remove un-needed columns
  dplyr::select(Domain, Phylum, Abundance) %>%
  # Group by phylum
  group_by(Phylum) %>%
  # Sum to get one instance of each
  dplyr::summarise(Abundance = sum(Abundance, na.rm = T)) %>%
  # Identify which are above/below threshold
  dplyr::mutate(Abundant_Phyla = dplyr::case_when(
    Abundance > phy.global.thresh ~ Phylum,
    TRUE ~ "Phyla < 5% Total")) %>%
  # Sum through any taxa under the threshold value
  group_by(Abundant_Phyla) %>%
  dplyr::summarise(Abundance = sum(Abundance, na.rm = T)) %>%
  # Rename the phyla column more simply
  dplyr::rename(Phylum = Abundant_Phyla) %>%
  # Calculate relativeAbundance
  dplyr::mutate(relativeabun = (Abundance / phy.global.abun) * 100) %>%
  # Return a dataframe
  as.data.frame()

# Do the same for family-level abundance
fam.global.df <- fams %>%
  dplyr::select(Domain, Family, Abundance) %>%
  group_by(Family) %>%
  dplyr::summarise(Abundance = sum(Abundance, na.rm = T)) %>%
  dplyr::mutate(Abundant_Families = dplyr::case_when(
    Abundance > fam.global.thresh ~ Family,
    TRUE ~ "Families < 5% Total")) %>%
  group_by(Abundant_Families) %>%
  dplyr::summarise(Abundance = sum(Abundance, na.rm = T)) %>%
  dplyr::rename(Family = Abundant_Families) %>%
  dplyr::mutate(relativeabun = (Abundance / fam.global.abun) * 100) %>%
  as.data.frame()

# Look at what that produced
head(phy.global.df)
head(fam.global.df)

# Get a custom color vector
phy.global.colors <- c("Phyla < 5% Total" = "white",# Phyla = blues
                   "Bacteroidetes" = "#023858", "Firmicutes" = "#3690c0",
                   "Proteobacteria" = "#a6bddb")
fam.global.colors <- c("Families < 5% Total" = "#ffffe5", # Family = oranges
                       "Desulfovibrionaceae" = "#662506", "Dysgonomonadaceae" = "#cc4c02",
                       "Lachnospiraceae" = "#ec7014", "Rikenellaceae" = "#fec44f",
                       "Ruminococcaceae" = "#fee391")

# Graph phylum level abundance
(phy.global.plt <- ggplot(phy.global.df, aes(x = 'x', y = relativeabun,
                           fill = reorder(Phylum, relativeabun),
                           color = 'x')) +
  geom_bar(stat = 'identity', position = 'stack') +
  scale_color_manual(values = 'black') +
  labs(x = 'Phyla', y = "Relative Abundance (%)") +
  scale_fill_manual(values = phy.global.colors) +
  # Remaining aesthetics things
  pref_theme + guides(color = 'none') +
  theme(legend.position = 'right',
        axis.text.x = element_blank(),
        legend.key = element_rect(color = "black")))

# Do the same for family
(fam.global.plt <- ggplot(fam.global.df, aes(x = 'x', y = relativeabun,
                                             fill = reorder(Family, relativeabun),
                                             color = 'x')) +
    geom_bar(stat = 'identity', position = 'stack') +
    scale_color_manual(values = 'black') +
    labs(x = 'Families', y = "Relative Abundance (%)") +
    scale_fill_manual(values = fam.global.colors) +
    # Remaining aesthetics things
    pref_theme + guides(color = 'none') +
    theme(legend.position = 'right',
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          legend.key = element_rect(color = "black")))

# Create combination graph
cowplot::plot_grid(phy.global.plt, fam.global.plt,
                   ncol = 2, nrow = 1, labels = "AUTO")

# Save it
ggplot2::ggsave(file = file.path("Figures", "Abundance-Superfigure.tiff"),
                device = 'tiff', width = 6, height = 5,
                plot = last_plot())

## -------------------------------------------------------------- ##
                # Figure 4 - Relative Abundance ####
## -------------------------------------------------------------- ##
## -------------------------------------------- ##
              # F4 Housekeeping ####
## -------------------------------------------- ##
# Look at the range of relative abundances
  ## Visually
psych::multi.hist(fams$relativeAbun)
psych::multi.hist(phyla$relativeAbun)
  ## Numerically
range(fams$relativeAbun)
range(phyla$relativeAbun)

# Identify the threshold value of "too low"
abun.thresh <- 5

# We have a good alternative to dropping low abundances altogether:
fams.v2 <- fams %>%
  # Remove unneeded columns
  dplyr::select(Sample.ID, Family, relativeAbun) %>%
  dplyr::mutate(
    # Change the family to "other" if it's less than the cutoff
    Family = ifelse(relativeAbun < abun.thresh,
                    yes = paste0("Family < ", abun.thresh, "% Total"),
                    no = Family)) %>%
  # Then, within sample sum through the relative abundances
  group_by(Sample.ID, Family) %>%
  dplyr::summarise(relativeAbun = sum(relativeAbun, na.rm = T)) %>%
  as.data.frame()

# Do the same for phylum
phyla.v2 <- phyla %>%
  dplyr::select(Sample.ID, Phylum, relativeAbun) %>%
  dplyr::mutate(Phylum = ifelse(relativeAbun < abun.thresh,
                                yes = paste0("Phyla < ", abun.thresh, "% Total"),
                                no = Phylum)) %>%
  group_by(Sample.ID, Phylum) %>%
  dplyr::summarise(relativeAbun = sum(relativeAbun, na.rm = T)) %>%
  as.data.frame()

# Check the ranges to see that worked
range(fams.v2$relativeAbun)
range(phyla.v2$relativeAbun)

## -------------------------------------------- ##
               # F4 Graphing ####
## -------------------------------------------- ##
# Family-level *relative* abundance figure
ggplot(fams.v2, aes(x = Sample.ID, y = relativeAbun,
                    fill = Family, color = 'x',
                    order = -relativeAbun)) +
  geom_bar(stat = 'identity', position = 'stack', width = 0.65) +
  scale_color_manual(values = 'black') +
  labs(x = 'Sample ID', y = "Relative Abundance (%)") +
  #scale_fill_manual(values = sp::bpy.colors()) +
  pref_theme + theme(axis.title.x = element_blank())
  #theme(legend.position = 'right', axis.text.x = element_text(size = 8))

# Re-level the Sample.ID factor
phyla.v2$Sample.ID <- factor(phyla.v2$Sample.ID,
                             levels = c(
                               # Adult midgut
                               "Amid7", "Amid12", "Amid39", "Amid40",
                               "Amid41", "Amid46", "Amid48", "Amid49",
                               # Adult hindgut
                               "hind4", "hind7", "hind12", "hind39",
                               "hind40", "hind41", "hind46", "hind48",
                               "hind49",
                               # Larval midgut
                               "Lmid2", "Lmid3", "Lmid5", "Lmid9",
                               # Larval ileum
                               "Ileum2", "Ileum3", "Ileum5", "Ileum7", "Ileum9",
                               "Ileum12", "Ileum17", "Ileum18", "Ileum19",
                               "Ileum20", "Ileum21",
                               # Larval paunch
                               "Paunch2", "Paunch3", "Paunch5", "Paunch9",
                               "Paunch12", "Paunch17", "Paunch18", "Paunch19",
                               "Paunch20", "Paunch21"))

# Set colors for each phylum
phyla.cols <- c("white", "#ffff99", "#6a3d9a",
                "#cab2d6","#ff7f00","#fdbf6f",
                "#e31a1c", "#fb9a99", "#33a02c",
                "#b2df8a", "#1f78b4", "#a6cee3")

# Phylum-level relative abundance
ggplot(phyla.v2, aes(x = Sample.ID, y = relativeAbun,
                     fill = reorder(Phylum, relativeAbun),
                     color = 'x')) +
  geom_bar(stat = 'identity', position = 'stack') +
  scale_color_manual(values = 'black') +
  labs(x = 'Sample ID', y = "Relative Abundance (%)") +
  scale_fill_manual(values = phyla.cols) +
  # Add vertical lines separating gut/life stages
  geom_vline(xintercept = 8.5) +
  geom_vline(xintercept = 17.5) +
  geom_vline(xintercept = 21.5) +
  geom_vline(xintercept = 32.5) +
  # Add text to top of graph defining each chunk (delinated by vertical lines)
  geom_text(label = "Adult Midgut", hjust = 'center', x = 4.25, y = 103) +
  geom_text(label = "Adult Hindgut", hjust = 'center', x = 13, y = 103) +
  geom_text(label = "Larval Mid.", hjust = 'center', x = 19.5, y = 103) +
  geom_text(label = "Ileum", hjust = 'center', x = 27, y = 103) +
  geom_text(label = "Paunch", hjust = 'center', x = 37.5, y = 103) +
  # Remaining aesthetics things
  pref_theme + guides(color = 'none') +
  theme(legend.position = 'right',
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 9),
        legend.key = element_rect(color = "black"))

# Save this!
ggplot2::ggsave(file.path("Figures", "Relative-Abundance-Superfigure.tiff"),
                width = 11, height = 5, plot = last_plot())

## -------------------------------------------------------------- ##
   # Figure 5 - Within Life/Gut Phyla Abundance Superfigure ####
## -------------------------------------------------------------- ##
## -------------------------------------------- ##
             # F5 Housekeeping ####
## -------------------------------------------- ##
# Want to exclude very low abundance phyla to clarify the figures

# Re-identify abundance threshold (%)
new.abun.thresh <- 5

# Create subsets with relative abundance threshold implemented
amid.phyl <- phyla %>%
  dplyr::mutate(
    Phylum = ifelse(relativeAbun < new.abun.thresh,
                    yes = paste0("Phyla < ", new.abun.thresh, "% Total"),
                    no = Phylum) ) %>%
  filter(Stage.Gut == "Adult midgut") %>%
  group_by(Stage.Gut, Phylum) %>%
  dplyr::summarise(
    abun = mean(Abundance, na.rm = T),
    se = sd(Abundance, na.rm = T) / sqrt(dplyr::n())
  ) %>%
  dplyr::rename(Abundance = abun) %>%
  as.data.frame()

ahind.phyl <- phyla %>%
  dplyr::mutate(
    Phylum = ifelse(relativeAbun < new.abun.thresh,
                    yes = paste0("Phyla < ", new.abun.thresh, "% Total"),
                    no = Phylum) ) %>%
  filter(Stage.Gut == "Adult hindgut") %>%
  group_by(Stage.Gut, Phylum) %>%
  dplyr::summarise(
    abun = mean(Abundance, na.rm = T),
    se = sd(Abundance, na.rm = T) / sqrt(dplyr::n())
  ) %>%
  dplyr::rename(Abundance = abun) %>%
  as.data.frame()

lmid.phyl <- phyla %>%
  dplyr::mutate(
    Phylum = ifelse(relativeAbun < new.abun.thresh,
                    yes = paste0("Phyla < ", new.abun.thresh, "% Total"),
                    no = Phylum) ) %>%
  filter(Stage.Gut == "Larval midgut") %>%
  group_by(Stage.Gut, Phylum) %>%
  dplyr::summarise(
    abun = mean(Abundance, na.rm = T),
    se = sd(Abundance, na.rm = T) / sqrt(dplyr::n())
  ) %>%
  dplyr::rename(Abundance = abun) %>%
  as.data.frame()

ileum.phyl <- phyla %>%
  dplyr::mutate(
    Phylum = ifelse(relativeAbun < new.abun.thresh,
                    yes = paste0("Phyla < ", new.abun.thresh, "% Total"),
                    no = Phylum) ) %>%
  filter(Stage.Gut == "Larval ileum") %>%
  group_by(Stage.Gut, Phylum) %>%
  dplyr::summarise(
    abun = mean(Abundance, na.rm = T),
    se = sd(Abundance, na.rm = T) / sqrt(dplyr::n())
  ) %>%
  dplyr::rename(Abundance = abun) %>%
  as.data.frame()

paunch.phyl <- phyla %>%
  dplyr::mutate(
    Phylum = ifelse(relativeAbun < new.abun.thresh,
                    yes = paste0("Phyla < ", new.abun.thresh, "% Total"),
                    no = Phylum) ) %>%
  filter(Stage.Gut == "Larval paunch") %>%
  group_by(Stage.Gut, Phylum) %>%
  dplyr::summarise(
    abun = mean(Abundance, na.rm = T),
    se = sd(Abundance, na.rm = T) / sqrt(dplyr::n())
  ) %>%
  dplyr::rename(Abundance = abun) %>%
  as.data.frame()

## -------------------------------------------- ##
     # Fig 5 Graphing - Low Abun Excl. ####
## -------------------------------------------- ##
# Make component graphs of phylum-level abundance within each life stage/gut bit
## Adult midgut
(amid <- ggplot(amid.phyl, aes(y = Abundance, x = reorder(Phylum, -Abundance),
                              fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  geom_text(label = "Adult Midgut", hjust = 'center',
            x = 3, y = 1475) +
  labs(x = "Phylum", y = "Abundance") +
  pref_theme + theme(legend.position = 'none',
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 9),
                     axis.title.y = element_text(size = 12),
                     axis.title.x = element_blank()) )

## Adult hindgut
(ahind <- ggplot(ahind.phyl, aes(y = Abundance, x = reorder(Phylum, -Abundance),
                                fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  geom_text(label = "Adult Hindgut", hjust = 'center',
            x = 3, y = 8500) +
  labs(x = "Phylum", y = "Abundance") +
  pref_theme + theme(legend.position = 'none',
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 9),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank()) )

## Larval midgut
(lmid <- ggplot(lmid.phyl, aes(y = Abundance, x = reorder(Phylum, -Abundance),
                              fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  geom_text(label = "Larval Midgut", hjust = 'center',
            x = 3.5, y = 1750) +
  labs(x = "Phyla", y = "Abundance") +
  pref_theme + theme(legend.position = 'none',
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 9),
                     axis.title = element_text(size = 12)) )

## Larval ileum
(lile <- ggplot(ileum.phyl, aes(y = Abundance, x = reorder(Phylum, -Abundance),
                               fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  geom_text(label = "Larval Ileum", hjust = 'center',
            x = 5, y = 4150) +
  ylim(c(0, 4250)) +
  labs(x = "Phyla", y = "Abundance") +
  pref_theme + theme(legend.position = 'none',
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 9),
                     axis.title = element_text(size = 12),
                     axis.title.y = element_blank()) )

## Larval paunch
(lpau <- ggplot(paunch.phyl, aes(y = Abundance, x = reorder(Phylum, -Abundance),
                                fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  geom_text(label = "Larval Paunch", hjust = 'center',
            x = 3, y = 13000) +
  labs(x = "Phyla", y = "Abundance") +
  pref_theme + theme(legend.position = 'none',
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 9),
                     axis.title = element_text(size = 12),
                     axis.title.y = element_blank()) )

# Assemble into larger composite graph
plot_grid(amid, ahind, NA, lmid, lile, lpau,
          ncol = 3, nrow = 2,
          labels = c("A", "B", NA, "C", "D", "E"))

# Option A:
ggplot2::ggsave(file.path("Figures", "Relative-Abundance-Subsets-Superfigure.tiff"),
                width = 11, height = 6.5, plot = last_plot())

## -------------------------------------------------------------- ##
                   # Figure 6 - Venn Diagram ####
## -------------------------------------------------------------- ##
# Look at beta diversity data
str(beta[1:20])

# Prepare it for use in the venn diagram
venn.df <- beta %>%
  # Pivot to long format
  pivot_longer(-Sample.ID:-Sex,
               names_to = 'Taxon',
               values_to = 'Abundance') %>%
  # Group by Stage.Gut and Taxon
  dplyr::group_by(Stage.Gut, Taxon) %>%
  # Sum across samples within the aforementioned groups
  dplyr::summarise(Abundance = sum(Abundance, na.rm = T)) %>%
  # Do some changes to work with the `ggvenn` function/package
  dplyr::mutate(
    Abundance = dplyr::case_when(
      is.na(Abundance) ~ FALSE,
      Abundance == 0 ~ FALSE,
      TRUE ~ TRUE)) %>%
  # Pivot to wide format with Stage.Gut as columns
  pivot_wider(id_cols = Taxon,
              names_from = Stage.Gut,
              values_from = Abundance )

# Check it out again
str(venn.df)

# Make Venn diagram for each and save them separately
## Adults
(adult_ggvenn <- ggvenn(data = venn.df,
       columns = c("Adult hindgut", "Adult midgut"),
       fill_color = c('blue', 'light blue'),
       text_size = 6, set_name_size = 4.5,
       show_percentage = F) )
ggplot2::ggsave(file.path("Graphs", "Venn-Diagram-ADULT.tiff"),
                width = 4, height = 4, plot = last_plot())

## Larvae
(larva_ggvenn <- ggvenn(data = venn.df,
        columns = c("Larval paunch", "Larval ileum", "Larval midgut"),
        fill_color = c('red', 'orange', 'yellow'),
        text_size = 4, set_name_size = 4.5,
        show_percentage = F) )
ggplot2::ggsave(file.path("Graphs", "Venn-Diagram-LARVAE.tiff"),
                width = 4, height = 4,  plot = last_plot())

# Make and save a composite figure
plot_grid(adult_ggvenn, larva_ggvenn,
          ncol = 1, nrow = 2, labels = c("A", "B"))

# Option A:
ggplot2::ggsave(file.path("Figures", "Venn-Diagram-Superfigure.tiff"),
                width = 5.5, height = 6, plot = last_plot())

# END ####
