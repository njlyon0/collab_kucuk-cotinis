## -------------------------------------------------------------- ##
           # Kucuk Cotinis Project - Visualization Code
## -------------------------------------------------------------- ##
# Code written by Nicholas J Lyon

# PURPOSE:
  ## Perform data visualization of the gut microbe communities found in Cotinis nitida

# Clear environment (always better to start with tabula rasa)
rm(list = ls())

# Necessary libraries
# install.packages("librarian")
librarian::shelf(tidyverse, vegan, ape, supportR)

## -------------------------------------------- ##
      # Data Retrieval & Housekeeping ####
## -------------------------------------------- ##
# Retrieve the relevant datasets
alpha <- read.csv(file.path("Data", "Tidy Data", "alpha-diversity-data.csv"))
beta <- read.csv(file.path("Data", "Tidy Data", "beta-diversity-data.csv"))
fams <- read.csv(file.path("Data", "Tidy Data", "family-abun.csv"))
phyla <- read.csv(file.path("Data", "Tidy Data", "phylum-abun.csv"))
wtd.frc.dist <- read.csv(file.path("Data", "wtd-unfrc-dist.csv")) %>%
  dplyr::select(-X)
uwt.frc.dist <- read.csv(file.path("Data", "unwtd-unfrc-dist.csv")) %>%
  dplyr::select(-X)

# Look at them to be sure nothing obvious is wrong
str(alpha)
str(beta)
str(fams)
str(phyla)
str(wtd.frc.dist)

# Create needed export directories
dir.create("Graphs", showWarnings = F)
dir.create("Figures", showWarnings = F)

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

# We also want to customize some plotting aesthetics for our ggplot plots that we can do here
all.cols <- c("Larval paunch" = "#a50026", "Larval ileum" = "#f46d43",
              "Larval midgut" = "#fee090",
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
supportR::pcoa_ord(mod = bc.pnts, groupcol = beta$Stage.Gut,
                leg_cont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
                colors = all.cols, leg_pos = "bottomleft")

# Once that looks good save it out using this dev.off bit
tiff(file = file.path("Graphs", "PCoA-bray-curtis.tiff"))
supportR::pcoa_ord(mod = bc.pnts, groupcol = beta$Stage.Gut,
                leg_cont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
                colors = all.cols, leg_pos = "bottomleft")
dev.off()

# Now make the Jaccard PCoA
supportR::pcoa_ord(mod = jacc.pnts, groupcol = beta$Stage.Gut,
                leg_cont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
                colors = all.cols, leg_pos = "bottomleft")

# Once that looks good save it out using this dev.off bit
tiff(file = file.path("Graphs", "PCoA-jaccard.tiff"))
supportR::pcoa_ord(mod = jacc.pnts, groupcol = beta$Stage.Gut,
                leg_cont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
                colors = all.cols, leg_pos = "bottomleft")
dev.off()

# Now make the Weighted Unifrac PCoA
supportR::pcoa_ord(mod = wtd.frc.pnts, groupcol = frc.data$Stage.Gut,
                leg_cont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
                colors = all.cols, leg_pos = "topright")

# Save the ordination
tiff(file = file.path("Graphs", "PCoA-weighted-unifrac.tiff"))
supportR::pcoa_ord(mod = wtd.frc.pnts, groupcol = frc.data$Stage.Gut,
                leg_cont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
                colors = all.cols, leg_pos = "topright")
dev.off()

# Finally, do the Unweighted Unifrac PCoA
supportR::pcoa_ord(mod = uwt.frc.pnts, groupcol = frc.data$Stage.Gut,
                leg_cont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
                colors = all.cols, leg_pos = "topleft")

# Save the ordination
tiff(file = file.path("Graphs", "PCoA-unweighted-unifrac.tiff"))
supportR::pcoa_ord(mod = uwt.frc.pnts, groupcol = frc.data$Stage.Gut,
                leg_cont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
                colors = all.cols, leg_pos = "topleft")
dev.off()

## -------------------------------------------- ##
        # Alpha Diversity Bar Graphs ####
## -------------------------------------------- ##
# Incl. Chao, ACE, Simpson, and Shannon diversity & Pielou's Evenness
  ## (Diversity ~ Life Stage + Gut Section, data = all)

# Summarize the response variables
alpha_plotdf <- alpha %>%
  group_by(Stage.Gut) %>%
  dplyr::summarise(
    shan = mean(Shannon, na.rm = T),
    shan_se = ( sd(Shannon, na.rm = T) / sqrt(dplyr::n()) ),
    piel = mean(Pielous, na.rm = T),
    piel_se = ( sd(Pielous, na.rm = T) / sqrt(dplyr::n()) ),
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
ggplot(alpha_plotdf, aes(y = shan, x = Stage.Gut,
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
  pref_theme + theme(legend.position = 'none')
ggsave(plot = last_plot(), file.path("Graphs", "bar-shannon.tiff"),
       width = 5, height = 5, units = 'in')

# Pielou's evenness bar graph
ggplot(alpha_plotdf, aes(y = piel, x = Stage.Gut,
                         fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = piel + piel_se, ymin = piel - piel_se), width = 0.2) +
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
ggsave(plot = last_plot(), file.path("Graphs", "bar-pielou.tiff"),
       width = 5, height = 5, units = 'in')

# Simpson Diversity
ggplot(alpha_plotdf, aes(y = simp, x = Stage.Gut,
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
  pref_theme + theme(legend.position = 'none')
ggsave(plot = last_plot(), file.path("Graphs", "bar-simpson.tiff"),
       width = 5, height = 5, units = 'in')

# ACE Diversity
ggplot(alpha_plotdf, aes(y = ace, x = Stage.Gut,
                          fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = ace + ace_se, ymin = ace - ace_se), width = 0.2) +
  geom_text(label = "ab", x = 0.7, y = 180, size = 5) +
  geom_text(label = "a", x = 1.7, y =  100, size = 5) +
  geom_text(label = "bc", x = 2.7, y = 740, size = 5) +
  geom_text(label = "b", x = 3.7, y =  590, size = 5) +
  geom_text(label = "c", x = 4.7, y =  975, size = 5) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "ACE") +
  pref_theme + theme(legend.position = 'none')
ggsave(plot = last_plot(), file.path("Graphs", "bar-ace.tiff"),
       width = 5, height = 5, units = 'in')

# Chao 1 Diversity
ggplot(alpha_plotdf, aes(y = chao, x = Stage.Gut,
                         fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = chao + chao_se, ymin = chao - chao_se), width = 0.2) +
  geom_text(label = "ab", x = 0.7, y = 180, size = 5) +
  geom_text(label = "a", x = 1.7, y =  100, size = 5) +
  geom_text(label = "bc", x = 2.7, y = 740, size = 5) +
  geom_text(label = "b", x = 3.7, y =  590, size = 5) +
  geom_text(label = "c", x = 4.7, y =  975, size = 5) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "Chao1 Index") +
  pref_theme + theme(legend.position = 'none')
ggsave(plot = last_plot(), file.path("Graphs", "bar-chao.tiff"),
       width = 5, height = 5, units = 'in')

## -------------------------------------------- ##
        # Family Abundance Bar Graphs ####
## -------------------------------------------- ##
# Summarize data
fam_plotdf <- fams %>%
  group_by(Stage.Gut, Family) %>%
  dplyr::summarise(
    abun = mean(Abundance, na.rm = T),
    se = ( sd(Abundance, na.rm = T) / sqrt(dplyr::n()) )) %>%
  dplyr::rename(Abundance = abun)

# Make and save each graph
## Adult midgut
ggplot(data = filter(fam_plotdf, Stage.Gut == "Adult midgut"),
       aes(y = Abundance, x = reorder(Family, -Abundance),
           fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Family", y = "Adult Midgut Abundance") +
  pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))
ggsave(plot = last_plot(), file.path("Graphs", "fam-abun-adult-midgut.tiff"),
       width = 6, height = 6, units = 'in')

## Adult hindgut
ggplot(data = filter(fam_plotdf, Stage.Gut == "Adult hindgut"),
       aes(y = Abundance, x = reorder(Family, -Abundance),
           fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Family", y = "Adult Hindgut Abundance") +
  pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))
ggsave(plot = last_plot(), file.path("Graphs", "fam-abun-adult-hindgut.tiff"),
       width = 6, height = 6, units = 'in')

## Larval midgut
ggplot(data = filter(fam_plotdf, Stage.Gut == "Larval midgut"),
       aes(y = Abundance, x = reorder(Family, -Abundance),
           fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Family", y = "Larval Midgut Abundance") +
  pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))
ggsave(plot = last_plot(), file.path("Graphs", "fam-abun-larval-midgut.tiff"),
       width = 6, height = 6, units = 'in')

## Larval paunch
ggplot(data = filter(fam_plotdf, Stage.Gut == "Larval paunch"),
       aes(y = Abundance, x = reorder(Family, -Abundance),
           fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Family", y = "Larval Paunch Abundance") +
  pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))
ggsave(plot = last_plot(), file.path("Graphs", "fam-abun-larval-paunch.tiff"),
       width = 6, height = 6, units = 'in')

## Larval ileum
ggplot(data = filter(fam_plotdf, Stage.Gut == "Larval ileum"),
       aes(y = Abundance, x = reorder(Family, -Abundance),
           fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Family", y = "Larval Ileum Abundance") +
  pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))
ggsave(plot = last_plot(), file.path("Graphs", "fam-abun-larval-ileum.tiff"),
       width = 6, height = 6, units = 'in')

## -------------------------------------------- ##
        # Phylum Abundance Bar Graphs ####
## -------------------------------------------- ##
# Summarize data
phyla_plotdf <- phyla %>%
  group_by(Stage.Gut, Phylum) %>%
  dplyr::summarise(
    abun = mean(Abundance, na.rm = T),
    se = ( sd(Abundance, na.rm = T) / sqrt(dplyr::n()) )) %>%
  dplyr::rename(Abundance = abun)

# Make and save each graph
## Adult midgut
ggplot(data = filter(phyla_plotdf, Stage.Gut == "Adult midgut"),
       aes(y = Abundance, x = reorder(Phylum, -Abundance),
           fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Phylum", y = "Adult Midgut Abundance") +
  pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))
ggsave(plot = last_plot(), file.path("Graphs", "phyla-abun-adult-midgut.tiff"),
       width = 6, height = 6, units = 'in')

## Adult hindgut
ggplot(data = filter(phyla_plotdf, Stage.Gut == "Adult hindgut"),
       aes(y = Abundance, x = reorder(Phylum, -Abundance),
           fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Phylum", y = "Adult Hindgut Abundance") +
  pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))
ggsave(plot = last_plot(), file.path("Graphs", "phyla-abun-adult-hindgut.tiff"),
       width = 6, height = 6, units = 'in')

## Larval midgut
ggplot(data = filter(phyla_plotdf, Stage.Gut == "Larval midgut"),
       aes(y = Abundance, x = reorder(Phylum, -Abundance),
           fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Phylum", y = "Larval Midgut Abundance") +
  pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))
ggsave(plot = last_plot(), file.path("Graphs", "phyla-abun-larval-midgut.tiff"),
       width = 6, height = 6, units = 'in')

## Larval paunch
ggplot(data = filter(phyla_plotdf, Stage.Gut == "Larval paunch"),
       aes(y = Abundance, x = reorder(Phylum, -Abundance),
           fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Phylum", y = "Larval Paunch Abundance") +
  pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))
ggsave(plot = last_plot(), file.path("Graphs", "phyla-abun-larval-paunch.tiff"),
       width = 6, height = 6, units = 'in')

## Larval ileum
ggplot(data = filter(phyla_plotdf, Stage.Gut == "Larval ileum"),
       aes(y = Abundance, x = reorder(Phylum, -Abundance),
           fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance - se), width = 0.2) +
  scale_fill_manual(values = all.cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Phylum", y = "Larval Ileum Abundance") +
  pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))
ggsave(plot = last_plot(), file.path("Graphs", "phyla-abun-larval-ileum.tiff"),
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
  dplyr::summarise(Abundance = sum(Abundance, na.rm = T)) %>%
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
ggsave(plot = last_plot(), file.path("Graphs", "phyla-abun-all.tiff"),
       width = 6, height = 6, units = 'in')

# END ####

