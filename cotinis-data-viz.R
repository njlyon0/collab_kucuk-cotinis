## -------------------------------------------------------------- ##
           # Kucuk Cotinis Project - Visualization Code
## -------------------------------------------------------------- ##
# Code written by Nicholas J Lyon

# PURPOSE:
  ## Perform data visualization of the gut microbe communities found in Cotinis nitida

# Necessary libraries
# install.packages("librarian")
librarian::shelf(tidyverse, vegan, ape, supportR)

# Clear environment (always better to start with tabula rasa)
rm(list = ls())

# Create needed export directory
dir.create(path = file.path("graphs"), showWarnings = F)

## -------------------------------------------- ##
      # Data Retrieval & Housekeeping ####
## -------------------------------------------- ##
# Retrieve the relevant datasets
alpha <- read.csv(file.path("data", "tidy_data", "alpha-diversity-data.csv"))
beta <- read.csv(file.path("data", "tidy_data", "beta-diversity-data.csv"))
fams <- read.csv(file.path("data", "tidy_data", "family-abun.csv"))
phyla <- read.csv(file.path("data", "tidy_data", "phylum-abun.csv"))
wtd_frc_dist <- read.csv(file.path("data", "wtd-unfrc-dist.csv")) %>%
  dplyr::select(-X)
uwt_frc_dist <- read.csv(file.path("data", "unwtd-unfrc-dist.csv")) %>%
  dplyr::select(-X)

# Look at them to be sure nothing obvious is wrong
dplyr::glimpse(alpha)
dplyr::glimpse(beta)
dplyr::glimpse(fams)
dplyr::glimpse(phyla)
dplyr::glimpse(wtd_frc_dist)

# Calculate the Bray Curtis and Jaccard distances
bc_dist <- vegdist(beta[-c(1:5)], method = 'bray')
jac_dist <- vegdist(beta[-c(1:5)], method = 'jaccard')

# We also want to order the factor levels of what we're using as groups for PCoAs
unique(beta$Stage.Gut)
beta$Stage.Gut <- factor(beta$Stage.Gut,
                              levels = c("Larval paunch", "Larval ileum", "Larval midgut",
                                         "Adult midgut", "Adult hindgut"))
unique(beta$Stage.Gut)

# Strip out the ID missing from the unifrac data
frc_data <- dplyr::filter(beta, Sample.ID != "Amid48")

# Roy (Kucuk) also wants some plots/ordinations that include larvae only
alpha_larv <- dplyr::filter(alpha, Lifestage == "larva")

# We also want to customize some plotting aesthetics for our ggplot plots that we can do here
all_cols <- c("Larval paunch" = "#a50026", "Larval ileum" = "#f46d43",
              "Larval midgut" = "#fee090",
              "Adult midgut" = "#abd9e9", "Adult hindgut" = "#4575b4")
sex_shps <- c("male" = 24, "female" = 25, "larva"= 21)
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
bc_pnts <- ape::pcoa(bc_dist)
jacc_pnts <- ape::pcoa(jac_dist)
wtd_frc_pnts <- ape::pcoa(wtd_frc_dist)
uwt_frc_pnts <- ape::pcoa(uwt_frc_dist)

# Now make the ordination so that you can look at it (i.e., so that it prints to the viewer in R)
supportR::pcoa_ord(mod = bc_pnts, groupcol = beta$Stage.Gut,
                   leg_cont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
                   colors = all_cols, leg_pos = "bottomleft")

# Once that looks good save it out using this dev.off bit
tiff(file = file.path("graphs", "PCoA-bray-curtis.tiff"))
supportR::pcoa_ord(mod = bc_pnts, groupcol = beta$Stage.Gut,
                leg_cont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
                colors = all_cols, leg_pos = "bottomleft")
dev.off()

# Now make the Jaccard PCoA
supportR::pcoa_ord(mod = jacc_pnts, groupcol = beta$Stage.Gut,
                leg_cont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
                colors = all_cols, leg_pos = "bottomleft")

# Once that looks good save it out using this dev.off bit
tiff(file = file.path("Graphs", "PCoA-jaccard.tiff"))
supportR::pcoa_ord(mod = jacc_pnts, groupcol = beta$Stage.Gut,
                leg_cont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
                colors = all_cols, leg_pos = "bottomleft")
dev.off()

# Now make the Weighted Unifrac PCoA
supportR::pcoa_ord(mod = wtd_frc_pnts, groupcol = frc_data$Stage.Gut,
                leg_cont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
                colors = all_cols, leg_pos = "topright")

# Save the ordination
tiff(file = file.path("Graphs", "PCoA-weighted-unifrac.tiff"))
supportR::pcoa_ord(mod = wtd_frc_pnts, groupcol = frc_data$Stage.Gut,
                leg_cont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
                colors = all_cols, leg_pos = "topright")
dev.off()

# Finally, do the Unweighted Unifrac PCoA
supportR::pcoa_ord(mod = uwt_frc_pnts, groupcol = frc_data$Stage.Gut,
                leg_cont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
                colors = all_cols, leg_pos = "topleft")

# Save the ordination
tiff(file = file.path("Graphs", "PCoA-unweighted-unifrac.tiff"))
supportR::pcoa_ord(mod = uwt_frc_pnts, groupcol = frc_data$Stage.Gut,
                leg_cont = c("Paunch", "Ileum", "L-Mid", "A-Mid", "A-Hind"),
                colors = all_cols, leg_pos = "topleft")
dev.off()

## -------------------------------------------- ##
        # Alpha Diversity Bar Graphs ####
## -------------------------------------------- ##
# Incl. Chao, ACE, Simpson, and Shannon diversity & Pielou's Evenness
  ## (Diversity ~ Life Stage + Gut Section, data = all)

# Summarize the response variables
alpha_plotdf <- alpha %>%
  dplyr::group_by(Stage.Gut) %>%
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
                                   levels = c("Adult hindgut", "Adult midgut",
                                              "Larval midgut", "Larval ileum",
                                              "Larval paunch")))

# Shannon diversity bar graph
ggplot(alpha_plotdf, aes(y = shan, x = Stage.Gut, fill = Stage.Gut, color = 'x')) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = shan + shan_se, ymin = shan - shan_se), width = 0.2) +
  geom_text(label = "b", x = 0.7, y =  3.8, size = 5) +
  geom_text(label = "a", x = 1.7, y =  2.3, size = 5) +
  geom_text(label = "c", x = 2.7, y =  6.15, size = 5) +
  geom_text(label = "bc", x = 3.7, y = 4.75, size = 5) +
  geom_text(label = "c", x = 4.7, y =  5.45, size = 5) +
  scale_fill_manual(values = all_cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "Shannon Diversity") +
  pref_theme + theme(legend.position = 'none')
## Export
ggsave(plot = last_plot(), file.path("graphs", "bar-shannon.tiff"),
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
  scale_fill_manual(values = all_cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "Pielou's Evenness") +
  pref_theme + theme(legend.position = 'none')
ggsave(plot = last_plot(), file.path("graphs", "bar-pielou.tiff"),
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
  scale_fill_manual(values = all_cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "Simpson Diversity") +
  pref_theme + theme(legend.position = 'none')
ggsave(plot = last_plot(), file.path("graphs", "bar-simpson.tiff"),
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
  scale_fill_manual(values = all_cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "ACE") +
  pref_theme + theme(legend.position = 'none')
ggsave(plot = last_plot(), file.path("graphs", "bar-ace.tiff"),
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
  scale_fill_manual(values = all_cols) +
  scale_color_manual(values = 'black') +
  labs(x = "Life Stage & Gut Region", y = "Chao1 Index") +
  pref_theme + theme(legend.position = 'none')
ggsave(plot = last_plot(), file.path("graphs", "bar-chao.tiff"),
       width = 5, height = 5, units = 'in')

## -------------------------------------------- ##
        # Family Abundance Bar Graphs ####
## -------------------------------------------- ##
# Summarize data
fam_plotdf <- supportR::summary_table(data = fams,
                                      groups = c("Stage.Gut", "Family"),
                                      response = "Abundance")

# Create a graph for each gut region (separately)
for(stage_gut in sort(unique(fam_plotdf$Stage.Gut))){

  # Beginning message
  message("Begun graph for life stage and gut section: ", stage_gut)

  # Subset data to appropriate level
  sub_fam_plotdf <- fam_plotdf %>%
    dplyr::filter(Stage.Gut == stage_gut)

  # Make plot
  q <- ggplot(data = sub_fam_plotdf, aes(y = mean, x = reorder(Family, -mean),
                                     fill = Stage.Gut, color = 'x')) +
    geom_bar(stat = 'identity') +
    geom_errorbar(aes(ymax = mean + std_error, ymin = mean - std_error), width = 0.2) +
    scale_fill_manual(values = all_cols) +
    scale_color_manual(values = 'black') +
    labs(x = "Family", y = paste0(stringr::str_to_title(stage_gut), " Abundance")) +
    pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))

  # Export
  ggsave(plot = q, file.path("graphs", paste0("fam-abun-", tolower(gsub(" ", "-", stage_gut)),
                                              ".tiff")),
         width = 6, height = 6, units = 'in') }

## -------------------------------------------- ##
        # Phylum Abundance Bar Graphs ####
## -------------------------------------------- ##
# Summarize data
phyla_plotdf <- supportR::summary_table(data = phyla,
                                        groups = c("Stage.Gut", "Phylum"),
                                        response = "Abundance")

# Create a graph for each gut region (separately)
for(stage_gut in sort(unique(phyla_plotdf$Stage.Gut))){

  # Beginning message
  message("Begun graph for life stage and gut section: ", stage_gut)

  # Subset data to appropriate level
  sub_phyla_plotdf <- phyla_plotdf %>%
    dplyr::filter(Stage.Gut == stage_gut)

  # Make plot
  q <- ggplot(data = sub_phyla_plotdf, aes(y = mean, x = reorder(Phylum, -mean),
                                         fill = Stage.Gut, color = 'x')) +
    geom_bar(stat = 'identity') +
    geom_errorbar(aes(ymax = mean + std_error, ymin = mean - std_error), width = 0.2) +
    scale_fill_manual(values = all_cols) +
    scale_color_manual(values = 'black') +
    labs(x = "Phylum", y = paste0(stringr::str_to_title(stage_gut), " Abundance")) +
    pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))

  # Export
  ggsave(plot = q,
         file.path("graphs", paste0("phyla-abun-", tolower(gsub(" ", "-", stage_gut)), ".tiff")),
         width = 6, height = 6, units = 'in') }

## -------------------------------------------- ##
       # Global Phylum Abundance Plot ####
## -------------------------------------------- ##
# Get a global version of that (i.e., one that ditches the sample information)
phyl_global <- phyla %>%
  # Get total sum per phylum
  dplyr::group_by(Phylum) %>%
  dplyr::summarise(Abundance = sum(Abundance, na.rm = T)) %>%
  dplyr::ungroup()

# Look at what that produced
dplyr::glimpse(phyl_global)

# Graph that
ggplot(phyl_global, aes(y = Abundance, x = reorder(Phylum, -Abundance),
                       fill = 'x', color = 'x')) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = "#41ab5d") +
  scale_color_manual(values = 'black') +
  labs(x = "Phylum", y = "Abundance") +
  pref_theme + theme(legend.position = 'none', axis.text.x = element_text(size = 8))

# Export
ggsave(plot = last_plot(), file.path("graphs", "phyla-abun-all.tiff"),
       width = 6, height = 6, units = 'in')

# END ####
