### Oxford Initial Analysis
### 2025.11.25
### Oliver G. Spacey

# Analysis of mistletoe distribution with a focus on aggregation across hosts, effects of host species and size on intensity, growth/survival and fruiting

# Clear environment
rm(list = ls())

# Set working directory

# Load required packages
library(tidyverse) # Data wrangling
library(ggplot2)   # Data visualisation
library(glmmTMB)   # Fit truncated distribution
library(RANN)      # Nearest neighbour calculation
library(ggbreak)   # Plotting with scale breaks
library(lme4)      # Generalised linear mixed models

# Set seed for reproducibilty
set.seed(16)

# Load data ---------------------------------------------------------------
# Load entire dataset in long format
trees_long_df <- read.csv("Oxford_Trees_Long.csv")

# Load tree metadata
trees_metadata_df <- read.csv("Oxford_Trees_Metadata_Wide.csv", header = TRUE) %>%
                     select(-starts_with("X"))

# Load coordinates of all trees
all_trees_coo_df <- read.csv("bluesky_trees_bng_with_coords.csv") %>%
                      rename(easting = easting_bng, northing = northing_bng, date = Date_1) %>%
                      select(-c(LAYER, Date_2, Date_3, Date_4, Copyright,
                                starts_with("OS"))) %>%
                      mutate(Index = 1:nrow(read.csv("bluesky_trees_bng_with_coords.csv")))
dim(all_trees_coo_df)

# Load coordinates of host trees
hst_trees_coo_df <- read.csv("host_trees_bng_with_coords.csv") %>%
                      rename(TreeID = Name,
                             easting = easting_bng,
                             northing = northing_bng) %>%
                      mutate(Index = 1:nrow(read.csv("host_trees_bng_with_coords.csv")))

# ADD IRENE'S DATA

# Wrangle data ------------------------------------------------------------
# Select relevant columns for analysis
trees_sel_long_df <- select(trees_long_df,
                              -Cir,
                              -DBH,
                              -starts_with("Base"),
                              -starts_with("Top"),
                              -starts_with("Height"),
                              -Crown_size,
                              -Entire_tree,
                              -starts_with("X"))

# ACCOUNT FOR MISSING MISTLETOES - 
# Measure uncertainty across sample?

# For now, take most recent census of each tree

# Select only most recent season where individuals were counted to get snapshot of population
trees_fil_long_df <- trees_sel_long_df %>%
  filter(Indiv_mst == 1) %>%           # keep only rows where individuals were counted
  group_by(TreeID) %>%
  filter(Season == max(Season)) %>%  # most recent counted season
  ungroup()

# Combine data frames into single data frame
combined_df <- trees_fil_long_df %>%
  full_join(trees_metadata_df, by = "TreeID")

#. Make intensity numeric
combined_df$Intensity <- as.numeric(combined_df$Intensity)

# Investigate correlation between Height, DBH and Genus
ggplot(data = combined_df, aes(x = DBH, y = Height, col = Genus)) +
  geom_point() +
  theme_bw()

# Match up host tree coordinates with those in all trees dataset
# Example input names (replace with your actual data frames)
# A <- data.frame(tree_id = ..., easting = ..., northing = ...)
# B <- data.frame(b_id = ..., easting = ..., northing = ...)

# Build coordinate matrices
hst_xy <- as.matrix(hst_trees_coo_df[, c("easting", "northing")])
all_xy <- as.matrix(all_trees_coo_df[, c("easting", "northing")])

# Find nearest neighbour in B for each A (fast kd-tree)
nn <- nn2(data = all_xy, query = hst_xy, k = 1)

# nn$nn.idx is an integer vector: for each host row, the index of nearest tree
nearest_idx <- as.integer(nn$nn.idx[,1])
distances <- as.numeric(nn$nn.dists[,1])  # Euclidean distance (same units as easting/northing)

# Output host-centric result: which trees each host is matched to, distance, and duplicate flag on B later
matched_hsts <- hst_trees_coo_df %>%
  mutate(
    matched_tree_row = nearest_idx,
    matched_tree_id = all_trees_coo_df$Index[nearest_idx],
    distance = distances
  )

# Flag assigned count and mark duplicates
assigned_counts <- table(matched_hsts$matched_tree_id)          # named table
matched_hsts <- matched_hsts %>%
  mutate(
    tree_assigned_count = as.integer(assigned_counts[as.character(matched_tree_id)]),
    duplicate_flag = ifelse(is.na(tree_assigned_count), FALSE, tree_assigned_count > 1)
  )

# All trees summary: for every tree, how many hosts matched, list of host ids, minimum distance
all_trees_summary <- matched_hsts %>%
  group_by(matched_tree_id) %>%
  summarise(
    n_hst_assigned = n(),
    hst_ids = paste(TreeID, collapse = ", "),
    min_distance = min(distance),
    .groups = "drop"
  ) %>%
  right_join(all_trees_coo_df %>% rename(matched_tree_id = Index), by = "matched_tree_id") %>%
  mutate(
    n_hst_assigned = replace_na(n_hst_assigned, 0),
    hst_ids = replace_na(hst_ids, ""),
    is_duplicate = n_hst_assigned > 1
  )

# Output objects:
# matched_hsts : one row per host (tree_id) with matched_tree_id, distance, tree_assigned_count, duplicate_flag
# all_trees_summary : one row per tree with how many hosts assigned and which hosts (if any)

matched_hsts
all_trees_summary

# Ignore duplicate trees for now, selecting closest nearest neighbour - CHANGE THIS LATER
unmatched_trees <- c("WAD2", "CHC38", "MAG52", "MAG55", "MAG75", "MAG79", "MER13", "MER14", "HIL4")

# Create negated "in" function
`%nin%` = Negate(`%in%`)

# Remove unwanted duplicates
matched_hsts_dup_rm <- matched_hsts %>%
                       filter(TreeID %nin% unmatched_trees)

# Recreate all trees summary: for every tree, how many hosts matched, list of host ids, minimum distance
all_trees_summary_dup_rm <- matched_hsts_dup_rm %>%
  group_by(matched_tree_id) %>%
  summarise(
    n_hst_assigned = n(),
    hst_ids = paste(TreeID, collapse = ", "),
    min_distance = min(distance),
    .groups = "drop"
  ) %>%
  right_join(all_trees_coo_df %>% rename(matched_tree_id = Index), by = "matched_tree_id")

# Estimate UNINFECTED HOSTS from all points - designate non-host species
# According to i-Tree Eco study, the 20 most common species constituted 78% of the trees in Oxford
# Of this 78% of trees, 9.27% are non-host species (including oaks and conifers)
# Assuming the same proportion of non-hosts, estimate that 11.86% of trees are non-hosts
p_nonhst <- 0.1186

# Designate trees with mistletoe as hosts
all_trees_df <- all_trees_summary_dup_rm %>%
             mutate(hst_status = if_else(!is.na(n_hst_assigned), "Host", "Non_host"))

# Calculate number of trees
n_trees <- nrow(all_trees_df)

# Calculate proportion of trees that are infested
n_inf <- as.numeric(all_trees_df %>%
                      filter(hst_status == "Host") %>%
                      summarise(count = n()))
p_inf <- n_inf / n_trees

# Calculate probability of an uninfected tree being a potential host (i.e., it is a known host species of Viscum album)
p_uninf <- 1 - p_inf
p_uninf_hst <- 1 - p_uninf * p_nonhst

# Calculate number of uninfected hosts that are potential hosts
n_uninf_hst <- floor((n_trees - n_inf) * p_uninf_hst)

# Distribution analysis ----------------------------------------------------
# Ignore uninfected hosts for now

# Plot discrete intensity in most recent year for each tree - density plot
ggplot(data = combined_df, aes(x = Intensity)) +
  geom_density(fill = "lightblue", alpha = 0.6) +
  theme_bw()

# Plot discrete intensity in most recent year for each tree - histogram
ggplot(data = combined_df, aes(x = Intensity)) +
  geom_histogram(fill = "lightblue", alpha = 0.6, binwidth = 1) +
  theme_bw()

# FIGURE OUT HOW TO DEAL WITH "TOO MANY TO COUNT"
# Use discrete classes of parasite load?

# Plot intensity by genus
ggplot(data = combined_df, aes(x = Genus, y = Intensity, col = Genus)) +
  geom_point() +
  geom_violin() +
  theme_bw()

# Plot intensity by height
ggplot(data = combined_df, aes(x = Height, y = Intensity)) +
  geom_point() +
  geom_smooth() +
  theme_bw()

# Plot intensity by DBH
ggplot(data = combined_df, aes(x = DBH, y = Intensity)) +
  geom_point() +
  geom_smooth() +
  theme_bw()

# Plot intensity by Height * DBH
ggplot(data = combined_df, aes(x = Height*DBH, y = Intensity)) +
  geom_point() +
  geom_smooth() +
  theme_bw()


# Fit to Poisson and NBDs -------------------------------------------------

# Add in uninfected hosts to test for prevalence - remove NAs
all_hsts_df <- combined_df %>%
  bind_rows(tibble(Intensity = rep(0, n_uninf_hst))) %>%
  drop_na(Intensity)

# Define intensities
y <- all_hsts_df$Intensity
n <- length(y)
k <- 0:max(y)

# Examine fit to zero-inflated Poisson model
fit_pois <- glmmTMB(
  Intensity ~ 1,
  ziformula = ~ 1,
  family = poisson,
  data = all_hsts_df
)
summary(fit_pois)

# Fit to zero-inflated negative binomial
fit_nb <- glmmTMB(
  Intensity ~ 1,
  ziformula = ~ 1,
  family = nbinom2,
  data = all_hsts_df
)
summary(fit_nb)

# Compare AIC
AIC(fit_pois, fit_nb)
# Truncated NB model fits better!

# SHOW PLOTS OF NB AND POISSON FIT
# Extract Poisson parameters
lambda_pois <- exp(fixef(fit_pois)$cond[1])
pi_pois     <- plogis(fixef(fit_pois)$zi[1])   # zero-inflation probability

# Extract negative binomial parameters
mu_nb    <- exp(fixef(fit_nb)$cond[1])
pi_nb    <- plogis(fixef(fit_nb)$zi[1])
theta_nb <- sigma(fit_nb)^-2   # NB2 size parameter

# Define PMFs for each distribution
# Zero-inflated Poisson
dzip <- function(k, lambda, pi) {
  ifelse(
    k == 0,
    pi + (1 - pi) * dpois(0, lambda),
    (1 - pi) * dpois(k, lambda)
  )
}


# Zero-inflated negative binomial
dzinb <- function(k, mu, theta, pi) {
  ifelse(
    k == 0,
    pi + (1 - pi) * dnbinom(0, mu = mu, size = theta),
    (1 - pi) * dnbinom(k, mu = mu, size = theta)
  )
}

# Compute expected frequencies under these distributions
pmf_df <- tibble(
  k = k,
  zip  = dzip(k,  lambda_pois, pi_pois) * n,
  zinb = dzinb(k, mu_nb, theta_nb, pi_nb) * n
)

# Plot fit of distributions to abundance data
ggplot(data = tibble(y = all_hsts_df$Intensity), aes(x = y)) +
  lims(x = c(1, max(y)), y = c(0, max(subset(pmf_df, k == 1)))) +
  geom_histogram(
    binwidth = 1,
    boundary = -0.5,
    fill = "grey85",
    color = "black"
  ) +
  geom_line(
    data = pmf_df,
    aes(x = k, y = zip, color = "Poisson"),
    linewidth = 1
  ) +
  geom_line(
    data = pmf_df,
    aes(x = k, y = zinb, color = "NB"),
    linewidth = 1
  ) +
  geom_point(
    data = pmf_df,
    aes(x = k, y = zip, color = "Poisson"),
    size = 2
  ) +
  geom_point(
    data = pmf_df,
    aes(x = k, y = zinb, color = "NB"),
    size = 2
  ) +
  scale_color_manual(
    values = c("Poisson" = "blue", "NB" = "red")
  ) +
  labs(
    x = "Intensity",
    y = "Frequency",
    color = "Model",
    title = "Observed counts with fitted zero-inflated distributions"
  ) +
  lims(y = c(0, 100)) +
  theme_minimal()

# Extract number of uninfected hosts from two distributions
subset(pmf_df, k == 0)

# Impact of height
ht_int_nb <- glmmTMB(
  Intensity ~ Height,
  family = truncated_nbinom2,
  data = combined_pos_df
)
summary(ht_int_nb)

# Impact of DBH


# Impact of hosts on population-level processes -----------------------------------------------------
# Test impact of host species and host size on recruitment, growth and fruiting
# Add metadata to long df
trees_sel_long_df <- full_join(trees_sel_long_df, trees_metadata_df)

# Convert intensity to numeric variable
trees_sel_long_df$Intensity <- as.numeric(trees_sel_long_df$Intensity)

# Remove genera with < 3 trees
# Find genera with < 3 trees or NA
rare_genera <- trees_metadata_df %>%
  count(Genus) %>%
  filter(n < 3) %>%
  pull(Genus)

# Calculate intensity growth rate (= change in intensity) for each tree each year
int_gro_df <- trees_sel_long_df %>%
  arrange(TreeID, Season) %>%            # chronological order per tree
  group_by(TreeID) %>%
  mutate(
    next_Census    = lead(Census),
    next_Season    = lead(Season),
    next_Intensity = lead(Intensity)
  ) %>%
  # keep rows where this season AND the next season were both measured (Census == 1)
  filter(Census == 1, next_Census == 1) %>%
  transmute(
    TreeID,
    Interval = paste0(Season, "->", next_Season),
    Intensity_growth = log(next_Intensity) - log(Intensity),
    Intensity = Intensity
  ) %>%
  ungroup() %>%
  full_join(trees_metadata_df) %>%
  subset(Genus %nin% rare_genera)

# Remove -Inf
int_gro_df <- int_gro_df %>%
  mutate(Intensity_growth = ifelse(is.finite(Intensity_growth), Intensity_growth, NA)) %>%
  subset(!is.na(Genus))


# Plot growth rates as function of genus and height
ggplot(data = int_gro_df, aes(x = Height, y = Intensity_growth, col = Genus)) +
  geom_point() +
  theme_bw()
ggplot(data = int_gro_df, aes(x = DBH, y = Intensity_growth, col = Genus)) +
  geom_point() +
  theme_bw()




# GLMMs of host genus, DBH and height on growth  
lmer(Intensity_growth ~ Genus + (1|TreeID), data = int_gro_df)

# GLMMs of host genus, DBH and height on growth  
lmer(Intensity_growth ~ Height + (1|TreeID), data = int_gro_df)

# Create combined data frame to test effects on fruiting
trees_fru_df <- full_join(trees_sel_long_df, trees_metadata_df) %>%
  subset(Genus %nin% rare_genera) %>%
  subset(!is.na(Genus))
trees_fru_df$PF <- as.numeric(trees_fru_df$PF)

# Plot proportion fruiting as function of genus and height
ggplot(data = trees_fru_df, aes(x = Genus, y = PF, col = Genus)) +
  geom_point() +
  geom_violin() +
  theme_bw()
ggplot(data = trees_fru_df, aes(x = Height, y = PF)) +
  geom_point() +
  geom_smooth() +
  theme_bw()
ggplot(data = trees_fru_df, aes(x = DBH, y = PF)) +
  geom_point() +
  geom_smooth() +
  theme_bw()

# GLMMs of host genus, DBH and height on fruiting  


# Test density-dependent effects; relationship between intensity and fruiting, intensity and growth
ggplot(data = trees_fru_df, aes(x = Intensity, y = PF)) +
  geom_point() +
  geom_smooth() +
  theme_bw()

ggplot(data = int_gro_df, aes(x = Intensity, y = Intensity_growth)) +
  geom_point() +
  geom_smooth() +
  theme_bw()

# Individual-based model --------------------------------------------------
# Define host-dependent processes - make height dependent IBM first

# Define steps in model

# Forecast with uncertainty


