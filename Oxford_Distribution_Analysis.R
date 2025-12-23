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

# Load data ---------------------------------------------------------------
# Load entire dataset in long format
trees_long_df <- read.csv("Oxford_Trees_Long.csv")

# Load tree metadata
trees_metadata_df <- read.csv("Oxford_Trees_Metadata_Wide.csv", header = TRUE) %>%
                     select(-starts_with("X"))

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

# Distribution analysis ----------------------------------------------------
# Plot discrete intensity in most recent year for each tree - density plot
ggplot(data = combined_df, aes(x = Intensity)) +
  geom_density(fill = "lightblue", alpha = 0.6) +
  theme_bw()

# Plot discrete intensity in most recent year for each tree - histogram
ggplot(data = combined_df, aes(x = Intensity)) +
  geom_histogram(fill = "lightblue", alpha = 0.6) +
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

# Remove 0 counts
combined_pos_df <- filter(combined_df, Intensity > 0)

# Add in non-infected hosts to test for prevalence

# Examine fit to truncated Poisson model
fit_pois <- glmmTMB(
  Intensity ~ 1,
  family = truncated_poisson,
  data = combined_pos_df
)
summary(fit_pois)

# Fit to negative binomial
fit_nb <- glmmTMB(
  Intensity ~ 1,
  family = truncated_nbinom2,
  data = combined_pos_df
)
summary(fit_nb)

# Compare AIC
AIC(fit_pois, fit_nb)
# Truncated NB model fits better!

# SHOW PLOTS OF NB AND POISSON FIT

# Impact of height
ht_int_nb <- glmmTMB(
  Intensity ~ Height,
  family = truncated_nbinom2,
  data = combined_pos_df
)
summary(ht_int_nb)



# Impact of hosts on population-level processes -----------------------------------------------------
# Test impact of host species and host size on recruitment, growth and fruiting


# Individual-based model --------------------------------------------------
# Define host-dependent processes

# Define steps in model

# Forecast with uncertainty


