### Oxford Distribution Analysis
### 2025.11.25
### Oliver G. Spacey, Erola Fenollosa Romani

# Analysis of mistletoe distribution with a focus on aggregation across hosts, effects of host species and size on intensity, growth/survival and fruiting

# Clear environment
rm(list = ls())

# Set working directory
getwd()

# Load required packages
library(tidyverse) # Data wrangling
library(ggplot2)   # Data visualisation
library(glmmTMB)   # Fit truncated distribution
library(RANN)      # Nearest neighbour calculation
library(ggbreak)   # Plotting with scale breaks
library(lme4)      # Generalised linear mixed models
library(agricolae) # Tukey-adjusted post-hoc comparisons
library(sjPlot)    # Plot model outputs
library(igraph)    # Network construction
library(MASS)      # Negative binomial GLM

# Set seed for reproducibilty
set.seed(16)

# Create negated "in" function
`%nin%` = Negate(`%in%`)

# Load data ---------------------------------------------------------------
# Load discrete intensity dataset in long format
trees_discrete_long_df <- read.csv("Oxford_Trees_Long.csv") %>%
  rename(discrete_intensity = Intensity) %>%
  dplyr::select(-starts_with("X"))

# Rewrite Season for consistency with other dataset and create distinct year x tree label
trees_discrete_long_df$Season <- gsub("_", "-", paste("20", trees_discrete_long_df$Season, sep = ""))
trees_discrete_long_df <- mutate(trees_discrete_long_df, Label = paste(trees_discrete_long_df$TreeID, "_", trees_discrete_long_df$Season, sep = ""))

# Make intensity numeric
trees_discrete_long_df$discrete_intensity <- as.numeric(trees_discrete_long_df$discrete_intensity)

# Load proportional intensity dataset (from imageJ analysis)
trees_imageJ_df <- read.csv("MASTER_data_wide.csv") %>%
                        dplyr::select(-c(TreeID, Tree))

# Load tree metadata
trees_metadata_df <- read.csv("Oxford_Trees_Metadata_Wide.csv", header = TRUE)

# Load coordinates of all trees
all_trees_coo_df <- read.csv("bluesky_trees_bng_with_coords.csv") %>%
                      rename(easting = easting_bng, northing = northing_bng, date = Date_1) %>%
                      dplyr::select(-c(LAYER, Date_2, Date_3, Date_4, Copyright,
                                starts_with("OS"))) %>%
                      mutate(Index = 1:nrow(read.csv("bluesky_trees_bng_with_coords.csv")))
dim(all_trees_coo_df)

# Load coordinates of host trees
hst_trees_coo_df <- read.csv("host_trees_bng_with_coords.csv") %>%
                      rename(TreeID = Name,
                             easting = easting_bng,
                             northing = northing_bng) %>%
                      mutate(Index = 1:nrow(read.csv("host_trees_bng_with_coords.csv")))

# Load carbon and nitrogen datasets
C_df <- read.csv("Results_Carbon_220506.csv") %>%
        rename("C_content" = X.C,
               "TreeID" = Sample.ID)
N_df <- read.csv("Results_Nitrogen_220531.csv") %>%
        rename("N_content" = X.N,
               "TreeID" = Sample.ID)
N_df$TreeID <- substr(N_df$TreeID, 1, nchar(N_df$TreeID) - 1)

# Load TWI data

# Load distance from roads data

# Load distance from path observed data

# Load seed experiment data

# Wrangle data ------------------------------------------------------------
# Select relevant columns for analysis from discrete dataset
trees_discrete_long_df <- dplyr::select(trees_discrete_long_df,
                              -Cir,
                              -DBH,
                              -starts_with("Base"),
                              -starts_with("Top"),
                              -starts_with("Height"),
                              -Crown_size,
                              -Entire_tree)

# For imageJ data, calculate:
# Estimated discrete intensity
trees_imageJ_df <- trees_imageJ_df %>%
  mutate(n_mistletoe_imageJ = rowSums(!is.na(dplyr::select(., starts_with("M") & ends_with("_Area")))))

# Total mistletoe area
trees_imageJ_df <- trees_imageJ_df %>%
  mutate(total_mistletoe_area = rowSums(dplyr::select(., ends_with("_Area") & starts_with("M")), na.rm = TRUE))

# Estimated proportional intensity
trees_imageJ_df <- trees_imageJ_df %>%
  mutate(proportional_intensity = total_mistletoe_area / MaskArea)

# Combine data frames into single data frame
combined_df <- trees_discrete_long_df %>%
  full_join(trees_imageJ_df, by = "Label") %>%
  full_join(trees_metadata_df, by = "TreeID") %>%
  full_join(C_df, by = "TreeID")%>%
  full_join(N_df, by = "TreeID")

  # ACCOUNT FOR MISSING MISTLETOES*** 
  # Measure uncertainty across sample?

# Count number "Missed" from .pdfs and number new recruits
# OR only use 2024-25 anyway as snapshot; use previous years for vital rates
# Assume none missed in that year? Or missed at a constant rate?

  
  # For now, take most recent census of each tree
  # Select only most recent season where individuals were counted to get snapshot of population - just one season is actually more accurate here
  trees_fil_discrete_df <- combined_df %>%
    filter(Indiv_mst == 1) %>%           # keep only rows where individuals were counted
    group_by(TreeID) %>%
    filter(Season == max(Season)) %>%  # most recent counted season
    ungroup()
  
# Compare discrete intensity measured from imageJ versus measured from pdfs
  ggplot(data = combined_df, aes(x = n_mistletoe_imageJ, y = discrete_intensity)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1) +
    theme_bw()

# Investigate correlation between Height, DBH and Genus
ggplot(data = combined_df, aes(x = DBH, y = Height, col = Genus)) +
  geom_point() +
  theme_bw()

# FIND WAY OF MATCHING TREES MORE ACCURATELY***
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
# Create dataframe to compare hosts and non-hosts
hst_vs_nhst_df <- all_trees_summary_dup_rm %>%
                 mutate(host = case_when(.$n_hst_assigned == 1 ~ "Host",
                                         .default = "Non-host"))
hst_vs_nhst_df$host <- as.factor(hst_vs_nhst_df$host)

# Plot height for hosts vs non-hosts
ggplot(hst_vs_nhst_df, aes(x = host, y = Max, fill = host)) +
  geom_violin(alpha = 0.2) +
  geom_jitter(alpha = 0.1, width = 0.1, aes(col = host)) +
  stat_summary(fun = mean,
    fun.min = function(x) mean(x) - sd(x)/sqrt(length(x)),
    fun.max = function(x) mean(x) + sd(x)/sqrt(length(x)),
    geom = "errorbar",
    width = 0.15) +
  stat_summary(fun = mean,
    geom = "point",
    size = 1) +
  labs(x = "Host or non-host", y = "Maximum height") +
  theme_bw()

ggplot(data = hst_vs_nhst_df, aes(x = host, y = Mean, col = host)) +
  geom_violin() +
  geom_jitter(alpha = 0.03) +
  labs(x = "Host or non-host", y = "Mean height") +
  theme_bw()

# Perform t-test
t.test(Max ~ host, data = hst_vs_nhst_df)

# Ignore uninfected hosts for now

# CHANGE TO 2024-25 SEASON - more reliable
# Plot discrete intensity in most recent year for each tree - density plot
ggplot(data = combined_df, aes(x = Intensity)) +
  geom_density(fill = "lightblue", alpha = 0.6) +
  theme_bw()

# Plot discrete intensity in most recent year for each tree - histogram
ggplot(data = combined_df, aes(x = Intensity)) +
  lims(x = c(0, max(combined_df$Intensity)), y = c(0, 100)) +
  geom_histogram(fill = "lightblue", alpha = 0.6, binwidth = 1) +
  theme_bw()

# FIGURE OUT WHERE TO DRAW BOUNDARY BETWEEN DISCRETE AND PROPORTIONAL INTENSITY
# REPORT BOTH

# Plot whether discrete mistletoes counted as a function of proportional intensity
ggplot(data = combined_df, aes(x = proportional_intensity, y = as.numeric(Indiv_mst))) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", family="binomial") +
  theme_bw()

# Plot discrete intensity against distance from path


# Remove genera with < 3 trees
# Find genera with < 3 trees or NA
rare_genera <- combined_df %>%
  count(Genus) %>%
  filter(n < 3) %>%
  pull(Genus)

# Subset to only genera with 3 or more observations
common_gen_df <- subset(combined_df, Genus %nin% rare_genera) %>%
  subset(!is.na(Genus))
common_gen_df$Genus <- as.factor(common_gen_df$Genus)

# Plot intensity by genus
ggplot(data = common_gen_df, aes(x = Genus, y = Intensity, col = Genus)) +
  geom_boxplot(outliers = FALSE) +
  labs(x = "Genus", y = "Mistletoe intensity") +
  geom_jitter(alpha = 0.5) +
  theme_bw()

# Perform ANOVA
gen_int_aov <- aov(Intensity ~ Genus, data = common_gen_df)
summary(gen_int_aov)
summary(lm(gen_int_aov))

# Tukey post-hoc test
get_int_tukey <- HSD.test(gen_int_aov, trt = 'Genus')
get_int_tukey$groups

# Perform glm
only_inf_df <- combined_df[which(combined_df$Intensity != 0), ]
ht_int_glm <- glm.nb(Intensity ~ Height, data = only_inf_df)
summary(ht_int_glm)
plot(ht_int_glm)

# Predict from glm
ht_glm_df <- cbind(only_inf_df[, c(5,12)], "resp" = predict(ht_int_glm, only_inf_df, type = "response", se.fit = TRUE)[1:2])
ht_glm_df <- cbind(ht_glm_df, "link" = predict(ht_int_glm, only_inf_df, type = "link", se.fit = TRUE)[1:2])

# Plot intensity by height with predicted
ggplot(data = ht_glm_df, aes(x = Height)) +
  geom_point(aes(y = Intensity)) +
  geom_line(aes(y = resp.fit), color = "red") +
  geom_line(aes(y = resp.fit + 2 * resp.se.fit), color = "red", linetype = "dashed") +
  geom_line(aes(y = resp.fit - 2 * resp.se.fit), color = "red", linetype = "dashed") +
  labs(x = "Tree height (m)", y = "Mistletoe intensity") +
  theme_bw()

# Perform basic lm
ht_int_lm <- lm(Intensity ~ Height, data = combined_df[which(combined_df$Intensity != 0), ])
summary(ht_int_lm)

# Plot intensity by DBH
ggplot(data = combined_df, aes(x = DBH, y = Intensity)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()

# Plot intensity by N content
ggplot(data = combined_df[combined_df$Genus %in% c("Acer", "Crataegus", "Malus", "Populus", "Salix", "Tilia"), ],
       aes(x = C_content, y = Intensity, col = Genus)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Wood carbon content", y = "Mistletoe intensity") + 
  theme_bw()

# Perform LM
C_lm <- lm(Intensity ~ C_content + Genus, data = combined_df[which(combined_df$Genus %in% c("Acer", "Crataegus", "Malus", "Populus", "Salix", "Tilia")), ])
summary(C_lm)

# Plot intensity by N content
ggplot(data = combined_df[combined_df$Genus %in% c("Acer", "Crataegus", "Malus", "Populus", "Salix", "Tilia"), ],
                                aes(x = N_content, y = Intensity, col = Genus)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Wood nitrogen content", y = "Mistletoe intensity") + 
  theme_bw()

# Perform LMM
# N_lm <- glm(PF ~ N_content + Genus, data = combined_df[which(combined_df$Genus %in% c("Acer", "Crataegus", "Malus", "Populus", "Salix", "Tilia")), ], family = "binomial")
# summary(N_lm)

# Explore proportional intensity

# Compare distribution of proportional intensity to a beta binomial distribution

# Plot proportional intensity against absolute intensity

# Plot proportional intensity against distance from path

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
  lims(x = c(0, max(y)), y = c(0, max(subset(pmf_df, k == 1)))) +
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
    size = 1
  ) +
  geom_point(
    data = pmf_df,
    aes(x = k, y = zinb, color = "NB"),
    size = 1
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

# Examine fit to zero-truncated models
# Examine fit to zero-truncated Poisson model
fit_ztpois <- glmmTMB(
  Intensity ~ 1,
  family = truncated_poisson,
  data = all_hsts_df[which(all_hsts_df$Intensity != 0),]
)
summary(fit_ztpois)

# Fit to zero-truncated negative binomial
fit_ztnb <- glmmTMB(
  Intensity ~ 1,
  family = truncated_nbinom2,
  data = all_hsts_df[which(all_hsts_df$Intensity != 0),]
)
summary(fit_ztnb)

# Extract Poisson parameters
lambda_ztpois <- exp(fixef(fit_ztpois)$cond)

# Extract negative binomial parameters
mu_ztnb    <- exp(fixef(fit_ztnb)$cond)
theta_ztnb <- sigma(fit_ztnb)   # NB2 dispersion parameter

# Define PMFs for each distribution
# Zero-truncated Poisson
dztpois <- function(x, lambda) {
  dpois(x, lambda) / (1 - dpois(0, lambda))
}

# Zero-inflated negative binomial
dztnb <- function(x, mu, theta) {
  dnbinom(x, mu = mu, size = theta) /
    (1 - dnbinom(0, mu = mu, size = theta))
}

# Create prediction grid
ints <- seq(
  from = min(all_hsts_df$Intensity[all_hsts_df$Intensity > 0]),
  to   = max(all_hsts_df$Intensity),
  by   = 1
)

# Compute expected frequencies under these distributions
pred_df <- data.frame(
  Intensity = ints,
  ztpois = dztpois(ints, lambda_ztpois),
  ztnb   = dztnb(ints, mu_ztnb, theta_ztnb)
)


# Plot fit of distributions to abundance data
ggplot(all_hsts_df[all_hsts_df$Intensity > 0, ],
       aes(x = Intensity)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 1,
                 fill = "grey80",
                 color = "black") +
  geom_line(data = pred_df,
            aes(y = ztpois),
            color = "blue",
            linewidth = 1) +
  geom_line(data = pred_df,
            aes(y = ztnb),
            color = "red",
            linewidth = 1) +
  labs(
    y = "Probability",
    title = "Zero-truncated model fits",
    subtitle = "Blue = Poisson, Red = Negative Binomial"
  ) +
  theme_bw()


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
ggplot(data = int_gro_df, aes(x = Genus, y = Intensity_growth, col = Genus)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter() +
  labs(y = "Δ log(intensity)") +
  theme_bw()

ggplot(data = int_gro_df, aes(x = Height, y = Intensity_growth)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(y = "Δ log(intensity)") +
  theme_bw()
ggplot(data = int_gro_df, aes(x = DBH, y = Intensity_growth)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(y = "Δ log(intensity)") +
  theme_bw()

# LMMs of host genus on growth - control for pseudorep as individual ID as random effect
int_gro_gen_glmm <- lmer(Intensity_growth ~ Genus + (1|TreeID), data = int_gro_df)
summary(int_gro_gen_glmm)
plot_model(int_gro_gen_glmm)

# LMMs of host height on growth - control for pseudorep as individual ID as random effect
int_gro_ht_glmm <- lmer(Intensity_growth ~ Height + (1|TreeID), data = int_gro_df)
summary(int_gro_ht_glmm)
plot_model(int_gro_ht_glmm)

# Create combined data frame to test effects on fruiting
trees_fru_df <- full_join(trees_sel_long_df, trees_metadata_df) %>%
  subset(Genus %nin% rare_genera) %>%
  subset(!is.na(Genus))
trees_fru_df$PF <- as.numeric(trees_fru_df$PF)

# Plot proportion fruiting as function of genus and height
ggplot(data = trees_fru_df, aes(x = Genus, y = PF, col = Genus)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter() +
  labs(y = "Proportion of mistletoes fruiting") +
  theme_bw()

ggplot(data = trees_fru_df, aes(x = Height, y = PF)) +
  geom_point() +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial")) +
  labs(y = "Proportion of mistletoes fruiting") +
  theme_bw()

ggplot(data = trees_fru_df, aes(x = DBH, y = PF)) +
  geom_point() +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial")) +
  labs(y = "Proportion of mistletoes fruiting") +
  theme_bw()

# GLMMs of host genus on fruiting  
fru_gen_glmm <- glmer(PF ~ Genus + (1|TreeID), data = trees_fru_df, family = "binomial")
summary(fru_gen_glmm)

# GLMMs of host height on fruiting  
fru_ht_glmm <- glmer(PF ~ Height + (1|TreeID), data = trees_fru_df, family = "binomial")
summary(fru_ht_glmm)

# Test density-dependent effects; relationship between intensity and fruiting, intensity and growth
ggplot(data = trees_fru_df, aes(x = Intensity, y = PF)) +
  geom_point() +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial")) +
  labs(y = "Proportion of mistletoes fruiting") +
  theme_bw()

ggplot(data = int_gro_df, aes(x = Intensity, y = Intensity_growth)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(y = "Δ log(intensity)") +
  theme_bw()

# Individual-based model --------------------------------------------------
# Define host-dependent processes - make height dependent IBM first

# Define steps in model

# Forecast with uncertainty

# Generate host community ------------------------------------------------
# ADD IN POTENTIAL HOSTS - REMOVE RANDOM - CAN MIX UP AND REPEAT LATER

# Get data frame of potential hosts, their coordinates and heights
hst_ht_df <- select(all_trees_df, hst_ids, easting, northing, Max, hst_status) %>%
             rename(TreeID = hst_ids, Height = Max) 

# Select random number of non-hosts to be potential hosts according to probabilities outlined earlier
hosts <- hst_ht_df[which(hst_ht_df$hst_status == "Host"), ]
nonhosts <- hst_ht_df[which(hst_ht_df$hst_status != "Host"), ]

nonhosts_keep <- nonhosts[sample(nrow(nonhosts), n_uninf_hst), ]

# Reassign indices
pot_hst_df <- rbind(hosts, nonhosts_keep)
pot_hst_df <- mutate(pot_hst_df, Index = 1:nrow(pot_hst_df))
                    
# Sample random x trees from dataframe for computational efficiency
pot_hst_df <- pot_hst_df[sample(nrow(pot_hst_df), 500), ]

# Get number of hosts
nhst <- nrow(pot_hst_df)

# FOR COMPUTATIONAL EFFICIENCY DO THIS IN FOR LOOP - YOU CAN DO THIS NOW IT'S LESS
## Weight host network by distance
# Create empty adjacency matrix
adj_mat <- matrix(NA, nrow = nhst, ncol = nhst)

# Calculate pythagorean distances between hosts and store in adjacency matrix - get from pot_hst_df
for(i in 1:nhst){
  i_lon <- as.numeric(pot_hst_df[i,2])
  i_lat <- as.numeric(pot_hst_df[i,3])
  for(j in 1:nhst){
    j_lon <- as.numeric(pot_hst_df[j,2])
    j_lat <- as.numeric(pot_hst_df[j,3])
    adj_mat[i,j] <- sqrt((i_lon - j_lon) ^ 2 + (i_lat - j_lat) ^ 2)
  }
}

# Construct network from adjacency matrix
# hst_ntw <- graph_from_adjacency_matrix(adjmatrix = adj_mat,
#                                        mode = "undirected",
#                                        weighted = TRUE)
# Plot host network
# plot(hst_ntw)

# Extract population growth rates
# Intensity growth rate
ggplot(data = int_gro_df, aes(x = Intensity_growth)) +
  geom_histogram() +
  theme_bw()

# Probability of no change
p_no_change <- length(which(int_gro_df$Intensity_growth == 0)) / length(which(!is.na(int_gro_df$Intensity_growth)))

# Given change, mean and sd growth rate
mean_igr <- mean(int_gro_df[which(int_gro_df$Intensity_growth != 0), ]$Intensity_growth)
sd_igr   <- sd(int_gro_df[which(int_gro_df$Intensity_growth != 0), ]$Intensity_growth)

# Set time steps
run_length <- 2000

# ASK: WHAT GENERATES THE EXPECTED LEVEL OF AGGREGATION?
  
# Create array to track status of hosts
status_ar <- array(NA, dim = c(nhst, run_length))

# Set all hosts to be uninfected
status_ar[, 1] <- 0

# Randomly select host to begin infection
status_ar[sample.int(nhst, 1), 1] <- 1

# Set probability of replacement - assume 1 every x years
p_replace <- 1/50

# Set transmission constant to define reciprocal relationship with distance - will be tweaked later
# Use from previous paper?
trnm_k <- 10

# Plot probability of transmission as an inverse function of distance - like Lavorel et al. 
# distances <- seq(0, max(adj_mat), length.out = 1000)
# p_trnms <- trnm_k / distances
# plot(x = distances, y = p_trnms)

# Set carrying capacity
cc <- 100

# Set seed
set.seed(16)

# Simulate spread and growth in for loop
for(time in 1:(run_length - 1)){
  # Set next year intensity as this year's to start
  status_ar[, time + 1] <- status_ar[, time]
  # Generate list of random host numbers
  hsts <- sample(nhst)
  # Run through hosts randomly
  for(hst in hsts){
  # Check if mistletoe present on host or not
    if(status_ar[hst, time] != 0){
      
  # Change in intensity according to intensity growth rate 
  # Probability of change - otherwise will stay the same
    if(runif(1, 0, 1) <= 1 - p_no_change){
      # Given change, probability of change - *e^IGR - round to nearest integer
      gro_rt <- rnorm(1, mean = mean_igr, sd = sd_igr)
      status_ar[hst, time + 1] <- round(status_ar[hst, time] * exp(gro_rt), digits = 0)
    }
      } else {
    
  # Transmission - reciprocal function of distance - constant with respect to height for now
      # Designate infected hosts
        infecteds <- which(status_ar[, time] > 0)
      # Create list for distances
        p_trnms_ls <- list()
      # Run through infected hosts
      for(inf in infecteds){
        # Calculate distance to target host
        # Designate host index
        # index <- infecteds[i]
        # hst_lon <- pot_hst_df[index, 2]
        # hst_lat <- pot_hst_df[index, 3]
        # tar_lon <- pot_hst_df[hst, 2]
        # tar_lat <- pot_hst_df[hst, 3]
        # distance <- sqrt((hst_lon - tar_lon) ^ 2 + (hst_lat - tar_lat) ^ 2)
        p_trnms_ls[i] <- trnm_k / adj_mat[hst, inf]
        }
        # Sum probabilities of transmission
        p_trnm_sum <- sum(unlist(p_trnms_ls))
      if(runif(1, 0, 1) <= p_trnm_sum){
        # Generate number of new recruits transmitted - could be Poisson, assume only 1 for now, trnm of a berry and successful establishment
        status_ar[hst, time + 1] <- 1 
        }
        
        }
    # DISTINGUISH SEED RAIN and MAKE NEW RECRUITS PROPORTIONAL TO FRUITING IN TREES
    # KEEP TACK OF NEW RECRUITS IN SEPARATE DATA FRAME
    
    # If reached over carrying capacity set to cc
    if(status_ar[hst, time] >= cc){
      status_ar[hst, time + 1] <- cc
    }
    
  # Probability of replacement - reset parasite load to 0 with set probability
    if(runif(1, 0, 1) <= p_replace){
      status_ar[hst, time + 1] <- 0
    }

  }
  
  # Print outputs each time step
  n_infected <- sum(status_ar[, time] != 0)
  cat("The number of infected hosts at t =", time, "is", n_infected)
}

# Plot output - number infected
times <- 1:run_length
infected_ls <- list()
iod_ls <- list()
for(i in 1:run_length){
  infected_ls[i] <- sum(status_ar[, i] != 0)
  iod_ls[i] <- var(status_ar[, i]) / mean(status_ar[, i])
}

plot(times, unlist(infected_ls), type = "l")
plot(times, unlist(iod_ls), type = "l")

# Plot outputs - plot mean/variance


# Simulate height-dependent infection - not dependent on intensity -------------------------------------
# If resampling
# Sample random x trees from dataframe for computational efficiency
pot_hst_df <- pot_hst_df[sample(nrow(pot_hst_df), 500), ]

# Get number of hosts
nhst <- nrow(pot_hst_df)

# FOR COMPUTATIONAL EFFICIENCY DO THIS IN FOR LOOP - YOU CAN DO THIS NOW IT'S LESS
## Weight host network by distance
# Create empty adjacency matrix
adj_mat <- matrix(NA, nrow = nhst, ncol = nhst)

# Calculate pythagorean distances between hosts and store in adjacency matrix - get from pot_hst_df
for(i in 1:nhst){
  i_lon <- as.numeric(pot_hst_df[i,2])
  i_lat <- as.numeric(pot_hst_df[i,3])
  for(j in 1:nhst){
    j_lon <- as.numeric(pot_hst_df[j,2])
    j_lat <- as.numeric(pot_hst_df[j,3])
    adj_mat[i,j] <- sqrt((i_lon - j_lon) ^ 2 + (i_lat - j_lat) ^ 2)
  }
}

# Construct network from adjacency matrix
# hst_ntw <- graph_from_adjacency_matrix(adjmatrix = adj_mat,
#                                        mode = "undirected",
#                                        weighted = TRUE)

# Plot host network
# plot(hst_ntw)

# Extract population growth rates
# Intensity growth rate
# ggplot(data = int_gro_df, aes(x = Intensity_growth)) +
#   geom_histogram() +
#   theme_bw()

# Probability of no change
p_no_change <- length(which(int_gro_df$Intensity_growth == 0)) / length(which(!is.na(int_gro_df$Intensity_growth)))

# Given change, mean and sd growth rate
mean_igr <- mean(int_gro_df[which(int_gro_df$Intensity_growth != 0), ]$Intensity_growth)
sd_igr   <- sd(int_gro_df[which(int_gro_df$Intensity_growth != 0), ]$Intensity_growth)

# Set time steps
run_length <- 2000

# ASK: WHAT GENERATES THE EXPECTED LEVEL OF AGGREGATION?

# Create array to track status of hosts
status_ar <- array(NA, dim = c(nhst, run_length))

# Set all hosts to be uninfected
status_ar[, 1] <- 0

# Randomly select host to begin infection
status_ar[sample.int(nhst, 1), 1] <- 1

# Set probability of replacement - assume 1 every x years
p_replace <- 1/50

# Set transmission constant to define reciprocal relationship with distance - will be tweaked later
# Use from previous paper?
trnm_k <- 10

# Plot probability of transmission as an inverse function of distance - like Lavorel et al. 
# distances <- seq(0, max(adj_mat), length.out = 1000)
# p_trnms <- trnm_k / distances
# plot(x = distances, y = p_trnms)

# Set carrying capacity
cc <- 100

# Set seed
set.seed(16)

# Calculate mean_ht
mean_ht <- mean(pot_hst_df$Height, na.rm = TRUE)

# Calculate max_ht
max_ht <- max(pot_hst_df$Height, na.rm = TRUE)

# Simulate spread and growth in for loop
for(time in 1:(run_length - 1)){
  # Set next year intensity as this year's to start
  status_ar[, time + 1] <- status_ar[, time]
  # Generate list of random host numbers
  hsts <- sample(nhst)
  # Run through hosts randomly
  for(hst in hsts){
    # Get host height
    hst_ht <- as.numeric(pot_hst_df[hst, 4])
    
    # Check if mistletoe present on host or not
    if(status_ar[hst, time] != 0){
      
      # Change in intensity according to intensity growth rate 
      # Probability of change - otherwise will stay the same
      if(runif(1, 0, 1) <= 1 - p_no_change){
        # Given change, probability of change - *e^IGR - round to nearest integer
        gro_rt <- rnorm(1, mean = mean_igr, sd = sd_igr)
        status_ar[hst, time + 1] <- round(status_ar[hst, time] * exp(gro_rt), digits = 0)
      }
    } else {
      
      # Transmission - reciprocal function of distance - constant with respect to height for now
      # Designate infected hosts
      infecteds <- which(status_ar[, time] > 0)
      # Create list for distances
      p_trnms_ls <- list()

      # Run through infected hosts
      for(inf in infecteds){
        # Calculate distance to target host
        # Designate host index
        # index <- infecteds[i]
        # hst_lon <- pot_hst_df[index, 2]
        # hst_lat <- pot_hst_df[index, 3]
        # tar_lon <- pot_hst_df[hst, 2]
        # tar_lat <- pot_hst_df[hst, 3]
        # distance <- sqrt((hst_lon - tar_lon) ^ 2 + (hst_lat - tar_lat) ^ 2)
        
        p_trnms_ls[i] <- trnm_k / adj_mat[hst, inf] * hst_ht / mean_ht
      }
      # Sum probabilities of transmission
      p_trnm_sum <- sum(unlist(p_trnms_ls))
      if(runif(1, 0, 1) <= p_trnm_sum){
        # Generate number of new recruits transmitted - could be Poisson, assume only 1 for now, trnm of a berry and successful establishment
        status_ar[hst, time + 1] <- 1 
      }
      
    }
    # DISTINGUISH SEED RAIN and MAKE NEW RECRUITS PROPORTIONAL TO FRUITING IN TREES
    # KEEP TACK OF NEW RECRUITS IN SEPARATE DATA FRAME
    
    # If reached over carrying capacity set to cc
    cc_ht <- round(cc * hst_ht / max_ht)
    if(status_ar[hst, time] >= cc_ht){
      status_ar[hst, time + 1] <- cc_ht
    }
    
    # Probability of replacement - reset parasite load to 0 with set probability
    if(runif(1, 0, 1) <= p_replace){
      status_ar[hst, time + 1] <- 0
    }
    
  }
  
  # Print outputs each time step
  n_infected <- sum(status_ar[, time] != 0)
  cat("The number of infected hosts at t =", time, "is", n_infected)
}

# Plot output - number infected
times <- 1:run_length
infected_ls <- list()
iod_ls <- list()
for(i in 1:run_length){
  infected_ls[i] <- sum(status_ar[, i] != 0)
  iod_ls[i] <- var(status_ar[, i]) / mean(status_ar[, i])
}

# Plot output - prevalence and index of dispersion for infected hosts
times <- 1:run_length
infected_ls <- list()
iod_ls <- list()
for(i in 1:run_length){
  infected_ls[i] <- sum(status_ar[, i] != 0) / nhst
  iod_ls[i] <- var(status_ar[which(status_ar[,i] != 0), i]) / mean(status_ar[which(status_ar[,i] != 0), i])
}

plot(times, unlist(infected_ls), type = "l", xlab = "Time", ylab = "Prevalence", ylim = c(0,1))
abline(h = mean_prev, lty = 2)
plot(times, unlist(iod_ls), type = "l", xlab = "Time", ylab = "σ^2/μ (infecteds only)", ylim = c(0, 65))
abline(h = mean_iod, lty = 2)


# Transmission dependent on intensity but not height ----------------------------------
# If resampling
# Sample random x trees from dataframe for computational efficiency
pot_hst_df <- pot_hst_df[sample(nrow(pot_hst_df), 1000), ]

# Get number of hosts
nhst <- nrow(pot_hst_df)

# FOR COMPUTATIONAL EFFICIENCY DO THIS IN FOR LOOP - YOU CAN DO THIS NOW IT'S LESS
## Weight host network by distance
# Create empty adjacency matrix
adj_mat <- matrix(NA, nrow = nhst, ncol = nhst)

# Calculate pythagorean distances between hosts and store in adjacency matrix - get from pot_hst_df
for(i in 1:nhst){
  i_lon <- as.numeric(pot_hst_df[i,2])
  i_lat <- as.numeric(pot_hst_df[i,3])
  for(j in 1:nhst){
    j_lon <- as.numeric(pot_hst_df[j,2])
    j_lat <- as.numeric(pot_hst_df[j,3])
    adj_mat[i,j] <- sqrt((i_lon - j_lon) ^ 2 + (i_lat - j_lat) ^ 2)
  }
}

# Construct network from adjacency matrix
# hst_ntw <- graph_from_adjacency_matrix(adjmatrix = adj_mat,
#                                        mode = "undirected",
#                                        weighted = TRUE)

# Plot host network
# plot(hst_ntw)


# Probability of no change
p_no_change <- length(which(int_gro_df$Intensity_growth == 0)) / length(which(!is.na(int_gro_df$Intensity_growth)))

# Given change, mean and sd growth rate
mean_igr <- mean(int_gro_df[which(int_gro_df$Intensity_growth != 0), ]$Intensity_growth)
sd_igr   <- sd(int_gro_df[which(int_gro_df$Intensity_growth != 0), ]$Intensity_growth)

# Set time steps
run_length <- 2000

# ASK: WHAT GENERATES THE EXPECTED LEVEL OF AGGREGATION?

# Create array to track status of hosts
status_ar <- array(NA, dim = c(nhst, run_length))

# Set all hosts to be uninfected
status_ar[, 1] <- 0

# Randomly select host to begin infection
status_ar[sample.int(nhst, 1), 1] <- 1

# Set probability of replacement - assume 1 every x years
p_replace <- 1/50

# Set transmission constant to define reciprocal relationship with distance - will be tweaked later
# Use from previous paper?
trnm_k <- 10

# Plot probability of transmission as an inverse function of distance - like Lavorel et al. 
# distances <- seq(0, max(adj_mat), length.out = 1000)
# p_trnms <- trnm_k / distances
# plot(x = distances, y = p_trnms)

# Set carrying capacity
cc <- 100

# Set seed
set.seed(16)

# Calculate mean height
mean_ht <- mean(pot_hst_df$Height, na.rm = TRUE)

# Calculate maximum height
max_ht <- max(pot_hst_df$Height, na.rm = TRUE)

# Simulate spread and growth in for loop
for(time in 1:(run_length - 1)){
  # Set next year intensity as this year's to start
  status_ar[, time + 1] <- status_ar[, time]
  # Generate list of random host numbers
  hsts <- sample(nhst)
  # Run through hosts randomly
  for(hst in hsts){
    # Get host height
    hst_ht <- as.numeric(pot_hst_df[hst, 4])
    
    # Check if mistletoe present on host or not
    if(status_ar[hst, time] != 0){
      
      # Change in intensity according to intensity growth rate 
      # Probability of change - otherwise will stay the same
      if(runif(1, 0, 1) <= 1 - p_no_change){
        # Given change, probability of change - *e^IGR - round to nearest integer
        gro_rt <- rnorm(1, mean = mean_igr, sd = sd_igr)
        status_ar[hst, time + 1] <- round(status_ar[hst, time] * exp(gro_rt), digits = 0)
      }
    } else {
      
      # Transmission - reciprocal function of distance - constant with respect to height for now
      # Designate infected hosts
      infecteds <- which(status_ar[, time] > 0)
      # Create list for distances
      p_trnms_ls <- list()
      
      # Run through infected hosts
      for(inf in infecteds){
        # Calculate distance to target host
        # Designate host index
        # index <- infecteds[i]
        # hst_lon <- pot_hst_df[index, 2]
        # hst_lat <- pot_hst_df[index, 3]
        # tar_lon <- pot_hst_df[hst, 2]
        # tar_lat <- pot_hst_df[hst, 3]
        # distance <- sqrt((hst_lon - tar_lon) ^ 2 + (hst_lat - tar_lat) ^ 2)
        
        p_trnms_ls[i] <- trnm_k / adj_mat[hst, inf] * log(status_ar[inf, time])
      }
      # Sum probabilities of transmission
      p_trnm_sum <- sum(unlist(p_trnms_ls))
      if(runif(1, 0, 1) <= p_trnm_sum){
        # Generate number of new recruits transmitted - could be Poisson, assume only 1 for now, trnm of a berry and successful establishment
        status_ar[hst, time + 1] <- 1 
      }
      
    }
    # DISTINGUISH SEED RAIN and MAKE NEW RECRUITS PROPORTIONAL TO FRUITING IN TREES
    # KEEP TACK OF NEW RECRUITS IN SEPARATE DATA FRAME
    
    # If reached over carrying capacity set to cc
    cc_ht <- round(cc * hst_ht / max_ht)
    if(status_ar[hst, time] >= cc_ht){
      status_ar[hst, time + 1] <- cc_ht
    }
    
    # Probability of replacement - reset parasite load to 0 with set probability
    if(runif(1, 0, 1) <= p_replace){
      status_ar[hst, time + 1] <- 0
    }
    
  }
  
  # Print outputs each time step
  n_infected <- sum(status_ar[, time] != 0)
  cat("The number of infected hosts at t =", time, "is", n_infected)
}

# Sample random x trees from dataframe to calculate number infected and index of dispersion
prev_ls <- list()
pop_iod_ls <- list()
for(i in 1:1000){
  random_hsts_df <- all_hsts_df[sample(nrow(all_hsts_df), nhst), ]
  prev_ls[i] <- sum(random_hsts_df$Intensity != 0) / nhst
  pop_iod_ls[i] <- var(subset(random_hsts_df, Intensity != 0)$Intensity) / mean(subset(random_hsts_df, Intensity != 0)$Intensity)
}

# Calculate mean prevalence and population index of dispersion
mean_prev <- mean(unlist(prev_ls), na.rm = TRUE)
mean_iod <- mean(unlist(pop_iod_ls), na.rm = TRUE)

# Plot output - prevalence and index of dispersion for infected hosts
times <- 1:run_length
infected_ls <- list()
iod_ls <- list()
for(i in 1:run_length){
  infected_ls[i] <- sum(status_ar[, i] != 0) / nhst
  iod_ls[i] <- var(status_ar[which(status_ar[,i] != 0), i]) / mean(status_ar[which(status_ar[,i] != 0), i])
}

plot(times, unlist(infected_ls), type = "l", xlab = "Time", ylab = "Prevalence", ylim = c(0, 1))
abline(h = mean_prev, lty = 2)
plot(times, unlist(iod_ls), type = "l", xlab = "Time", ylab = "σ^2/μ (infecteds only)", ylim = c(0, 65))
abline(h = mean_iod, lty = 2)

# Transmission dependent on height and intensity ------------------------
# If resampling
# Sample random x trees from dataframe for computational efficiency
pot_hst_df <- pot_hst_df[sample(nrow(pot_hst_df), nhst), ]

# Get number of hosts
nhst <- nrow(pot_hst_df)

# FOR COMPUTATIONAL EFFICIENCY DO THIS IN FOR LOOP - YOU CAN DO THIS NOW IT'S LESS
## Weight host network by distance
# Create empty adjacency matrix
adj_mat <- matrix(NA, nrow = nhst, ncol = nhst)

# Calculate pythagorean distances between hosts and store in adjacency matrix - get from pot_hst_df
for(i in 1:nhst){
  i_lon <- as.numeric(pot_hst_df[i,2])
  i_lat <- as.numeric(pot_hst_df[i,3])
  for(j in 1:nhst){
    j_lon <- as.numeric(pot_hst_df[j,2])
    j_lat <- as.numeric(pot_hst_df[j,3])
    adj_mat[i,j] <- sqrt((i_lon - j_lon) ^ 2 + (i_lat - j_lat) ^ 2)
  }
}

# Construct network from adjacency matrix
# hst_ntw <- graph_from_adjacency_matrix(adjmatrix = adj_mat,
#                                        mode = "undirected",
#                                        weighted = TRUE)

# Plot host network
# plot(hst_ntw)

# Extract population growth rates
# Intensity growth rate
# ggplot(data = int_gro_df, aes(x = Intensity_growth)) +
#   geom_histogram() +
#   theme_bw()

# Probability of no change
p_no_change <- length(which(int_gro_df$Intensity_growth == 0)) / length(which(!is.na(int_gro_df$Intensity_growth)))

# Given change, mean and sd growth rate
mean_igr <- mean(int_gro_df[which(int_gro_df$Intensity_growth != 0), ]$Intensity_growth)
sd_igr   <- sd(int_gro_df[which(int_gro_df$Intensity_growth != 0), ]$Intensity_growth)

# Set time steps
run_length <- 2000

# ASK: WHAT GENERATES THE EXPECTED LEVEL OF AGGREGATION?

# Create array to track status of hosts
status_ar <- array(NA, dim = c(nhst, run_length))

# Set all hosts to be uninfected
status_ar[, 1] <- 0

# Randomly select host to begin infection
status_ar[sample.int(nhst, 1), 1] <- 1

# Set probability of replacement - assume 1 every x years
p_replace <- 1/50

# Set transmission constant to define reciprocal relationship with distance - will be tweaked later
# Use from previous paper?
trnm_k <- 10

# Plot probability of transmission as an inverse function of distance - like Lavorel et al. 
# distances <- seq(0, max(adj_mat), length.out = 1000)
# p_trnms <- trnm_k / distances
# plot(x = distances, y = p_trnms)

# Set carrying capacity
cc <- 100

# Set seed
set.seed(16)

# Calculate mean height
mean_ht <- mean(pot_hst_df$Height, na.rm = TRUE)

# Calculate maximum height
max_ht <- max(pot_hst_df$Height, na.rm = TRUE)

# Simulate spread and growth in for loop
for(time in 1:(run_length - 1)){
  # Set next year intensity as this year's to start
  status_ar[, time + 1] <- status_ar[, time]
  # Generate list of random host numbers
  hsts <- sample(nhst)
  # Run through hosts randomly
  for(hst in hsts){
    # Get host height
    hst_ht <- as.numeric(pot_hst_df[hst, 4])
    
    # Check if mistletoe present on host or not
    if(status_ar[hst, time] != 0){
      
      # Change in intensity according to intensity growth rate 
      # Probability of change - otherwise will stay the same
      if(runif(1, 0, 1) <= 1 - p_no_change){
        # Given change, probability of change - *e^IGR - round to nearest integer
        gro_rt <- rnorm(1, mean = mean_igr, sd = sd_igr)
        status_ar[hst, time + 1] <- round(status_ar[hst, time] * exp(gro_rt), digits = 0)
      }
    } else {
      
      # Transmission - reciprocal function of distance - constant with respect to height for now
      # Designate infected hosts
      infecteds <- which(status_ar[, time] > 0)
      # Create list for distances
      p_trnms_ls <- list()

      # Run through infected hosts
      for(inf in infecteds){
        # Calculate distance to target host
        # Designate host index
        # index <- infecteds[i]
        # hst_lon <- pot_hst_df[index, 2]
        # hst_lat <- pot_hst_df[index, 3]
        # tar_lon <- pot_hst_df[hst, 2]
        # tar_lat <- pot_hst_df[hst, 3]
        # distance <- sqrt((hst_lon - tar_lon) ^ 2 + (hst_lat - tar_lat) ^ 2)
        
        p_trnms_ls[i] <- trnm_k / adj_mat[hst, inf] * hst_ht / mean_ht * log(status_ar[inf, time])
      }
      # Sum probabilities of transmission
      p_trnm_sum <- sum(unlist(p_trnms_ls))
      if(runif(1, 0, 1) <= p_trnm_sum){
        # Generate number of new recruits transmitted - could be Poisson, assume only 1 for now, trnm of a berry and successful establishment
        status_ar[hst, time + 1] <- 1 
      }
      
    }
    # DISTINGUISH SEED RAIN and MAKE NEW RECRUITS PROPORTIONAL TO FRUITING IN TREES
    # KEEP TACK OF NEW RECRUITS IN SEPARATE DATA FRAME
    
    # If reached over carrying capacity set to cc
    cc_ht <- round(cc * hst_ht / max_ht)
    if(status_ar[hst, time] >= cc_ht){
      status_ar[hst, time + 1] <- cc_ht
    }
    
    # Probability of replacement - reset parasite load to 0 with set probability
    if(runif(1, 0, 1) <= p_replace){
      status_ar[hst, time + 1] <- 0
    }

  }
  
  # Print outputs each time step
  n_infected <- sum(status_ar[, time] != 0)
  cat("The number of infected hosts at t =", time, "is", n_infected)
}

# Plot output - prevalence and index of dispersion for infected hosts
times <- 1:run_length
infected_ls <- list()
iod_ls <- list()
for(i in 1:run_length){
  infected_ls[i] <- sum(status_ar[, i] != 0) / nhst
  iod_ls[i] <- var(status_ar[which(status_ar[,i] != 0), i]) / mean(status_ar[which(status_ar[,i] != 0), i])
}

plot(times, unlist(infected_ls), type = "l", xlab = "Time", ylab = "Prevalence", ylim = c(0, 1))
abline(h = mean_prev, lty = 2)
plot(times, unlist(iod_ls), type = "l", xlab = "Time", ylab = "σ^2/μ (infecteds only)", ylim = c(0, 65))
abline(h = mean_iod, lty = 2)



# Set parameters ----------------------------------------------------------
# Set turnover rate - set to be constant, calculate from felling rate in dataset
tno_rt <- 
  
  # Set initial rate of new infections/transmission probability - number of trees infested in t+1 not infested in t
  trm_rt <- 
  # Probability of seed being picked up from tree A, landing on tree B, germinating and establishing new infection
  # This parameter will be tweaked later for best fit
  # Use value estimated from IPM paper for now
  
  # Set parasite infestation growth rate - mean increase in parasite load, accounting for missed individuals
  gro_rt  <-
  # This parameter will later be dependent on host species 
  
  # Set maximum intensity - from data - maximum observed * 110%
  max_it  <- 
  
  # Homogeneous model
  params <- c(tno_rt, # host turnover rate
              trm_rt, # transmission probability (probability of new infection) trm_probab = trm_rt / distance
              gro_rt, # rate of parasite infestation growth
              
  )



# Run model ---------------------------------------------------------------
# Set initial conditions

# Each time-step

# Randomly assign intensities for trees with "too many to count"

# Output number of parasites on each host

# Plot mean, variation and k parameter
