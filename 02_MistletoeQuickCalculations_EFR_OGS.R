# EXPLORATORY ANALYSIS

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)

# Set working directory
getwd()

# 1 Import data
df <- read_csv("MASTER_data_wide.csv")

# 2. Structure checks
glimpse(df)

# Number of trees
n_distinct(df$TreeID)

# Number of years
n_distinct(df$Year)

# Rows per year
df %>% count(Year)


# 3. Count mistletoes per tree
df <- df %>%
  mutate(
    n_mistletoe = rowSums(!is.na(select(., starts_with("M") & ends_with("_Area"))))
  )

hist(df$n_mistletoe)


# 4. Total mistletoe load per tree
df <- df %>%
  mutate(
    total_mistletoe_area = rowSums(select(., ends_with("_Area") & starts_with("M")), na.rm = TRUE)
  )

hist(df$total_mistletoe_area)

# Scatterplot of number of mistletoe and total mistletoe area
ggplot(data = df, aes(x = n_mistletoe, y = total_mistletoe_area)) +
  geom_point() +
  theme_bw()

# 5. Tree size distributions
ggplot(df, aes(TreeArea)) +
  geom_histogram(bins = 40) +
  theme_minimal() +
  labs(title = "Distribution of Tree Area")

# 6. Mistletoe abundance distribution
ggplot(df, aes(n_mistletoe)) +
  geom_histogram(binwidth = 1) +
  theme_minimal() +
  labs(title = "Number of Mistletoes per Tree")

# 7. Mistletoe load vs tree size
ggplot(df, aes(TreeArea, total_mistletoe_area)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  theme_minimal() +
  labs(title = "Mistletoe Load vs Tree Size")

ggplot(df, aes(TreeArea, n_mistletoe, group = TreeID)) +
  geom_line(alpha = 0.2) +
  geom_point(alpha = 0.3) +
  theme_minimal() +
  labs(title = "Within-tree relationship: Size vs Mistletoe")

# 8. Relationship with crown
ggplot(df, aes(MaskArea, total_mistletoe_area)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  theme_minimal() +
  labs(title = "Mistletoe Load vs Crown Area")

# Calculate proportional intensity - area of mistletoe as proportion of crown area
df <- df %>%
  mutate(
    proportional_intensity = total_mistletoe_area / MaskArea
    )

# Histogram of proportional intensity
ggplot(df, aes(proportional_intensity)) +
  geom_histogram(binwidth = 0.01) +
  theme_minimal() +
  labs(title = "Proportional intensity")

# Plot discrete intensity against proportional intensity
ggplot(df, aes(n_mistletoe, proportional_intensity)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  theme_minimal() +
  labs(title = "Discrete intensity vs Proportional intensity")

# 9. Trends over time
df %>%
  group_by(Year) %>%
  summarise(
    mean_mistletoe = mean(n_mistletoe, na.rm = TRUE),
    mean_load = mean(total_mistletoe_area, na.rm = TRUE)
  ) %>%
  ggplot(aes(Year, mean_mistletoe, group = 1)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Mean Mistletoe Count Over Time")


df %>%
  group_by(Year) %>%
  summarise(
    mean_mistletoe = mean(n_mistletoe, na.rm = TRUE),
    sd_mistletoe = sd(n_mistletoe, na.rm = TRUE)
  ) %>%
  ggplot(aes(Year, mean_mistletoe, group = 1)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(
    ymin = mean_mistletoe - sd_mistletoe,
    ymax = mean_mistletoe + sd_mistletoe
  ), width = 0.1) +
  theme_minimal() +
  labs(title = "Average Mistletoe Count Over Time")









#Spaghetti
df_multi <- df %>%
  group_by(Tree) %>%
  filter(n() > 1) %>%
  ungroup()

ggplot(df_multi, aes(x = Year, y = n_mistletoe, group = Tree)) +
  geom_line(alpha = 0.2) +
  geom_point(alpha = 0.3) +
  theme_minimal() +
  labs(
    title = "Mistletoe Count Trajectories per Tree",
    y = "Number of mistletoes"
  )


df_change <- df %>%
  arrange(Tree, Year) %>%
  group_by(Tree) %>%
  mutate(
    delta_mistletoe = n_mistletoe - lag(n_mistletoe)
  ) %>%
  ungroup() %>%
  filter(!is.na(delta_mistletoe))

ggplot(df_change, aes(delta_mistletoe)) +
  geom_histogram(binwidth = 1) +
  theme_minimal() +
  labs(title = "Change in Mistletoe Count Between Years") + geom_vline(xintercept = 0, col= "red")



# 10. Height distribution of mistletoes
height_long <- df %>%
  select(TreeID, Year, starts_with("M") & ends_with("TreeHeightRelative")) %>%
  pivot_longer(
    cols = -c(TreeID, Year),
    values_to = "height_rel"
  ) %>%
  filter(!is.na(height_rel))

ggplot(height_long, aes(height_rel)) +
  geom_histogram(bins = 30) +
  theme_minimal() +
  labs(title = "Relative Height Distribution of Mistletoes")


# 11. Spatial distribution
xy_long <- df %>%
  select(TreeID, Year, matches("^M\\d+_(X|Y)$")) %>%
  
  pivot_longer(
    cols = -c(TreeID, Year),
    names_to = c("M_id", "coord"),
    names_pattern = "(M\\d+)_(X|Y)",
    values_to = "value"
  ) %>%
  
  pivot_wider(
    names_from = coord,
    values_from = value
  ) %>%
  
  filter(!is.na(X) & !is.na(Y))

ggplot(xy_long, aes(X, Y)) +
  geom_point(alpha = 0.3) +
  theme_minimal() +
  labs(title = "Spatial Distribution of Mistletoes")



height_long <- df %>%
  select(TreeID, Year, matches("^M\\d+_TreeHeightRelative$")) %>%
  
  pivot_longer(
    cols = -c(TreeID, Year),
    names_to = c("M_id"),
    names_pattern = "(M\\d+)_TreeHeightRelative",
    values_to = "height_rel"
  ) %>%
  
  filter(!is.na(height_rel))

ggplot(height_long, aes(height_rel)) +
  geom_histogram(bins = 30) +
  theme_minimal() +
  labs(
    title = "Distribution of Mistletoe Relative Height",
    x = "Relative height (0 = base, 1 = top)",
    y = "Count"
  )


mistletoe_long <- df %>%
  select(
    TreeID, Year,
    matches("^M\\d+_(TreeHeightRelative|Area)$")
  ) %>%
  
  pivot_longer(
    cols = -c(TreeID, Year),
    names_to = c("M_id", "variable"),
    names_pattern = "(M\\d+)_(.*)",
    values_to = "value"
  ) %>%
  
  pivot_wider(
    names_from = variable,
    values_from = value
  ) %>%
  
  filter(!is.na(TreeHeightRelative) & !is.na(Area))


ggplot(mistletoe_long, aes(
  x = 1,  # all trees stacked (we'll jitter)
  y = TreeHeightRelative,
  size = Area,
  colour = Area
)) +
  
  geom_jitter(
    width = 0.2,
    height = 0,
    alpha = 0.6
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    name = "Relative height (0 = base, 1 = top)"
  ) +
  
  scale_size_continuous(range = c(1, 6)) +
  
  scale_colour_viridis_c() +
  
  theme_minimal() +
  
  labs(
    title = "Vertical Distribution of Mistletoes on Trees",
    x = NULL,
    colour = "Area",
    size = "Area"
  ) +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# By year
ggplot(mistletoe_long, aes(
  x = 1,
  y = TreeHeightRelative,
  size = Area,
  colour = Area
)) +
  
  geom_jitter(width = 0.2, alpha = 0.6) +
  
  facet_wrap(~Year) +
  
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_viridis_c() +
  scale_size_continuous(range = c(1, 5)) +
  
  theme_minimal() +
  
  labs(
    title = "Mistletoe Vertical Distribution by Year",
    y = "Relative height"
  ) +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


# with density
ggplot(mistletoe_long, aes(
  x = factor(1),
  y = TreeHeightRelative
)) +
  
  # Half violin (density)
  geom_violin(
    fill = "grey80",
    alpha = 0.6,
    width = 0.5,
    trim = TRUE
  ) +
  
  # Points
  geom_jitter(aes(
    size = Area,
    colour = Area
  ),
  width = 0.15,
  alpha = 0.6
  ) +
  
  facet_wrap(~Year) +
  
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_viridis_c() +
  scale_size_continuous(range = c(1, 5)) +
  
  theme_minimal() +
  
  labs(
    title = "Mistletoe Vertical Distribution by Year",
    y = "Relative height"
  ) +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )





# 12. Data quality checks
# Trees with no mistletoes
df %>% filter(n_mistletoe == 0) #  17 trees

# Missing mask values
df %>% filter(is.na(MaskArea)) # 13 trees do not have crown measurements

# Extreme values
summary(df$TreeArea)
hist(df$TreeArea)
hist(df$total_mistletoe_area)



