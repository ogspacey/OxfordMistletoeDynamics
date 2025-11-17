### Oxford Tree Wide
### 2025.11.17
### Oliver G. Spacey

# Script to extract mean height and DBH from long file and output wide data frame with measurements for each tree

# Clear environment
rm(list = ls())

# Set working directory

# Load required packages
library(tidyverse)

# Load data ---------------------------------------------------------------
# Load entire dataset in long format
oxford_trees_df <- read.csv("Oxford_Trees_Long.csv")

# Wrangle data ------------------------------------------------------------
# Select relevant columns - TreeID, Height and DBH
ID_ht_df <- select(oxford_trees_df, TreeID, Height_Tri, DBH)

# Convert heights and DBHs to numeric variables
ID_ht_df$Height_Tri <- as.numeric(ID_ht_df$Height_Tri)
ID_ht_df$DBH        <- as.numeric(ID_ht_df$DBH)

# Calculate mean height and DBH for each tree
means_df <- ID_ht_df %>%
  group_by(TreeID) %>%
  summarise(Height = mean(Height_Tri, na.rm = TRUE),
            DBH = mean(DBH, na.rm = TRUE))

# Export data -------------------------------------------------------------
# Write .csv 
write.csv(means_df, file = "Oxford_Trees_Wide.csv")
