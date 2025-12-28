### Oxford IBM
### 2025.12.28
### Oliver G. Spacey
### Code to generate 

# Pre-amble ---------------------------------------------------------------
# Clear environment
rm(list = ls())

# Load packages 
library(igraph)    # Network construction
library(tidyverse) # Data wrangling

# Set working directory as required

# Load data --------------------------------------------------------------
# Load raw census data
raw_cen_df <- read.csv("Oxford_Trees_Long.csv")

# Load host coordinate data
hst_coo_df <- read.csv("host_trees_bng_with_coords.csv") %>%
              rename(TreeID = Name,
                     lon = easting_bng,
                     lat = northing_bng)%>%
              mutate(Index = 1:nrow(read.csv("host_trees_bng_with_coords.csv")))

# Join raw data and coordinates, and add indexing column
oxf_hst_df <- full_join(raw_cen_df, hst_coo_df) 

# Add in data for non-hosts of host species - randomly select

# Generate host community ------------------------------------------------
# Get number of hosts
n <- nrow(hst_coo_df)

## Weight host network by distance
# Create empty adjacency matrix
adj_mat <- matrix(NA, nrow = n, ncol = n)

# Calculate pythagorean distances between hosts and store in adjacency matrix
for(i in 1:n){
  i_lon <- hst_coo_df[i,2]
  i_lat <- hst_coo_df[i,3]
  for(j in 1:n){
    j_lon <- hst_coo_df[j,2]
    j_lat <- hst_coo_df[j,3]
    adj_mat[i,j] <- sqrt((i_lon - j_lon) ^ 2 + (i_lat - j_lat) ^ 2)
  }
}

# Construct network from adjacency matrix
hst_ntw <- graph_from_adjacency_matrix(adjmatrix = adj_mat,
                                       mode = "undirected",
                                       weighted = TRUE)

# Plot host network
plot(hst_ntw)

# Set parameters ----------------------------------------------------------
# GO FROM HERE - GET THESE PARAMETERS
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
  
# Set maximum intensity - from data - maximum observed * 10%
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


