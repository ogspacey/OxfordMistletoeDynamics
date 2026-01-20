## CSV BINDER ##
# This code combines all .csv files from the folder of choice, and outputs this as one
# large .csv, which will be named after the Tree Code + "_data".

library(tidyverse)
library(rstudioapi)

# -1- Ask user to pick the folder ----------------------------------------------------------
wd <- selectDirectory(caption = "Select folder containing chosen CSVs: ")
setwd(wd)


# -2- Read and combine all CSVs ------------------------------------------------------------
csv_files <- list.files(path = wd, pattern = "\\.csv$", full.names = TRUE) %>%
  str_sort(numeric = TRUE)

# Creates a new dataframe which combines rows from each CSV based on the column
df <- lapply(csv_files, function(file) {
  read_csv(file, na = c("#N/A"), show_col_types = FALSE)   # treat #N/A as NA
}) %>%
  bind_rows()


# -3- Output the filename using the tree code prefix ---------------------------------------
prefix <- str_sub(basename(csv_files[1]), 1, 3) 
output_path <- file.path(wd, str_c(prefix, "_data.csv"))


# -4- Save combined dataframe --------------------------------------------------------------
write_csv(df, output_path)

cat("New CSV saved to: ", output_path, "\n")
