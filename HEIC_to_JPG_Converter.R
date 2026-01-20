## HEIC TO JPG CONVERTER ##
# This code allows the user to convert a file containing .HEIC images to .JPG images,
# and outputs a folder containing .JPG images next to the .HEIC images.

library(magick)
library(rstudioapi)   # provides function to let user pick directory

# -1- Ask user to pick the folder  ----------------------------------------------------
wd <- selectDirectory(caption = "Select the folder containing your chosen images:")
setwd(wd)

# -2- Get all HEIC files --------------------------------------------------------------
heic_files <- list.files(pattern = "\\.HEIC$", ignore.case = TRUE)


# -3- Define and create output folder -------------------------------------------------
output_folder <- "JPG_images"
if (!dir.exists(output_folder)) dir.create(output_folder)


# -4- Function to convert one file ----------------------------------------------------
convert_to_jpg <- function(file) {
  img <- image_read(file)                        # Read HEIC
  img_jpg <- image_convert(img, format = "jpg")  # Convert to JPG
  jpg_file <- file.path(output_folder, paste0(tools::file_path_sans_ext(file), ".jpg"))
  image_write(img_jpg, path = jpg_file)          # Write JPG
  cat("Converted:", file, "->", jpg_file, "\n")
}


# -5- Apply function to all HEIC files ------------------------------------------------
lapply(heic_files, convert_to_jpg)

cat("All files converted!\n")
