## HEIC/PNG TO JPG CONVERTER ##
# This code allows the user to convert a file containing .HEIC images to .JPG images,
# and outputs a folder containing .JPG images next to the .HEIC images.

library(magick)
library(svDialogs)   # provides function to construct dialog boxes for GUI


# -1- Ask user to pick the folder --------------------------------------------------------------
wd <- dlg_dir(title = "Select the folder containing your chosen images")$res
setwd(wd)


# -2- Allow user to select HEIC or PNG to convert ----------------------------------------------
png_or_heic <- dlg_list(c("HEIC", "PNG"), title = "Choose format to convert to JPG")$res


# -2.1- Create a function which combines the selected images and creates an output folder ------
convert_all <- function(ext, output_folder) {
  files <- list.files(pattern = paste0("\\.", ext, "$"), ignore.case = TRUE)
  
  
  # Create an output folder with the converted images if it does not exist
  if (!dir.exists(output_folder)) dir.create(output_folder)
  
  
  # -2.2- Create a function which converts an image to JPG -------------------------------------
  convert_to_jpg <- function(file) {
    img <- image_read(file)
    img_jpg <- image_convert(img, format = "jpg")
    jpg_file <- file.path(output_folder, paste0(tools::file_path_sans_ext(file), ".jpg"))
    image_write(img_jpg, path = jpg_file)
    cat("Converted:", file, "->", jpg_file, "\n")
  }
  
  
  # -2.3- Apply the convert_to_jpg function to each image in 'files' ---------------------------
  lapply(files, convert_to_jpg)
  cat("All", toupper(ext), "files converted!\n")
}


# -3- Apply the convert_all function based from file ext of choice made by user ----------------
if (png_or_heic == "HEIC") {
  convert_all("HEIC", "HEIC_to_JPG")

} else if (png_or_heic == "PNG") {
  convert_all("PNG", "PNG_to_JPG")
  
} else
  cat("Invalid input. Please run the code again.\n")
