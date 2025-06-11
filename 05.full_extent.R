# Load libraries
library(raster)

# Raster with full extent
r <- raster("data/hfp_bd.tif")
r[] <- 0

# List of rasters
tifs <- list.files(path = "output/range/", pattern = ".tif", recursive = TRUE, full.names = TRUE)

for (i in tifs) {
  print(i)
  
  # Changing file name
  file_name <- gsub("output/range/", "", i)
  file_name <- gsub(".tif", "", file_name)
  
  # Import species raster
  ras <- raster(i)
  
  # Reproject to match full extent raster
  proj_ras <- projectRaster(ras, r, method = "ngb")  # nearest neighbour for binary
  
  # Replace NA (outside convex hull) with 0 (absence)
  proj_ras[is.na(proj_ras)] <- 0
  
  # Optional: ensure it's binary
  proj_ras[proj_ras > 0] <- 1
  
  # Exporting raster
  writeRaster(proj_ras, paste0("output/range_full/", file_name, ".tif"), 
              NAflag=-9999, overwrite = TRUE)
  
}
