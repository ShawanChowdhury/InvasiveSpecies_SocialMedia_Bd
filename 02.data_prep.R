# Load libraries
library(tidyverse)
library(ggplot2)
library(countrycode)
library(CoordinateCleaner)
library(tidyverse)
library(rworldmap)
library(raster)
library(terra)
library(dismo)
library(geodata)
library(sf)
library(wdpar)

##############################
# Cleaning GBIF data
# Reading GBIF data
gbif_data <- data.table::fread("data/gbif/gbif.csv")

# Removing blank cells
gbif_data <- gbif_data[!(is.na(gbif_data$species) | gbif_data$species == ""),]
gbif_data <- gbif_data[!(is.na(gbif_data$decimalLongitude) | gbif_data$decimalLongitude == ""),]
gbif_data <- gbif_data[!(is.na(gbif_data$decimalLatitude) | gbif_data$decimalLatitude == ""),]

# Removing duplicated records
gbif_data_dedup <- gbif_data[!duplicated(gbif_data),]

# Cleaning memory
rm(gbif_data)

write_csv(gbif_data_dedup, "data/cleanedRecords_GBIF.csv")

##############################
# Combining datasets
# Import datasets
gbif_data <- read_csv("data/cleanedRecords_GBIF.csv")
fb_data <- read_csv("data/occ_fb.csv")

# Facebook data is from 2012, so we are removing GBIF data before 2012
gbif_data <- gbif_data %>% 
  filter(year > 2011)

# Subsetting by columns
gbif_data <- gbif_data %>% 
  select(class, order, family, species, decimalLongitude, decimalLatitude, day, month, year)

fb_data <- fb_data %>% 
  select(class, order, family, species, decimalLongitude, decimalLatitude, day, month, year)

# Adding data source
gbif_data$source <- "GBIF"
fb_data$source <- "Facebook"

combined_data <- rbind(fb_data, gbif_data)

write_csv(combined_data, "data/com_data.csv")

rm(fb_data, gbif_data)

##############################
# Import data
occ <- read_csv("data/com_data.csv")

# Data by source and taxa
occ_taxa_source <- occ %>% 
  group_by(source, taxa) %>% 
  summarise(n = NROW(species))

# Data by source
occ_source <- occ %>% 
  group_by(source) %>% 
  summarise(n = NROW(species))

# Number of occurrence records by species
sp_occ <- occ %>% 
  group_by(taxa, species, source) %>% 
  summarise(n = NROW(species))

# Plot
ggplot(sp_occ, aes(taxa, n)) +
  geom_boxplot() +  # suppress default black outliers
  geom_jitter(aes(col = source), width = 0.1, alpha = 0.5) +
  theme_classic() +
  scale_color_manual(values = c("Facebook" = "tomato", "GBIF" = "deepskyblue")) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 10)
  ) +
  xlab("") +
  ylab("Number of occurrence records") + 
  scale_y_log10()

ggsave("output/figures/sp_bias.png")

# Sp occ wide
sp_occ_wide <- occ %>% 
  group_by(taxa, species, source) %>% 
  summarise(n = NROW(species)) %>% 
  pivot_wider(names_from = "source", values_from = "n") %>%
  mutate(
    dominant_source = case_when(
      is.na(Facebook) & is.na(GBIF) ~ NA_character_,
      is.na(Facebook) ~ "GBIF",
      is.na(GBIF) ~ "Facebook",
      Facebook > GBIF ~ "Facebook",
      GBIF > Facebook ~ "GBIF",
      Facebook == GBIF ~ "Equal"
    )
  )




##############################
# Cropping variables to bd extent
# Import datasets
elev <- rast("data/elev/elev.tif")
r <- rast("output/range/Acacia_mangium_Facebook_CH.tif") # Test species

# Bangladesh boundary
bd <- geodata::gadm('GADM', country='BGD', level=1)

# Reproject to match
elev_re <- project(elev, r, method = "bilinear")
bd_proj <- project(bd, elev_re)

# Cropping layer
elev_crop <- crop(elev_re, bd_proj)
elev_crop <- mask(elev_crop, bd_proj)

###########################
# BIO1 = Annual Mean Temperature (°C × 10)
# BIO12 = Annual Precipitation (mm)
# Import a sample species
r <- rast("output/range/Acacia_mangium_Facebook_CH.tif") # Test species

# Download only Bangladesh bioclimatic data
bio_bd <- worldclim_country(country = "BGD", var = "bio", res = 0.5 / 60, path = "data/")

# Extract layers of interest
temp_bd <- bio_bd[[1]] / 10         # BIO1: divide by 10 to convert to °C
rain_bd <- bio_bd[[12]]           # BIO12

# Download elevation for Bangladesh
elev_bd <- worldclim_country(country = "BGD", var = "elev", res = 0.5 / 60, path = "data/")

# Stack and rename
climate_bd <- c(temp_bd, rain_bd, elev_bd)
names(climate_bd) <- c("mean_temp", "annual_precip", "elevation")

# Bangladesh boundary
bd <- geodata::gadm('GADM', country='BGD', level=1)

# Reproject to match
climate_bd_re <- project(climate_bd, r, method = "bilinear")
bd_proj <- project(bd, climate_bd_re)

# Cropping layer
climate_bd_re_crop <- crop(climate_bd_re, bd_proj)
climate_bd_re_crop <- mask(climate_bd_re_crop, bd_proj)

# # Plot and check
# plot(climate_bd_re_crop[[1]])
# plot(r, add = TRUE)

# Export rasters
writeRaster(climate_bd_re_crop, "data/climate_bd.tif", overwrite = TRUE, filetype = "GTiff")

########################
# HFP data [source: https://wcshumanfootprint.org/data-access; 2020-01-01]
# Import datasets
hfp <- rast("data/hfp.tif")
r <- rast("output/range/Acacia_mangium_Facebook_CH.tif") # Test species

# Bangladesh boundary
bd <- geodata::gadm('GADM', country='BGD', level=1)

# Reproject to match
hfp_re <- project(hfp, r, method = "bilinear")
bd_proj <- project(bd, hfp_re)

# Cropping layer
hfp_crop <- crop(hfp_re, bd_proj)
hfp_crop <- mask(hfp_crop, bd_proj)

# Export rasters
writeRaster(hfp_crop, "data/hfp_bd.tif", overwrite = TRUE)

########################
# Built areas [data source: https://human-settlement.emergency.copernicus.eu/download.php]
# Import datasets
built_areas <- rast("data/built_areas/built_areas.tif")
r <- rast("output/range/Acacia_mangium_Facebook_CH.tif") # Test species

# Bangladesh boundary
bd <- geodata::gadm('GADM', country='BGD', level=1)

# Reproject to match
built_areas_re <- project(built_areas, r, method = "bilinear")
bd_proj <- project(bd, built_areas_re)

# Cropping layer
built_areas_crop <- crop(built_areas_re, bd_proj)
built_areas_crop <- mask(built_areas_crop, bd_proj)

# Export rasters
writeRaster(built_areas_crop, "data/built_areas_bd.tif", overwrite = TRUE)

############################
# Species summary
# Import datasets
occ <- read_csv("data/com_data.csv")

sp <- occ %>% 
  group_by(species, taxa, source) %>% 
  summarise(n = NROW(species)) %>% 
  pivot_wider(names_from = "source",
              values_from = "n")

