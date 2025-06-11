# Loading libraries
library(sp)
library(raster)
library(sf)
library(tidyverse)
library(fasterize)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)

# Import datasets
occ <- read_csv("data/com_data.csv")

# Get country polygons (scale = 'medium' gives a good balance of detail)
bd <- ne_countries(scale = "large", country = "Bangladesh", returnclass = "sf")
crs <- "+proj=cea +lon_0=0 +lat_ts=30 +datum=WGS84 +units=m +no_defs"
bd_proj <- st_transform(bd, crs = crs)
bd_sp <- as(bd_proj, "Spatial")  # For raster masking

# Use bbox of Bangladesh to create raster template
bd_bbox <- st_bbox(bd_sp)
bd_ext <- extent(c(bd_bbox["xmin"], bd_bbox["xmax"], bd_bbox["ymin"], bd_bbox["ymax"]))
bd_ras_template <- raster(bd_ext, res = 1000, crs = crs)

# Species list
species <- unique(occ$species)

# species <- c("Mikania micrantha", "Carassius auratus", "Argemone mexicana") # Problematic species
# i <- c("Mikania micrantha")
# j <- "Overall"

# Empty output dataframe
df <- data.frame()
df_rem <- data.frame()

for (i in species) try({
  print(i)
  speciesname <- gsub(" ", "_", i)
  
  sp_data_filt <- occ %>% filter(species == i)
  
  # Now generate Overall entry for this species
  sp_data_ov <- sp_data_filt %>%
    mutate(source = "Overall")
  
  # Combine Facebook, GBIF, and Overall for this species only
  sp_data_all <- bind_rows(sp_data_filt, sp_data_ov)
  
  # Get sources for this species
  data_sources <- unique(sp_data_all$source)

  for (j in data_sources) {
      
    sp_data_source <- sp_data_all %>% 
      filter(source == j)
    
    cat("Processing:", i, "| Source:", j, "| Records:", nrow(sp_data_source), "\n")
    
      if (NROW(sp_data_source) > 2) {
        
        # Convert to sf and project
        occ_sf <- st_as_sf(sp_data_source, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
        occ_sf_proj <- st_transform(occ_sf, crs = crs)
        
        # Convex hull
        multipoint <- st_combine(occ_sf_proj)
        hull <- st_convex_hull(multipoint)
        hull_sf <- st_sf(geometry = hull)
        
        # Rasterise onto Bangladesh-sized template
        hull_ras <- fasterize(hull_sf, bd_ras_template, field = NULL, background = NA)
        
        # Mask to Bangladesh
        hull_ras_masked <- mask(hull_ras, bd_sp)
        
        if (sum(!is.na(values(hull_ras_masked))) == 0) {
          cat("No valid cells for", i, "-", j, "\n")
          next
        }
        
        # Save raster
        raster_file <- paste0("output/range/", speciesname, "_", j, "_CH.tif")
        writeRaster(hull_ras_masked,
                    filename = raster_file,
                    NAflag = -9999,
                    overwrite = TRUE)
        
        if (!file.exists(raster_file)) {
          cat("Failed to write raster for", i, "-", j, "\n")
        }
        
        # Area in km² (1 cell = 1 km²)
        area_km2 <- cellStats(raster::Which(hull_ras_masked > 0), "sum")
        
        # Append to output dataframe
        df <- rbind(df_rem,
                    data.frame(
                      species = i,
                      source = j,
                      area_km2 = area_km2
                    ))
      }
    }
}, silent = FALSE)

# Pivot wider
df_wide <- df %>% 
  pivot_wider(names_from = source, values_from = area_km2) %>% 
  mutate(across(where(is.numeric), ~replace_na(., 0)),
         ov_contrib = (Overall - GBIF),
         fb_higher = (Facebook - GBIF))

# Export output
write_csv(df_wide, "output/ch_summary.csv")

#################################
## Changes in habitat area
# Importing datasets
ch <- read_csv("output/ch_summary.csv")
occ <- read_csv("data/com_data.csv")

# Extracting taxonomic information
taxa <- occ %>% 
  dplyr::select(taxa, species) %>% 
  dplyr::distinct(taxa, species, .keep_all = TRUE)

# Merging with the CH data
ch_taxa <- dplyr::left_join(ch, taxa, by = "species")

# Create an “Overall” row
ch_taxa_ov <- ch_taxa %>%
  mutate(taxa = "Overall")

# Combine with your original dataset
ch_taxa_com <- bind_rows(ch_taxa, ch_taxa_ov)

# Categorise the contributions into three bins
ch_taxa_com <- ch_taxa_com %>%
  mutate(contrib_group = case_when(
    ov_contrib == 0 ~ "GBIF",
    ov_contrib > 0 ~ "Facebook"
  ))

# Exporting output
write_csv(ch_taxa_com, "output/ch_taxa_com.csv")

ch_taxa_com$taxa <- factor(
  ch_taxa_com$taxa,
  levels = c("Birds", "Fishes", "Insects", "Molluscs", "Plants", "Overall")
)

ch_taxa_com <- ch_taxa_com %>%
  mutate(contrib_group = factor(contrib_group, levels = c("GBIF", "Facebook")))

# Plot
ggplot(ch_taxa_com, aes(x = ov_contrib, fill = contrib_group, color = contrib_group)) +
  geom_histogram(binwidth = 5000, position = "dodge") +  # Try 10000 or 2000 if 5000 too wide/narrow
  facet_wrap(~taxa, scales = "fixed") +  # Free x-scale helps with different taxa ranges
  scale_fill_manual(values = c(
    "Facebook" = "tomato",
    "GBIF" = "steelblue"
  )) +
  scale_color_manual(values = c(
    "Facebook" = "tomato",
    "GBIF" = "steelblue"
  )) +
  xlab("Change in range size") +
  ylab("Number of species") +
  theme_classic() +
  theme(legend.title = element_blank(),
        strip.text = element_text(margin = margin(b = 5)),   # adds space below strip text
        panel.spacing = unit(1.5, "lines"))                   # adds space between facets)

# Exporting image
ggsave("output/figures/ch_diff.png")

# Summary by taxa
ch_taxa_com_sum <- ch_taxa_com %>% 
  group_by(taxa) %>% 
  summarise(avg_change = mean(ov_contrib))

# Exporting output
write_csv(ch_taxa_com_sum, "output/ch_taxa_com_sum.csv")
