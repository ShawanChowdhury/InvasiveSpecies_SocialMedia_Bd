# 1. Load libraries
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(raster)
library(sp)
library(tmap)
library(rnaturalearth)
library(rnaturalearthdata)

# Import datasets
occ <- read_csv("data/com_data.csv")

# Creating the overall source
occ_ov <- occ %>% 
  mutate(source = "Overall")

occ_all <- rbind(occ, occ_ov)

# Separating GBIF data
gbif_data <- occ_all %>% 
  filter(source == "GBIF")
fb_data <- occ_all %>% 
  filter(source == "Facebook")
ov_data <- occ_all %>% 
  filter(source == "Overall")

# Set CRS (coordinate reference system)
crs_proj <- CRS("+proj=utm +zone=46 +datum=WGS84 +units=km +no_defs") # UTM for Bangladesh

# Convert to spatial points
gbif_sf <- st_as_sf(gbif_data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
fb_sf <- st_as_sf(fb_data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
com_sf   <- st_as_sf(ov_data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# Reproject all to UTM
gbif_utm     <- st_transform(gbif_sf, crs = crs_proj)
fb_utm     <- st_transform(fb_sf, crs = crs_proj)
com_utm <- st_transform(com_sf, crs = crs_proj)

# Create 5x5 km grid for Bangladesh
# Get bounding box covering all points
bbox <- st_bbox(com_utm)

# Create grid
grid <- st_make_grid(
  st_as_sfc(bbox),
  cellsize = c(1, 1),
  square = TRUE
) %>% st_sf(grid_id = 1:length(.))

# Spatial join: map points to grid cells
gbif_grid     <- st_join(gbif_utm, grid, join = st_within) %>% filter(!is.na(grid_id))
fb_grid     <- st_join(fb_utm, grid, join = st_within) %>% filter(!is.na(grid_id))
com_grid <- st_join(com_utm, grid, join = st_within) %>% filter(!is.na(grid_id))

# Count unique grid cells
gbif_cells     <- unique(gbif_grid$grid_id) # 1129
fb_cells     <- unique(fb_grid$grid_id) # 245
com_cells <- unique(com_grid$grid_id) # 1129

# Calculate new grid cells added by Facebook data
new_cells <- setdiff(com_cells, gbif_cells) # 182

# Optional: visualise
grid$source <- ifelse(grid$grid_id %in% new_cells, "New from FB", 
                      ifelse(grid$grid_id %in% gbif_cells, "GBIF only", NA))

##############
# Present in both
grid$source[grid$grid_id %in% intersect(gbif_cells, com_cells)] <- "Both"

# GBIF-only
grid$source[grid$grid_id %in% setdiff(gbif_cells, fb_cells)] <- "GBIF only"

# Facebook-only
grid$source[grid$grid_id %in% setdiff(com_cells, gbif_cells)] <- "Facebook only"

# Ensure 'source' is clean and not a factor
grid$source <- as.character(grid$source)
grid$source[is.na(grid$source)] <- "Absent"

# Get Bangladesh polygon
bangladesh <- ne_countries(scale = "large", country = "Bangladesh", returnclass = "sf")
bangladesh_utm <- st_transform(bangladesh, crs = st_crs(grid))  # match CRS

grid$source <- factor(grid$source, levels = c("Both", "Facebook only", "GBIF only", "Absent"))

ggplot() +
  geom_sf(data = grid, aes(fill = source), color = NA) +  # fill without borders
  geom_sf(data = bangladesh_utm, fill = NA, color = "grey80", size = 0.4) +  # outline on top
  scale_fill_manual(
    values = c(
      "Both" = "grey",
      "Facebook only" = "tomato",
      "GBIF only" = "deepskyblue3",
      "Absent" = "white"
    ),
    name = ""
  ) +
  coord_sf() +
  theme_void() +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave("output/figures/grid_all_taxa.png")

#################################
# Import datasets
occ <- read_csv("data/com_data.csv")

# Creating the overall source
occ_ov <- occ %>% 
  mutate(source = "Overall")

occ_all <- rbind(occ, occ_ov)

# Set CRS (coordinate reference system)
crs_proj <- CRS("+proj=utm +zone=46 +datum=WGS84 +units=km +no_defs") # UTM for Bangladesh

# Convert to spatial points
all_sf   <- st_as_sf(occ_all, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# Reproject all to UTM
all_utm <- st_transform(all_sf, crs = crs_proj)

# Create 1x1 km grid for Bangladesh
# Get bounding box covering all points
bbox <- st_bbox(all_utm)

# Create grid
grid <- st_make_grid(
  st_as_sfc(bbox),
  cellsize = c(1, 1),
  square = TRUE
) %>% st_sf(grid_id = 1:length(.))

# Spatial join: map points to grid cells
all_grid <- st_join(all_utm, grid, join = st_within) %>% filter(!is.na(grid_id))

# Create a grid-taxonomy-source presence table
presence_tbl <- all_grid %>%
  st_drop_geometry() %>%
  distinct(grid_id, taxa, source) %>%
  mutate(present = TRUE) %>%
  pivot_wider(names_from = source, values_from = present, values_fill = FALSE)

# Assign presence category
presence_tbl <- presence_tbl %>%
  mutate(category = case_when(
    Facebook & GBIF ~ "Both",
    Facebook & !GBIF ~ "Facebook only",
    !Facebook & GBIF ~ "GBIF only",
    !Facebook & !GBIF ~ "Absent"
  ))

# Join with grid geometry
grid_taxa_sf <- presence_tbl %>%
  left_join(grid, by = "grid_id") %>%
  st_as_sf()

# Create overall presence by grid + source
overall_tbl <- all_grid %>%
  st_drop_geometry() %>%
  distinct(grid_id, source) %>%
  mutate(present = TRUE) %>%
  pivot_wider(names_from = source, values_from = present, values_fill = FALSE) %>%
  mutate(category = case_when(
    Facebook & GBIF ~ "Both",
    Facebook & !GBIF ~ "Facebook only",
    !Facebook & GBIF ~ "GBIF only",
    !Facebook & !GBIF ~ "Absent"
  )) %>%
  mutate(taxa = "Overall")  # add 'Overall' as a fake taxon

# Join with grid to create spatial object
overall_sf <- overall_tbl %>%
  left_join(grid, by = "grid_id") %>%
  st_as_sf()

# Combine with your existing grid_taxa_sf
grid_taxa_all <- bind_rows(grid_taxa_sf, overall_sf)

# Determine order: existing taxa, then "Overall"
ordered_taxa <- c(
  sort(unique(grid_taxa_sf$taxa)),  # this gives current taxon names in alphabetical order
  "Overall"
)

# Apply the order to the combined dataset
grid_taxa_all$taxa <- factor(grid_taxa_all$taxa, levels = ordered_taxa)

grid_taxa_all_filtered <- grid_taxa_all %>%
  filter(category != "Absent")

# Get Bangladesh polygon
bangladesh <- ne_countries(scale = "large", country = "Bangladesh", returnclass = "sf")
bangladesh_utm <- st_transform(bangladesh, crs = st_crs(grid))  # match CRS

ggplot() +
  geom_sf(data = grid_taxa_all_filtered, aes(fill = category), color = NA) +
  geom_sf(data = bangladesh_utm, fill = NA, color = "grey80", size = 0.3) +
  facet_wrap(~taxa) +
  scale_fill_manual(
    values = c("Both" = "grey",
               "Facebook only" = "tomato",
               "GBIF only" = "deepskyblue3"),
    na.translate = FALSE  # hides any remaining NA in legend
  ) +
  coord_sf() +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    strip.text = element_text(size = 12)
  )

ggsave("output/figures/grid_taxa.png")

# Pie Chart
# Count number of grid cells per category per taxa
# Prepare data: counts and proportions
donut_data <- grid_taxa_all_filtered %>%
  st_drop_geometry() %>%
  group_by(taxa, category) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(taxa) %>%
  mutate(percentage = count / sum(count) * 100)

# Export table
write_csv(donut_data, "output/donut_data.csv")

# Donut chart
ggplot(donut_data, aes(x = 2, y = percentage, fill = category)) +
  geom_bar(stat = "identity", width = 1, colour = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +  # creates the donut hole by pushing bars outward
  facet_wrap(~taxa) +
  scale_fill_manual(
    values = c("Both" = "grey",
               "Facebook only" = "tomato",
               "GBIF only" = "deepskyblue3")
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 12),
    panel.spacing = unit(1, "lines")
  ) +
  labs(fill = "Data source")

# Export figure
ggsave("output/figures/grid_taxa_donut.png")

# grid_taxa_all_filtered2 <- grid_taxa_all %>%
#   filter(taxa != "Overall")
# 
# ggplot() +
#   geom_sf(data = grid_taxa_all_filtered2, aes(fill = taxa), color = NA) +
#   geom_sf(data = bangladesh_utm, fill = NA, color = "grey80", size = 0.3) +
#   facet_wrap(~category) +
#   scale_fill_viridis_d(option = "A") +
#   coord_sf() +
#   theme_void() +
#   theme(
#     legend.position = "right",
#     legend.title = element_blank(),
#     strip.text = element_text(size = 12)
#   )
# 
# ggsave("output/figures/grid_taxa.png")

# Calculating summary by taxa
grid_taxa_sum <- grid_taxa_all %>% 
  group_by(taxa, category) %>% 
  summarise(n = NROW(grid_id)) %>% 
  mutate(total_cell = sum(n), prop_cell = (n/total_cell)*100)

# Removing geometry
grid_taxa_sum$. <- NULL 

# Exporting output
write_csv(grid_taxa_sum, "output/grid_taxa_sum.csv")

#################################################
# Species-wise grid details
# Create a grid-taxonomy-source presence table
presence_tbl <- all_grid %>%
  st_drop_geometry() %>%
  distinct(grid_id, taxa, species, source) %>%
  mutate(present = TRUE) %>%
  pivot_wider(names_from = source, values_from = present, values_fill = FALSE)

# Assign presence category
presence_tbl <- presence_tbl %>%
  mutate(category = case_when(
    Facebook & GBIF ~ "Both",
    Facebook & !GBIF ~ "Facebook only",
    !Facebook & GBIF ~ "GBIF only",
    !Facebook & !GBIF ~ "Absent"
  ))

# Join with grid geometry
grid_taxa_sf <- presence_tbl %>%
  left_join(grid, by = "grid_id") %>%
  st_as_sf()

# Create overall presence by grid + source
overall_tbl <- all_grid %>%
  st_drop_geometry() %>%
  distinct(grid_id, source) %>%
  mutate(present = TRUE) %>%
  pivot_wider(names_from = source, values_from = present, values_fill = FALSE) %>%
  mutate(category = case_when(
    Facebook & GBIF ~ "Both",
    Facebook & !GBIF ~ "Facebook only",
    !Facebook & GBIF ~ "GBIF only",
    !Facebook & !GBIF ~ "Absent"
  )) %>%
  mutate(taxa = "Overall")  # add 'Overall' as a fake taxon

# Join with grid to create spatial object
overall_sf <- overall_tbl %>%
  left_join(grid, by = "grid_id") %>%
  st_as_sf()

# Combine with your existing grid_taxa_sf
grid_taxa_all <- bind_rows(grid_taxa_sf, overall_sf)

# Calculating summary by species
grid_sp_sum <- grid_taxa_all %>% 
  group_by(taxa, species, category) %>% 
  summarise(n = NROW(grid_id)) %>% 
  mutate(total_cell = sum(n), prop_cell = (n/total_cell)*100)

# Removing geometry
grid_sp_sum$. <- NULL 

# Exporting output
write_csv(grid_sp_sum, "output/grid_sp_sum.csv")
