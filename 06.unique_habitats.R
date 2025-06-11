# Load libraries
library(raster)
library(terra)
library(tidyverse)
library(stringr)
library(purrr)
library(broom)
library(dplyr)

raster_dir <- "output/range_full"
tifs <- list.files(raster_dir, pattern = "\\.tif$", full.names = TRUE)

# Importing predictor variables
built <- rast("data/built_areas_bd.tif")
hfp <- rast("data/hfp_bd.tif")
clim <- rast("data/climate_bd.tif") # temperature, rainfall, and elevation

# stacking rasters
vars <- c(built, hfp, clim)
names(vars)

# Import species
sp <- read_csv("output/ch_summary.csv")
sp_list <- unique(sp$species)

df <- data.frame()

for (i in sp_list) try({
  print(i)

  # i <- "Carassius auratus"

  # Changing name structure
  sp_name <- gsub(" ", "_", i)

  # Subset species
  tifs_species <- tifs[stringr::str_detect(tifs, sp_name)]
  gbif <- tifs_species[stringr::str_detect(tifs_species, "GBIF")]
  ov <- tifs_species[stringr::str_detect(tifs_species, "Overall")]

  # Import ov raster
  ov_ras <- rast(ov)
  # gbif_ras <- rast(gbif)

  # Check if GBIF raster exists, otherwise create empty raster
  if (length(gbif) == 0) {
    gbif_ras <- ov_ras
    values(gbif_ras) <- 0  # Or NA, depending on your logic
  } else {
    gbif_ras <- rast(gbif)
  }

  # Calculating the difference
  diff <- (gbif_ras - 2*(ov_ras))
  names(diff) <- "diff_dist"
  
  # table(getValues(diff))
  
  # stacking rasters
  dist_vars <- c(diff, vars)
  names(dist_vars)
  
  # Converting to a dataframe
  diff_df <- as.data.frame(dist_vars, xy = TRUE, na.rm = TRUE)
  colnames(diff_df)[1:2] <- c("lon", "lat")
  
  # Adding species details
  diff_df <- diff_df %>% 
    mutate(species = sp_name) %>% 
    filter(diff_dist != 0)

  # Combining with the empty dataframe
  df <- rbind(df, diff_df)

}, silent = FALSE)

# Export output
write_csv(df, "output/diff_dist.csv")

############################################
# Import data file
diff_dist <- read_csv("output/diff_dist.csv")
occ <- read_csv("data/com_data.csv")

# Extracting taxonomic information
taxa <- occ %>% 
  dplyr::select(taxa, species) %>% 
  dplyr::distinct(taxa, species, .keep_all = TRUE)

# Renaming species names to match the dataset
taxa$species <- gsub(" ", "_", taxa$species)

# Merging with the CH data
diff_dist_taxa <- dplyr::left_join(diff_dist, taxa, by = "species")

# Convert to categorical for clarity
diff_dist_taxa$diff_dist <- factor(diff_dist$diff_dist, 
                              levels = c(-2, -1), 
                              labels = c("Facebook", "GBIF"))

colnames(diff_dist_taxa)[4:8] <- c("Built areas", "Human footprint", "Temperature",
                              "Rainfall", "Elevation")

# Reshape data for facet plotting
diff_long <- diff_dist_taxa %>%
  pivot_longer(cols = c("Built areas", "Human footprint", "Temperature",
                        "Rainfall", "Elevation"),
               names_to = "predictor", values_to = "value")

# Violin + boxplot overlay
# ggplot(diff_long, aes(x = diff_dist, y = value, fill = diff_dist)) +
#   geom_violin(trim = FALSE, alpha = 0.4) +
#   geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
#   facet_wrap(~predictor, scales = "free_y") +
#   theme_classic() +
#   labs(x = "Data source", y = "Value", fill = "Source") +
#   scale_fill_manual(values = c("Facebook" = "tomato", "GBIF" = "deepskyblue3"))
# 
# ggplot(diff_long, aes(x = diff_dist, y = value, fill = diff_dist)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.8) +
#   facet_grid(predictor ~ taxa, scales = "free_y") +
#   theme_classic(base_size = 12) +
#   labs(x = "Data source", y = "Value", fill = "Source") +
#   scale_fill_manual(values = c("Facebook" = "tomato", "GBIF" = "deepskyblue3")) +
#   theme(strip.text.x = element_text(size = 10, face = "bold"),
#         strip.text.y = element_text(size = 10))

diff_summary <- diff_long %>%
  group_by(taxa, diff_dist, predictor) %>%
  summarise(mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = diff_dist, values_from = mean) %>%
  mutate(diff = GBIF - Facebook)

ggplot(diff_summary, aes(x = predictor, y = diff)) +
  geom_bar(stat = "identity", position = position_dodge(), fill = "skyblue") +
  coord_flip() +
  facet_wrap(~ taxa, scales = "free_x") +
  theme_classic() +
  labs(y = "Mean difference", x = NULL)

# Exporting output
ggsave("output/figures/taxa_source_vars_dist.png")

###############################
# Summary
summary_table <- diff_long %>%
  group_by(taxa, diff_dist, predictor) %>%
  summarise(median = median(value, na.rm = TRUE),
            mean = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE),
            n = n(), .groups = "drop")

summary_stats <- diff_long %>%
  group_by(taxa, diff_dist, predictor) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            median = median(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE),
            n = n(), .groups = "drop")

summary_long <- summary_stats %>%
  pivot_longer(cols = c(mean, median, sd, n),
               names_to = "stat", values_to = "stat_value")

# Export output
write_csv(summary_long, "output/summary_long.csv")

####################################
# Difference with original data
# Import occ data
occ <- read_csv("data/com_data.csv")

# Importing predictor variables
built <- rast("data/built_areas_bd.tif")
hfp <- rast("data/hfp_bd.tif")
clim <- rast("data/climate_bd.tif") # temperature, rainfall, and elevation

# stacking rasters
vars <- c(built, hfp, clim)
names(vars)

# Convert to spatial points
occ_spat <- vect(occ, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")  # WGS84

# Reproject to match raster CRS
occ_proj <- project(occ_spat, crs(vars))  # Match raster CRS

# Extract predictor values
vals <- extract(vars, occ_proj)

# Join with original data
occ_vars <- bind_cols(occ, vals %>% select(-ID))

colnames(occ_vars)[12:16] <- c("Built areas", "Human footprint", "Temperature",
                                   "Rainfall", "Elevation")

# Pivot to long format for plotting
occ_long <- occ_vars %>%
  pivot_longer(cols = c("Built areas", "Human footprint", "Temperature",
                        "Rainfall", "Elevation"),
               names_to = "predictor", values_to = "value")

# Exporting output
write_csv(occ_long, "data/occ_long_var.csv")

# Summary difference by source
diff_summary <- occ_long %>%
  group_by(taxa, source, predictor) %>%
  summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = source, values_from = mean_val) %>%
  mutate(diff = Facebook - GBIF)

# Bar plot of differences
ggplot(diff_summary, aes(x = predictor, y = diff)) +
  geom_col(fill = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(y = "Facebook – GBIF (mean)", x = "Predictor") + facet_wrap(~taxa)


# Reshape to long format
occ_long <- occ_vars %>%
  pivot_longer(cols = c("Built areas", "Human footprint", "Temperature",
                        "Rainfall", "Elevation"),
               names_to = "predictor", values_to = "value")

# Calculate range and median per taxa and source
summary_stats <- occ_long %>%
  group_by(taxa, source, predictor) %>%
  summarise(min_val = min(value, na.rm = TRUE),
            max_val = max(value, na.rm = TRUE),
            med_val = median(value, na.rm = TRUE),
            .groups = "drop")


summary_stats_log <- summary_stats %>%
  mutate(min_val_log = log10(min_val + 1),
         max_val_log = log10(max_val + 1),
         med_val_log = log10(med_val + 1))


ggplot(summary_stats_log, aes(x = predictor, group = source)) +
  geom_linerange(aes(ymin = min_val_log, ymax = max_val_log, colour = source),
                 position = position_dodge(width = 0.5), linewidth = 1) +
  geom_point(aes(y = med_val_log, fill = source),
             shape = 21, colour = "black",
             position = position_dodge(width = 0.5), size = 3) +
  facet_wrap(~ taxa, scales = "free_y") +
  scale_colour_manual(values = c("Facebook" = "tomato", "GBIF" = "deepskyblue3")) +
  scale_fill_manual(values = c("Facebook" = "tomato", "GBIF" = "deepskyblue3")) +
  theme_classic(base_size = 13) +
  labs(x = "", y = "", colour = "Source", fill = "Source") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("output/figures/taxa_source_vars_data.png")

###################################
# Linear model
# Import data
occ_long <- read_csv("data/occ_long_var.csv")

# Standardise values within each predictor
occ_long_scaled <- occ_long %>%
  group_by(predictor) %>%
  mutate(value_scaled = scale(value)) %>%
  ungroup()

occ_long_scaled$source <- factor(occ_long_scaled$source, levels = c("GBIF", "Facebook"))

# Fit linear models for each predictor
effect_by_taxa <- occ_long_scaled %>%
  group_by(predictor, taxa) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(value_scaled ~ source, data = .x)),
    result = map(model, tidy)
  ) %>%
  unnest(result) %>%
  filter(term == "sourceFacebook")  # effect = difference Facebook vs GBIF

effect_by_taxa <- effect_by_taxa %>%
  mutate(signif = ifelse(p.value < 0.05, "Significant", "Not significant"))

# Export stat table
write_csv(effect_by_taxa, "output/effect_taxa.csv")

# Plot the effect sizes
ggplot(effect_by_taxa, aes(x = estimate, y = fct_reorder(predictor, estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(xmin = estimate - std.error, 
                      xmax = estimate + std.error,
                      colour = signif)) +
  scale_colour_manual(values = c("Significant" = "darkgoldenrod1", 
                                 "Not significant" = "skyblue")) +
  facet_wrap(~ taxa, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(
    x = "Effect size (Facebook vs GBIF)",
    y = "",
    colour = "Effect"
  )

# Export figure
ggsave("output/figures/taxa_source_vars_data_stat.png")
