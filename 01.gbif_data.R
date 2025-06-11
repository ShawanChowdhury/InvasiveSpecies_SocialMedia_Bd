# Loading required packages
library(tidyverse)
library(rgbif)

# Species name
species <- read_csv("data/species_list.csv")
species <- species$species

# match the names 
gbif_taxon_keys <- species %>% 
  name_backbone_checklist() %>% # match to backbone 
  filter(!matchType == "NONE") %>% # get matched names
  pull(usageKey) 

# download the data
occ_download(
  pred_in("taxonKey", gbif_taxon_keys), # important to use pred_in
  pred("hasCoordinate", TRUE),
  pred("occurrenceStatus","PRESENT"),
  pred("country", "BD"),
  format = "SIMPLE_CSV",
  user= "shawan_zl",pwd= "nabolakothaFB89",email= "shawan1094061@gmail.com"
)

# Check download status
occ_download_wait('0002574-250525065834625')

# Citation
# GBIF Occurrence Download https://www.gbif.org/occurrence/download/0002603-250525065834625 Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2025-05-27
