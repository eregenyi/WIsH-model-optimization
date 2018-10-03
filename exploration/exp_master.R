"
###############
# Information #
###############

These scripts aim to explore longitudinal viral and bacterial metagenomic data (in a phyloseq format)
- contig lenth
- diversity and richness of samples
- contig similarity (co-occurence)

"

# Clear work space
rm(list = ls())

# Clear the terminal
cat("\014")

# Getlocation
location <- (rstudioapi::getActiveDocumentContext()$path)
location <- gsub("exp_master.R", "", location)

# Load data from the intermediate folder (In case you have already had some work done)
intermediate_data_files <- dir(path = paste(location, "/data/intermediate/", sep =""), pattern = "(?i)[.]", recursive = TRUE, full.names = TRUE, all.files = TRUE)
lapply(intermediate_data_files,  load, .GlobalEnv)

# Source subscripts
#source(paste(location, "src/installation.R", sep = "")) # This script installs every library on your computer that you need, you only need to run it once
source(paste(location, "src/setup.R", sep = "")) # Load needed libraries & data (input data files have to be in the data/input folder)
source(paste(location, "src/functions.R", sep = ""))
source(paste(location, "src/preprocessing.R", sep = "")) # Explore in general: sequence length etc.
source(paste(location, "src/div_rich.R", sep = "")) 
source(paste(location, "src/co_occurence_viral.R", sep = "")) #note: occurrence is 2 r-s
source(paste(location, "src/co_occurence_phage_host.R", sep = "")) #note: occurrence is 2 r-s

