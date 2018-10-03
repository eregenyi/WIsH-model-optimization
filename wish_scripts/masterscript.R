"
######################
# Module Information #
######################

The present scripts were written to process output files of WIsH (Who Is the Host).
The goal of the script is to assess the prediction accuracy of models of different order,
using various p-value cutoffs and strategies that could help improve the accuracy. 
The resulting visualizations aim to aid decision making on:
- which order model is the best for the particular dataset 
- what is a sensible p-value cut-off
- what strategy to follow, to achieve the best accuracy
  (e.g. using the best prediction, the majority vote of the 10 best predictions,
  or use the lowest common ancestor of the best X predictions - if ceratain loss of 
  the prediction rank can be tolerated)
- the clustering section explores possible grouping structures in the log-likelihood matrix

DISCLAIMER
- halfway through I realized that using taxize may keep it clean and avoid the need of local databases (and updates),
  but it also makes some parts of the code very slow (it's a package that retrieves taxonomic information directly from databases)
  If you have decent internet and your test dataset is the order of hundreds of observations, it is probably going to run overnight. (Y)
  Otherwise maybe consider modifing my script so taxonomic data is retrieved from a local database
- taxize can have bugs. If it prompts you for input, reinstall it (the choice should be automated by the rows = 2 attribute)
  If any other taxize functions don't work as expected, try reinstalling the package, or ask them at their github page (in issues), 
  to me, they answered immediately. 

"

# Clear work space
rm(list = ls())

# Clear the terminal
cat("\014")

# Get masterscript location
location <- (rstudioapi::getActiveDocumentContext()$path)
location <- gsub("masterscript.R", "", location)

# Load data from the intermediate folder (In case you have already had some work done)
intermediate_data_files <- dir(path = paste(location, "/data/intermediate/", sep =""), pattern = "(?i)[.]", recursive = TRUE, full.names = TRUE, all.files = TRUE)
lapply(intermediate_data_files,  load, .GlobalEnv)

# Source subscripts
#source(paste(location, "subscripts/installation.R", sep = "")) # This script installs every library on your computer that you need. you only need to run it once
source(paste(location, "subscripts/setup.R", sep = "")) # Load needed libraries 
source(paste(location, "subscripts/functions.R", sep = "")) # Run the script with the custom functions
source(paste(location, "subscripts/import_data.R", sep = "")) # Import datasets (initial datasets: input data files have to be in the data/input folder)
source(paste(location, "subscripts/preprocessing.R", sep = "")) # Create a dataframe with the contig name, the p-value, and the 10 best predictons 
source(paste(location, "subscripts/exploration.R", sep = "")) # Explore in general: sequence length etc.
source(paste(location, "subscripts/best_prediction_and_majority_vote.R", sep = "")) # Evaluates the performance of models using the best prediction and the majority vote
source(paste(location, "subscripts/lca.R", sep = "")) # Evaluates the performance of models 
source(paste(location, "subscripts/analysis_on_the_full_dataset.R", sep = "")) # Analysis on the full dataset

