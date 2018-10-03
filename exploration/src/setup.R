
print("setup.R is running...")

##################
# Load libraries #
##################

library(phyloseq)
library(ggplot2)
library(sets)
library(xlsx)

###############
# Import data #
###############

input_data_files <- dir(path = paste(location, "/data/input/", sep =""), pattern = "(?i)[.]", 
                        recursive = TRUE, full.names = TRUE, all.files = TRUE)
lapply(input_data_files,  load, .GlobalEnv)

#################
# End of script #
#################

print("setup.R is done.")