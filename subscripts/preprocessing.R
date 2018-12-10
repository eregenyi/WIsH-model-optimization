print("preprocessing.R is running...")

###############
# Format data #
###############

# prok_list
prok_list <- data.frame(prok_list)
names(prok_list) <- c("file_name", "accession_number", "genus", "species", "strain")
prok_list$genus[prok_list$genus == "Candidatus"] <- prok_list$species[prok_list$genus == "Candidatus"]
prok_list$species[prok_list$genus == prok_list$species] <- prok_list$strain[prok_list$genus == prok_list$species]

# phages_list
names(phages_list) <- c("file_name", "contig_name")

# k6 predictions and likelihood data
names(ll_k6) <- file_to_contig_names(names(ll_k6), phages_list)
row.names(ll_k6) <- gsub("_genomic", "", row.names(ll_k6), fixed = TRUE)
predictions_k6$contig <- file_to_contig_names(predictions_k6$contig, phages_list)

# k7
names(ll_k7) <- file_to_contig_names(names(ll_k7), phages_list)
row.names(ll_k7) <- gsub("_genomic", "", row.names(ll_k7), fixed = TRUE)
names(predictions_k7) <- c("contig", "best_hit", "likelihood", "p-value")
predictions_k7$contig <- file_to_contig_names(predictions_k7$contig, phages_list)

# k8
names(ll_k8) <- file_to_contig_names(names(ll_k8), phages_list)
row.names(ll_k8) <- gsub("_genomic", "", row.names(ll_k8), fixed = TRUE)
names(predictions_k8) <- c("contig", "best_hit", "likelihood", "p-value")
predictions_k8$contig <- file_to_contig_names(predictions_k8$contig, phages_list)

##########################
# Get top 10 predictions #
##########################

contigs <- blastn_host$contig

# k6
pred_top_10_k6 <- parse_into_top_10_predictions(contigs, ll_k6, predictions_k6, prok_list)

# k7
pred_top_10_k7 <- parse_into_top_10_predictions(contigs, ll_k7, predictions_k7, prok_list)

# k8
pred_top_10_k8 <- parse_into_top_10_predictions(contigs, ll_k8, predictions_k8, prok_list)

#############################
# Save intermediate results #
#############################

save(list = c("pred_top_10_k6", "pred_top_10_k7", "pred_top_10_k8"), file = paste(location, "data/intermediate/preprocessing/pred_top_10.RData", sep = ""))

#################
# End of script #
#################

print("preprocessing.R is done.")
