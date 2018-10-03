# script to check which hosts are predicted in nodes appearing in the first week

###################
# Predicted hosts #
###################

# Hosts predicted for contigs appearing in the samples of the first week of life
# Only contigs that are more abundant than 0.0001 are considered

#viruses and predictions
View(sample_data(vir_ps))
vir_ps_first_week <- subset_samples(vir_ps, X.weeks == 1)
vir_ps_first_week <- prune_taxa(taxa_sums(vir_ps_first_week)>=0.0001, vir_ps_first_week)
first_week_contigs <- data.frame(contigs_first_week = row.names(tax_table(vir_ps_first_week)))

fwc_predictions <- merge(x = first_week_contigs, y = pred_lca_5_p_full_k7,
                         by.x = "contigs_first_week", by.y = "contig")
View(fwc_taxonomy)
plot_pie_of_preds(df = fwc_predictions, column = 3, trsh = 3)
fwc_taxonomy <- t(sapply(as.character(fwc_predictions[,3]), function(x) get_fop(x)))
View(fwc_taxonomy)


# first week predictions at a phylum level:
plot_pie_of_preds(df = as.data.frame(na.omit(fwc_taxonomy[,4])), 
                  column = 1, trsh = 3, lab = 'Phylum level predictions \n in first week samples')
# first week predictions at an order level:
plot_pie_of_preds(df = as.data.frame(na.omit(fwc_taxonomy[,3])), 
                  column = 1, trsh = 3, lab = 'Order level predictions \n in first week samples')
# first week predictions at a family level:
plot_pie_of_preds(df = as.data.frame(na.omit(fwc_taxonomy[,2])), 
                  column = 1, trsh = 3, lab = 'Family level predictions \n in first week samples')
# first week predictions at a geus level (only the ones that were predicted up to a genus level)
pred_gen <- fwc_predictions[fwc_predictions$lca_rank == "genus",]
plot_pie_of_preds(df = pred_gen, column = 3, lab = 'Genus level predictions \n in first week samples' )

############ 
# Bacteria #
############

# The actual bacterial presence/abundance in the samples of the first week of life
# Only prokaryotes that are more abundant than 0.0001 are considered

first_week_samples <- row.names(sample_data(vir_ps_first_week))
first_week_samples <- gsub(".", "-", first_week_samples, fixed = TRUE)
bac_ps_first_week <- subset_samples(bac_ps, unique_ID_old %in% first_week_samples)
bac_ps_first_week <- prune_taxa(taxa_sums(bac_ps_first_week)>=0.0001, bac_ps_first_week)
first_week_bac_contigs <- as.data.frame(tax_table(bac_ps_first_week))
dim(first_week_bac_contigs)
#View(first_week_bac_contigs)
#View(sample_data(bac_ps_first_week))


# Phyum level prokaryotic presence in first week samples
plot_pie_of_preds(df = as.data.frame(na.omit(first_week_bac_contigs[,2])), 
                  column = 1, trsh = 3, lab = 'Phylum level presence\n in first week samples')
# Order level prokaryotic presence in first week samples
plot_pie_of_preds(df = as.data.frame(na.omit(first_week_bac_contigs[,4])), 
                  column = 1, trsh = 3, lab = 'Order level presence\n in first week samples')
# Family level prokaryotic presence in first week samples
plot_pie_of_preds(df = as.data.frame(na.omit(first_week_bac_contigs[,5])), 
                  column = 1, trsh = 3, lab = 'Family level presence\n in first week samples')
# Genus level presence in first week samples
plot_pie_of_preds(df = as.data.frame(na.omit(first_week_bac_contigs[,6])), column = 1, trsh = 3,
                  lab = 'Genus level presence\n in first week samples')

########
# Save #
########

save(list = c("fwc_predictions", "fwc_taxonomy", "first_week_bac_contigs"), 
     file = paste(location, "data/intermediate/final_analysis/first_week.RData", sep = ""))
