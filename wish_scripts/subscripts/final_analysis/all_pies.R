###################
# Predicted hosts #
###################

# Hosts predicted for contigs appearing in the samples of the first week of life

#viruses and predictions
View(sample_data(vir_ps))
vir_ps_pruned <- prune_taxa(taxa_sums(vir_ps)>=0.0001, vir_ps)
vir_contigs <- data.frame(contigs = row.names(tax_table(vir_ps_pruned)))

vir_ps_predictions <- merge(x = vir_contigs, y = pred_lca_5_p_full_k7,
                         by.x = "contigs", by.y = "contig")
View(vir_ps_predictions)
plot_pie_of_preds(df = vir_ps_predictions, column = 3, trsh = 3)
vir_ps_taxonomy <- t(sapply(as.character(vir_ps_predictions[,3]), function(x) get_fop(x)))
View(vir_ps_taxonomy)


# predictions at a phylum level:
plot_pie_of_preds(df = as.data.frame(na.omit(vir_ps_taxonomy[,4])), 
                  column = 1, trsh = 3, lab = 'Phylum level predictions \n in all samples')
# predictions at an order level:
plot_pie_of_preds(df = as.data.frame(na.omit(vir_ps_taxonomy[,3])), 
                  column = 1, trsh = 3, lab = 'Order level predictions \n in all samples')
# predictions at a family level:
plot_pie_of_preds(df = as.data.frame(na.omit(vir_ps_taxonomy[,2])), 
                  column = 1, trsh = 3, lab = 'Family level predictions \n in all samples')
# predictions at a geus level (only the ones that were predicted up to a genus level)
pred_gen_all <- vir_ps_predictions[vir_ps_predictions$lca_rank == "genus",]
plot_pie_of_preds(df = pred_gen_all, column = 3, lab = 'Genus level predictions \n in all samples' )

############ 
# Bacteria #
############

# The actual bacterial presence/abundance in the samples of the first week of life
# Only prokaryotes that are more abundant than 0.0001 are considered

bac_ps_pruned <- prune_taxa(taxa_sums(bac_ps)>=0.0001, bac_ps)
bac_contigs <- as.data.frame(tax_table(bac_ps_pruned))
#dim(bac_contigs)
#View(bac_contigs)
#View(sample_data(bac_ps_first_week))


# Phyum level prokaryotic presence
phy <- plot_pie_of_preds(df = as.data.frame(na.omit(bac_contigs[,2])), 
                  column = 1, trsh = 3, lab = 'Phylum level presence\n in first week samples')
# Order level prokaryotic 
ord <- plot_pie_of_preds(df = as.data.frame(na.omit(bac_contigs[,4])), 
                  column = 1, trsh = 3, lab = 'Order level presence\n in first week samples')
# Family level prokaryotic
fam <- plot_pie_of_preds(df = as.data.frame(na.omit(bac_contigs[,5])), 
                  column = 1, trsh = 3, lab = 'Family level presence\n in first week samples')
# Genus level presence 
gen <- plot_pie_of_preds(df = as.data.frame(na.omit(bac_contigs[,6])), column = 1, trsh = 3,
                  lab = 'Genus level presence\n in first week samples')

########
# Save #
########

save(list = c("vir_ps_predictions", "vir_ps_taxonomy", "bac_contigs"), 
     file = paste(location, "data/intermediate/final_analysis/full_pies.RData", sep = ""))
