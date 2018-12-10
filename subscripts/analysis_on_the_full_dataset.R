# Based on the results of the previous analyses on my data, I chose a model of order 7, 
# and decided to use a p-value cutoff of 0.06, along with using the LCA of the 5 best predictions.

print("analysis_on_the_full_dataset.R is running...")

#########################
# Get final predictions #
#########################

# Create data frame with top 10 predictions
pred_top_10_full_k7 <- parse_into_top_10_predictions(predictions_k7$contig, ll_k7, predictions_k7, prok_list)

# Apply a p-value cutoff
pred_top_10_full_k7_p  <- pred_top_10_full_k7[pred_top_10_full_k7$p.value < 0.05,]

# Get lcas
pred_lca_5 <- merge_into_lca(pred_top_10_full_k7_p, col_taxa = 3:7, id = TRUE)
pred_lca_5_full_k7 <- data.frame(pred_top_10_full_k7_p[,1:2], pred_lca_5)

# Did it work:
sum(rowSums(is.na(pred_lca_5_p_full_k7))) # 3 na -> for some reason here taxize got stuck. Assign it manually
which( 0 < rowSums(is.na(pred_lca_5_p_full_k7))) # obesrvation 478 (rowname: 1445)
# check if I reassign the right observation:
to_reassign <- which( 0 < rowSums(is.na(pred_lca_5_p_full_k7)))
View(pred_top_10_full_k7_p[to_reassign,])
View(pred_lca_5_p_full_k7[to_reassign,])
# Compute manually:
pred_lca_5_p_full_k7[to_reassign,3:5] <- merge_into_lca(pred_top_10_full_k7_p[to_reassign,], col_taxa = 3:7, id = TRUE)

# Did it work now?
sum(rowSums(is.na(pred_lca_5_p_full_k7))) # YUSS!

############
# Analysis #
############

# number of sequences
nrow(pred_lca_5_p_full_k7) 

# plot of taxonomic rank distribution among the final predictions
plot_lca_5_full_k7 <- plot_bar_tax_ranks(pred_lca_5_p_full_k7,  col_rank = 4, colour = col_k7)
plot_lca_5_full_k7 + ylim(c(0,40))

#######################################
# Link prokaryotes to the final model #
#######################################
# determine host-composition for genus, family, and oder level predictions in the final dataset

# 1. genus level predictions
pred_genus_lev <- pred_lca_5_p_full_k7[pred_lca_5_p_full_k7$lca_rank == "genus",]
plot_pie_of_preds(pred_genus_lev, lab = "Pie Chart of Genus Level Predictions")
# Get family, order and phylum (the last one is not needed, but I already have a function for that)
unique_perc_gen <- sum_perc_of_unique(pred_genus_lev[,3])
unique_gen_fop <- merge(x = unique_perc_gen, y = get_tax(unique_perc_gen[,1]), by.x = "unique_elements", by.y = "unique_elements")

# 2. family level predictions
pred_fam_lev <- pred_lca_5_p_full_k7[pred_lca_5_p_full_k7$lca_rank == "family",]
# just family level predictions
plot_pie_of_preds(pred_fam_lev, lab = "Pie Chart of Family Level Predictions")
# family level predictions + genus level considered at a family level
unique_perc_fam <- sum_perc_of_unique(pred_fam_lev[,3])
unique_count_fam <- collapse_and_sum(unique_gen_fop, 5, 3)
merged_fam_gen <- merge_and_sum(dfa= unique_perc_fam, dfb=unique_count_fam, c1a = 1, c2a = 3, c1b = 1, c2b = 2)
plot_pie_perc(merged_fam_gen, 1, 3, lab = "Pie Chart of Family and Genus, cosidered at an order level")

unique_fam_fop <- merge(x = merged_fam_gen, y = get_tax(merged_fam_gen[,1]), by.x = "family", by.y = "family")

# 3. order level predictions
pred_ord_lev <- pred_lca_5_p_full_k7[pred_lca_5_p_full_k7$lca_rank == "order",]
# just order level predictions
plot_pie_of_preds(pred_ord_lev, lab = "Pie Chart of Order Level Predictions")
# order level + family and genus level considered at an order level
unique_perc_ord <- sum_perc_of_unique(pred_ord_lev[,3])
unique_count_ord <- collapse_and_sum(unique_fam_fop, 6, 2)
merged_fam_gen_ord <- merge_and_sum(dfa= unique_perc_ord, dfb=unique_count_ord, c1a = 1, c2a = 3, c1b = 1, c2b = 2)
plot_pie_perc(merged_fam_gen_ord, 1, 3, lab = "Pie Chart of Order, Family and Genus, cosidered at an order level")

unique_ord_fop <- merge(x = merged_fam_gen_ord, y = get_tax(merged_fam_gen_ord[,1]), by.x = "family", by.y = "order")

# 4. phylum level predictions
pred_phy_lev <- pred_lca_5_p_full_k7[pred_lca_5_p_full_k7$lca_rank == "phylum",]
# just order level predictions
plot_pie_of_preds(pred_phy_lev, lab = "Pie Chart of Phylum Level Predictions")
# order level + family and genus level considered at an order level
unique_perc_phy <- sum_perc_of_unique(pred_phy_lev[,3])
unique_count_phy <- collapse_and_sum(unique_fam_fop, 7, 2)
merged_fam_gen_ord_phy <- merge_and_sum(dfa= unique_perc_phy, dfb=unique_count_phy, c1a = 1, c2a = 3, c1b = 1, c2b = 2)
plot_pie_perc(merged_fam_gen_ord_phy, 1, 3, lab = "Pie Chart of Phylum, Order, Family and Genus, cosidered at a Phylum level")


############
# Heatmaps #
############

## Bacteria

#first get the phages likelihood profiles from the ll_k7 matrix
retain <- which(names(ll_k7) %in% pred_lca_5_p_full_k7$contig)
ll_p_k7 <- ll_k7[,retain]
dim(ll_p_k7)

#plot heatmap: bacterial genera cluster together! and this is only a subset!
plot_hmp(ll_p_k7, seed = 1) #seed is a "base" for the random number generator. This was, your results will be reproducible, but random
plot_hmp(ll_p_k7, seed = 2)
plot_hmp(ll_p_k7, seed = 3)
plot_hmp(ll_p_k7, seed = 4)

## Phages
#take a subset for which the phages are known (genus/family/order) -> do they cluster together? 

# by genera

ll_t_p_k7 <- t(ll_p_k7)
retain_phage <- which(row.names(ll_t_p_k7) %in% test_phage_clustering$contig)
ll_t_p_k7 <- ll_t_p_k7[retain_phage,]
plot_hmp(ll_t_p_k7, seed = 1, prok_l = test_phage_clustering, col_fname = 1, col_gen = 4, sample_row = FALSE, sample_col = FALSE, sample_size_row = 1, sample_size_col = 1)

# by family

ll_t_p_k7 <- t(ll_p_k7)
retain_phage <- which(row.names(ll_t_p_k7) %in% test_phage_clustering$contig)
ll_t_p_k7 <- ll_t_p_k7[retain_phage,]
plot_hmp(ll_t_p_k7, seed = 1, prok_l = test_phage_clustering, col_fname = 1, col_gen = 3, sample_row = FALSE, sample_col = FALSE, sample_size_row = 1, sample_size_col = 1)

ll_p_k7 <- ll_k7[,retain]
ll_t_p_k7_f <- t(ll_p_k7)
retain_phage <- which(row.names(ll_t_p_k7_f) %in% kt_family$V1)
ll_t_p_k7_f <- ll_t_p_k7_f[retain_phage,]
plot_hmp(ll_t_p_k7_f, seed = 1, prok_l = kt_family, col_fname = 1, col_gen = 5, sample_row = TRUE, sample_col = FALSE, sample_size_row = 100, sample_size_col = 1)

# by host genera and host family
ll_p_k7 <- ll_k7[,retain]
ll_t_p_k7_h <- t(ll_p_k7)
retain_phage <- which(row.names(ll_t_p_k7_h) %in% blastn_host$contig)
ll_t_p_k7_h <- ll_t_p_k7_h[retain_phage,]
family <- tax_name(levels(blastn_host$blastn_host), db = "ncbi", get = "family" , rows = 2)[2:3]
blastn_host <- merge(x = blastn_host, y = family, by.x = "blastn_host", by.y = "query")
blastn_host <- blastn_host[,c(2,1,3)]
plot_hmp(ll_t_p_k7_h, seed = 1, prok_l = blastn_host, col_fname = 1, col_gen = 2, sample_row = FALSE, sample_col = FALSE, sample_size_row = 1, sample_size_col = 1)
plot_hmp(ll_t_p_k7_h, seed = 1, prok_l = blastn_host, col_fname = 1, col_gen = 3, sample_row = FALSE, sample_col = FALSE, sample_size_row = 1, sample_size_col = 1)


####################
# Cluster analysis #
####################
# Since most of the data is not annotated, unsupervised methods are needed.
# as in an OTU table, and we have a much greater "sampe size", 
# so it could be fine to use a parametric approach

source(paste(location, "subscripts/final_analysis/pca.R", sep = ""))

# clustering? (NMDA) --> this should be done on a dataset/subset where the phages are known to a family/genus/order level, to see if its useful


########
# Save #
########

save(list = c("pred_top_10_full_k7", "pred_lca_5_full_k7", "pred_lca_5_p_full_k7", "plot_lca_5_full_k7", 
              "plot_hmp_k7_1", "plot_hmp_k7_2", "plot_hmp_k7_3", "plot_hmp_k7_4", "ll_p_k7", "test_phage_clustering",
              "unique_gen_fop", "unique_gen_fop",
              "unique_perc_gen", "merged_fam_gen", "merged_fam_gen_ord"), 
     file = paste(location, "data/intermediate/final_analysis/final_analysis.RData", sep = ""))

# write.table(pred_lca_5_p_full_k7, file = paste(location, "data/output/lca5_order7.txt", sep = ""))

#################
# End of script #
#################

print("analysis_on_the_full_dataset.R is done.")

