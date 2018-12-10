# This script evaluates the performance of different order models using
# 1) the best prediction
# 2) the majority vote of the 10 best predictions

print("best_prediction_and_majority_vote.R is now running...")

###################
# Best prediction #
###################

# k6
best_pred_k6 <- merge(x = pred_top_10_k6[,1:3], y = blastn_host, by.x = "contig", by.y = "contig")
best_pred_fop_k6 <- get_host_fop_pred_and_true(best_pred_k6, 3, 4)

# k7
best_pred_k7 <- merge(x = pred_top_10_k7[,1:3], y = blastn_host, by.x = "contig", by.y = "contig")
best_pred_fop_k7 <- get_host_fop_pred_and_true(best_pred_k7, 3, 4)

# k8
best_pred_k8 <- merge(x = pred_top_10_k8[,1:3], y = blastn_host, by.x = "contig", by.y = "contig")
best_pred_fop_k8 <- get_host_fop_pred_and_true(best_pred_k8, 3, 4)

# accuracies 
bp_acc_no_p <- get_bp_acc_table(best_pred_fop_k6, best_pred_fop_k7, best_pred_fop_k8)
bp_acc_p_006 <- get_bp_acc_table(best_pred_fop_k6, best_pred_fop_k7, best_pred_fop_k8, p = 0.06)
bp_acc_p_005 <- get_bp_acc_table(best_pred_fop_k6, best_pred_fop_k7, best_pred_fop_k8, p = 0.05)
bp_acc_p_001 <- get_bp_acc_table(best_pred_fop_k6, best_pred_fop_k7, best_pred_fop_k8, p = 0.01)
bp_acc <- data.frame(bp_acc_no_p, bp_acc_p_006, bp_acc_p_005, bp_acc_p_001)
names(bp_acc) <- c("order6_no_p", "order7_no_p", "order8_no_p", "order6_p0.06", "order7_p0.06", "order8_p0.06",
                  "order6_p0.05", "order7_p0.05", "order8_p0.05", "order6_p0.01", "order7_p0.01", "order8_p0.01")
rownames(bp_acc) <- c("genus", "family", "order", "phylum")
View(bp_acc)



#################
# Majority vote #
#################

# k6
maj_pred_k6_temp <- data.frame(pred_top_10_k6[,1:2], pred_host_gen = apply(pred_top_10_k6[,3:12], 1, function(x) get_mode(x)))
maj_pred_k6 <- merge(x = maj_pred_k6_temp, y = blastn_host, by.x = "contig", by.y = "contig")
maj_pred_fop_k6 <- get_host_fop_pred_and_true(maj_pred_k6, 3, 4)

# k7
maj_pred_k7_temp <- data.frame(pred_top_10_k7[,1:2], pred_host_gen = apply(pred_top_10_k7[,3:12], 1, function(x) get_mode(x)))
maj_pred_k7 <- merge(x = maj_pred_k7_temp, y = blastn_host, by.x = "contig", by.y = "contig")
maj_pred_fop_k7 <- get_host_fop_pred_and_true(maj_pred_k7, 3, 4)

# k8
maj_pred_k8_temp <- data.frame(pred_top_10_k8[,1:2], pred_host_gen = apply(pred_top_10_k8[,3:12], 1, function(x) get_mode(x)))
maj_pred_k8 <- merge(x = maj_pred_k8_temp, y = blastn_host, by.x = "contig", by.y = "contig")
maj_pred_fop_k8 <- get_host_fop_pred_and_true(maj_pred_k8, 3, 4)

################
# Create plots #
################

col_k6 <- viridis(3)[1]
col_k7 <- viridis(3)[2]
col_k8 <- plasma(5)[4]

#---> plot the best predictions 

  # genus, order, family and phylum level accuracies vs contig length
g_bp_l_6 <- plot_al_gfop(best_pred_fop_k6, color = col_k6) 
g_bp_l_7 <- plot_al_gfop(best_pred_fop_k7, color = col_k7) 
g_bp_l_8 <- plot_al_gfop(best_pred_fop_k8, color = col_k8)
  # genus, order, family and phylum level accuracies vs p-value cutoff
g_bp_p_6 <- plot_ap_gfop(best_pred_fop_k6, color = col_k6) 
g_bp_p_7 <- plot_ap_gfop(best_pred_fop_k7, color = col_k7) 
g_bp_p_8 <- plot_ap_gfop(best_pred_fop_k8, color = col_k8) 
  # number of remaining observations in the test set vs p-value cutoff
g_bp_n_6 <- plot_nobs_vs_cutoff(best_pred_fop_k6, color = col_k6)
g_bp_n_7 <- plot_nobs_vs_cutoff(best_pred_fop_k7, color = col_k7)
g_bp_n_8 <- plot_nobs_vs_cutoff(best_pred_fop_k8, color = col_k8)

#---> plot the majority votes

  # genus, order, family and phylum level accuracies vs contig length
g_mv_l_6 <- plot_al_gfop(maj_pred_fop_k6, color = col_k6) 
g_mv_l_7 <- plot_al_gfop(maj_pred_fop_k7, color = col_k7) 
g_mv_l_8 <- plot_al_gfop(maj_pred_fop_k8, color = col_k8) 
  # genus, order, family and phylum level accuracies vs p-value cutoff
g_mv_p_6 <- plot_ap_gfop(maj_pred_fop_k6, color = col_k6) 
g_mv_p_7 <- plot_ap_gfop(maj_pred_fop_k7, color = col_k7) 
g_mv_p_8 <- plot_ap_gfop(maj_pred_fop_k8, color = col_k8) 
  # number of remaining observations in the test set vs p-value cutoff
g_mv_n_6 <- plot_nobs_vs_cutoff(maj_pred_fop_k6, color = col_k6)
g_mv_n_7 <- plot_nobs_vs_cutoff(maj_pred_fop_k7, color = col_k7)
g_mv_n_8 <- plot_nobs_vs_cutoff(maj_pred_fop_k8, color = col_k8)

###################
# Check plots out #
###################

grid.arrange(g_bp_l_6, g_bp_l_7, g_bp_l_8,
             g_bp_p_6, g_bp_p_7, g_bp_p_8,
             g_bp_n_6, g_bp_n_7, g_bp_n_8,
             nrow = 3, ncol = 3)

grid.arrange(g_mv_l_6, g_mv_l_7, g_mv_l_8,
             g_mv_p_6, g_mv_p_7, g_mv_p_8,
             g_mv_n_6, g_mv_n_7, g_mv_n_8,
             nrow = 3, ncol = 3)

########
# Save #
########

print("Saving results of best_prediction_and_majority_vote.R...")

# data files
save(list = c("best_pred_fop_k6", "best_pred_fop_k7", "best_pred_fop_k8"), 
     file = paste(location, "data/intermediate/best_and_majority_prediction_results/best_pred.RData", sep = ""))
save(list = c("maj_pred_fop_k6", "maj_pred_fop_k7", "maj_pred_fop_k8"), 
     file = paste(location, "data/intermediate/best_and_majority_prediction_results/maj_pred.RData", sep = ""))

# plots
save(list = c("g_bp_l_6", "g_bp_l_7", "g_bp_l_8", "g_bp_p_6", "g_bp_p_7", "g_bp_p_8", "g_bp_n_6", "g_bp_n_7", "g_bp_n_8"), 
     file = paste(location, "data/intermediate/best_and_majority_prediction_results/plots_best_pred.RData", sep = ""))
save(list = c("g_mv_l_6", "g_mv_l_7", "g_mv_l_8", "g_mv_p_6", "g_mv_p_7", "g_mv_p_8", "g_mv_n_6", "g_mv_n_7", "g_mv_n_8"), 
     file = paste(location, "data/intermediate/best_and_majority_prediction_results/plots_maj_pred.RData", sep = ""))

#################
# End of script #
#################

print("best_prediction_and_majority_vote.R is done.")