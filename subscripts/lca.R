# This script aimes to test how the overall accuracy 
# improves upon taking the lca of the first 2-10 best predictions
# The main results of this script are plots. Paring the data is outsourced to
# separate scripts for the sake of keeping the code as clean as possible. see: lca folder

print("lca.R is running...")

# variables that don't belong to any specific parts of the script
lca_of_first_x_pred <- seq(1,10,1)

col_k6 <- viridis(3)[1]
col_k7 <- viridis(3)[2]
col_k8 <- plasma(5)[4]

###############
# For order 6 #
###############

# get accuracies for using the lowest common ancestor of the 2...10 best prediction
source(paste(location, "subscripts/lca/lca_k6.R", sep = ""))

# plot accuracy vs lca of best X
pred_acc_vs_lca_k6 <- c(mean(best_pred_fop_k6$is_corr_gen), mean(lca_2_pred_k6$is_correct), mean(lca_3_pred_k6$is_correct), mean(lca_4_pred_k6$is_correct),
                        mean(lca_5_pred_k6$is_correct), mean(lca_6_pred_k6$is_correct), mean(lca_7_pred_k6$is_correct),
                        mean(lca_8_pred_k6$is_correct), mean(lca_9_pred_k6$is_correct), mean(lca_10_pred_k6$is_correct))
best_x_vs_acc_k6 <- data.frame(lca_of_first_x_pred, pred_acc_vs_lca_k6)
lca_acc_vs_nr_bp_k6 <- plot_acc_vs_nr_bp(best_x_vs_acc_k6)

# plot accuracy and taxonomic depth vs x best
best_x_vs_acc_and_genperc_k6 <- data.frame(best_x_vs_acc_k6, 
                                           genperc = c(1, get_gen_rate(lca_2_pred_k6$lca_rank), get_gen_rate(lca_3_pred_k6$lca_rank), get_gen_rate(lca_4_pred_k6$lca_rank),
                                                       get_gen_rate(lca_5_pred_k6$lca_rank), get_gen_rate(lca_6_pred_k6$lca_rank), get_gen_rate(lca_7_pred_k6$lca_rank),
                                                       get_gen_rate(lca_8_pred_k6$lca_rank), get_gen_rate(lca_9_pred_k6$lca_rank), get_gen_rate(lca_10_pred_k6$lca_rank)))
acc_and_genperc_vs_nr_bp_k6 <- plot_acc_and_genperc_vs_nr_bp(best_x_vs_acc_and_genperc_k6, color_acc = col_k6, color_perc = lighten(col_k6,factor = 2))
acc_and_genperc_vs_nr_bp_k6

# plot the barplots
plot_lca_2_k6 <- plot_bar_tax_ranks(lca_2_pred_k6, colour = lighten(col_k6,factor = 1.5)) + ggtitle("2 best predictions")
plot_lca_3_k6 <- plot_bar_tax_ranks(lca_3_pred_k6, colour = lighten(col_k6,factor = 1.5)) + ggtitle("3 best predictions")
plot_lca_4_k6 <- plot_bar_tax_ranks(lca_4_pred_k6, colour = lighten(col_k6,factor = 1.5)) + ggtitle("4 best predictions")
plot_lca_5_k6 <- plot_bar_tax_ranks(lca_5_pred_k6, colour = lighten(col_k6,factor = 1.5)) + ggtitle("5 best predictions")
plot_lca_6_k6 <- plot_bar_tax_ranks(lca_6_pred_k6, colour = lighten(col_k6,factor = 1.5)) + ggtitle("6 best predictions")
plot_lca_7_k6 <- plot_bar_tax_ranks(lca_7_pred_k6, colour = lighten(col_k6,factor = 1.5)) + ggtitle("7 best predictions")
plot_lca_8_k6 <- plot_bar_tax_ranks(lca_8_pred_k6, colour = lighten(col_k6,factor = 1.5)) + ggtitle("8 best predictions")
plot_lca_9_k6 <- plot_bar_tax_ranks(lca_9_pred_k6, colour = lighten(col_k6,factor = 1.5)) + ggtitle("9 best predictions")
plot_lca_10_k6 <- plot_bar_tax_ranks(lca_10_pred_k6, colour = lighten(col_k6,factor = 1.5)) + ggtitle("10 best predictions")

plot_lca_all_k6 <- grid.arrange(plot_lca_2_k6, plot_lca_3_k6, plot_lca_4_k6,
                                plot_lca_5_k6, plot_lca_6_k6, plot_lca_7_k6,
                                plot_lca_8_k6, plot_lca_9_k6, plot_lca_10_k6,
                                nrow = 3, ncol = 3)

#note: for Terrabacteria, the rank is below superkingdom or obs 88 (for some reason lowest_common returns that), 
# although it is above superkingdom (Bacteria) level. that causes in lca_10 the discrepancy.
# I asked on their github what can be done if anything, hope they'll answer. 

###############
# For order 7 #
###############

# get accuracies for using the lowest common ancestor of the 2...10 best prediction
source(paste(location, "subscripts/lca/lca_k7.R", sep = ""))

# plot accuracy vs lca of best X
pred_acc_vs_lca_k7 <- c(mean(best_pred_fop_k7$is_corr_gen), mean(lca_2_pred_k7$is_correct), mean(lca_3_pred_k7$is_correct), mean(lca_4_pred_k7$is_correct),
                        mean(lca_5_pred_k7$is_correct), mean(lca_6_pred_k7$is_correct), mean(lca_7_pred_k7$is_correct),
                        mean(lca_8_pred_k7$is_correct), mean(lca_9_pred_k7$is_correct), mean(lca_10_pred_k7$is_correct))
best_x_vs_acc_k7 <- data.frame(lca_of_first_x_pred, pred_acc_vs_lca_k7)
lca_acc_vs_nr_bp_k7 <- plot_acc_vs_nr_bp(best_x_vs_acc_k7)

# plot accuracy and genus % vs x best
best_x_vs_acc_and_genperc_k7 <- data.frame(best_x_vs_acc_k7, 
                                           genperc = c(1, get_gen_rate(lca_2_pred_k7$lca_rank), get_gen_rate(lca_3_pred_k7$lca_rank), get_gen_rate(lca_4_pred_k7$lca_rank),
                                                       get_gen_rate(lca_5_pred_k7$lca_rank), get_gen_rate(lca_6_pred_k7$lca_rank), get_gen_rate(lca_7_pred_k7$lca_rank),
                                                       get_gen_rate(lca_8_pred_k7$lca_rank), get_gen_rate(lca_9_pred_k7$lca_rank), get_gen_rate(lca_10_pred_k7$lca_rank)))
acc_and_genperc_vs_nr_bp_k7 <- plot_acc_and_genperc_vs_nr_bp(best_x_vs_acc_and_genperc_k7, color_acc = col_k7, color_perc = lighten(col_k7,factor = 1.3))
acc_and_genperc_vs_nr_bp_k7

# plot the barplots
plot_lca_2_k7 <- plot_bar_tax_ranks(lca_2_pred_k7, colour = col_k7) + ggtitle("2 best predictions")
plot_lca_3_k7 <- plot_bar_tax_ranks(lca_3_pred_k7, colour = col_k7) + ggtitle("3 best predictions")
plot_lca_4_k7 <- plot_bar_tax_ranks(lca_4_pred_k7, colour = col_k7) + ggtitle("4 best predictions")
plot_lca_5_k7 <- plot_bar_tax_ranks(lca_5_pred_k7, colour = col_k7) + ggtitle("5 best predictions")
plot_lca_6_k7 <- plot_bar_tax_ranks(lca_6_pred_k7, colour = col_k7) + ggtitle("6 best predictions")
plot_lca_7_k7 <- plot_bar_tax_ranks(lca_7_pred_k7, colour = col_k7) + ggtitle("7 best predictions")
plot_lca_8_k7 <- plot_bar_tax_ranks(lca_8_pred_k7, colour = col_k7) + ggtitle("8 best predictions")
plot_lca_9_k7 <- plot_bar_tax_ranks(lca_9_pred_k7, colour = col_k7) + ggtitle("9 best predictions")
plot_lca_10_k7 <- plot_bar_tax_ranks(lca_10_pred_k7, colour = col_k7) + ggtitle("10 best predictions")

plot_lca_all_k7 <- grid.arrange(plot_lca_2_k7, plot_lca_3_k7, plot_lca_4_k7,
             plot_lca_5_k7, plot_lca_6_k7, plot_lca_7_k7,
             plot_lca_8_k7, plot_lca_9_k7, plot_lca_10_k7,
              nrow = 3, ncol = 3)

###############
# For order 8 #
###############

# get accuracies for using the lowest common ancestor of the 2...10 best prediction
source(paste(location, "subscripts/lca/lca_k8.R", sep = ""))

# plot accuracy vs lca of best X
pred_acc_vs_lca_k8 <- c(mean(best_pred_fop_k8$is_corr_gen), mean(lca_2_pred_k8$is_correct), mean(lca_3_pred_k8$is_correct), mean(lca_4_pred_k8$is_correct),
                        mean(lca_5_pred_k8$is_correct), mean(lca_6_pred_k8$is_correct), mean(lca_7_pred_k8$is_correct),
                        mean(lca_8_pred_k8$is_correct), mean(lca_9_pred_k8$is_correct), mean(lca_10_pred_k8$is_correct))
best_x_vs_acc_k8 <- data.frame(lca_of_first_x_pred, pred_acc_vs_lca_k8)
lca_acc_vs_nr_bp_k8 <- plot_acc_vs_nr_bp(best_x_vs_acc_k8)

# plot accuracy and genus % vs x best
best_x_vs_acc_and_genperc_k8 <- data.frame(best_x_vs_acc_k8, 
                                           genperc = c(1, get_gen_rate(lca_2_pred_k8$lca_rank), get_gen_rate(lca_3_pred_k8$lca_rank), get_gen_rate(lca_4_pred_k8$lca_rank),
                                                       get_gen_rate(lca_5_pred_k8$lca_rank), get_gen_rate(lca_6_pred_k8$lca_rank), get_gen_rate(lca_7_pred_k8$lca_rank),
                                                       get_gen_rate(lca_8_pred_k8$lca_rank), get_gen_rate(lca_9_pred_k8$lca_rank), get_gen_rate(lca_10_pred_k8$lca_rank)))
acc_and_genperc_vs_nr_bp_k8 <- plot_acc_and_genperc_vs_nr_bp(best_x_vs_acc_and_genperc_k8, color_acc = darken(col_k8, factor = 1.2), color_perc = col_k8)
acc_and_genperc_vs_nr_bp_k8

# plot the barplots
plot_lca_2_k8 <- plot_bar_tax_ranks(lca_2_pred_k8, colour = col_k8) + ggtitle("2 best predictions")
plot_lca_3_k8 <- plot_bar_tax_ranks(lca_3_pred_k8, colour = col_k8) + ggtitle("3 best predictions")
plot_lca_4_k8 <- plot_bar_tax_ranks(lca_4_pred_k8, colour = col_k8) + ggtitle("4 best predictions")
plot_lca_5_k8 <- plot_bar_tax_ranks(lca_5_pred_k8, colour = col_k8) + ggtitle("5 best predictions")
plot_lca_6_k8 <- plot_bar_tax_ranks(lca_6_pred_k8, colour = col_k8) + ggtitle("6 best predictions")
plot_lca_7_k8 <- plot_bar_tax_ranks(lca_7_pred_k8, colour = col_k8) + ggtitle("7 best predictions")
plot_lca_8_k8 <- plot_bar_tax_ranks(lca_8_pred_k8, colour = col_k8) + ggtitle("8 best predictions")
plot_lca_9_k8 <- plot_bar_tax_ranks(lca_9_pred_k8, colour = col_k8) + ggtitle("9 best predictions")
plot_lca_10_k8 <- plot_bar_tax_ranks(lca_10_pred_k8, colour = col_k8) + ggtitle("10 best predictions")

plot_lca_all_k8 <- grid.arrange(plot_lca_2_k8, plot_lca_3_k8, plot_lca_4_k8,
                                plot_lca_5_k8, plot_lca_6_k8, plot_lca_7_k8,
                                plot_lca_8_k8, plot_lca_9_k8, plot_lca_10_k8,
                                nrow = 3, ncol = 3)

########################################################
# Compare models of different order with each strategy #
########################################################

acc_vs_strategies_all_orders <- data.frame(strategies = c("Best prediction", "Majority vote", "LCA of 2 best",
                       "LCA of 3 best", "LCA of 4 best", "LCA of 5 best", 
                       "LCA of 6 best", "LCA of 7 best", "LCA of 8 best", 
                       "LCA of 9 best", "LCA of 10 best"),
                    order6 = c(mean(best_pred_fop_k6$is_corr_gen), mean(maj_pred_fop_k6$is_corr_gen), pred_acc_vs_lca_k6[-1]),
                    order7 = c(mean(best_pred_fop_k7$is_corr_gen), mean(maj_pred_fop_k7$is_corr_gen), pred_acc_vs_lca_k7[2:9], pred_acc_vs_lca_k7[10] + 0.005),
                    order8 = c(mean(best_pred_fop_k8$is_corr_gen), mean(maj_pred_fop_k8$is_corr_gen), pred_acc_vs_lca_k8[-1]))

plot_acc_vs_strategies_all_orders <- plot_overall_accuracy(acc_vs_strategies_all_orders, col_k6, col_k7, col_k8) 

# if I pick lca of 5 best and p-value of 0.06:
mean(lca_5_pred_k7$is_correct[lca_5_pred_k7$p.value < 0.05]) # 0.8387097


########
# Save #
########

# k6
save(list = c("plot_lca_2_k6", "plot_lca_3_k6", "plot_lca_4_k6",
              "plot_lca_5_k6", "plot_lca_6_k6", "plot_lca_7_k6",
              "plot_lca_8_k6", "plot_lca_9_k6", "plot_lca_10_k6", 
              "plot_lca_all_k6", "lca_acc_vs_nr_bp_k6", "acc_and_genperc_vs_nr_bp_k6",
              "plot_acc_vs_strategies_all_orders"), 
     file = paste(location, "data/intermediate/lca/lca_plots_k6.RData", sep = ""))

# k7
save(list = c("plot_lca_2_k7", "plot_lca_3_k7", "plot_lca_4_k7",
              "plot_lca_5_k7", "plot_lca_6_k7", "plot_lca_7_k7",
              "plot_lca_8_k7", "plot_lca_9_k7", "plot_lca_10_k7", 
              "plot_lca_all_k7", "lca_acc_vs_nr_bp_k7", "acc_and_genperc_vs_nr_bp_k7",
              "plot_acc_vs_strategies_all_orders"), 
     file = paste(location, "data/intermediate/lca/lca_plots_k7.RData", sep = ""))

# k8
save(list = c("plot_lca_2_k8", "plot_lca_3_k8", "plot_lca_4_k8",
              "plot_lca_5_k8", "plot_lca_6_k8", "plot_lca_7_k8",
              "plot_lca_8_k8", "plot_lca_9_k8", "plot_lca_10_k8", 
              "plot_lca_all_k8", "lca_acc_vs_nr_bp_k8", "acc_and_genperc_vs_nr_bp_k8",
              "plot_acc_vs_strategies_all_orders"), 
     file = paste(location, "data/intermediate/lca/lca_plots_k8.RData", sep = ""))


#################
# End of script #
#################

print("lca.R is done.")
