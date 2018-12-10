
print("\t Running lca_k8.R...")

##########
# For k8 #
##########

#lca of first 2
lca_2_pred_k8 <- merge(x = data.frame(pred_top_10_k8[,1:2], merge_into_lca(pred_top_10_k8, col_taxa = 3:4, id = FALSE)), 
                       y = blastn_host, by.x = "contig", by.y = "contig")
lca_2_pred_k8$is_correct <- is_ancestor(lca_2_pred_k8, col_ancestor = 3, col_target = 5)

#lca of first 3
lca_3_pred_k8 <- merge(x = data.frame(pred_top_10_k8[,1:2], merge_into_lca(pred_top_10_k8, col_taxa = 3:5, id = FALSE)), 
                       y = blastn_host, by.x = "contig", by.y = "contig")
lca_3_pred_k8$is_correct <- is_ancestor(lca_3_pred_k8, col_ancestor = 3, col_target = 5)

#lca of first 4
lca_4_pred_k8 <- merge(x = data.frame(pred_top_10_k8[,1:2], merge_into_lca(pred_top_10_k8, col_taxa = 3:6, id = FALSE)), 
                       y = blastn_host, by.x = "contig", by.y = "contig")
lca_4_pred_k8$is_correct <- is_ancestor(lca_4_pred_k8, col_ancestor = 3, col_target = 5)

#lca of first 5 
lca_5_pred_k8 <- merge(x = data.frame(pred_top_10_k8[,1:2], merge_into_lca(pred_top_10_k8, col_taxa = 3:7, id = FALSE)), 
                       y = blastn_host, by.x = "contig", by.y = "contig")
lca_5_pred_k8$is_correct <- is_ancestor(lca_5_pred_k8, col_ancestor = 3, col_target = 5)

#lca of first 6
lca_6_pred_k8 <- merge(x = data.frame(pred_top_10_k8[,1:2], merge_into_lca(pred_top_10_k8, col_taxa = 3:8, id = FALSE)), 
                       y = blastn_host, by.x = "contig", by.y = "contig")
lca_6_pred_k8$is_correct <- is_ancestor(lca_6_pred_k8, col_ancestor = 3, col_target = 5)

#lca of first 7
lca_7_pred_k8 <- merge(x = data.frame(pred_top_10_k8[,1:2], merge_into_lca(pred_top_10_k8, col_taxa = 3:9, id = FALSE)), 
                       y = blastn_host, by.x = "contig", by.y = "contig")
lca_7_pred_k8$is_correct <- is_ancestor(lca_7_pred_k8, col_ancestor = 3, col_target = 5)

#lca of first 8 
lca_8_pred_k8 <- merge(x = data.frame(pred_top_10_k8[,1:2], merge_into_lca(pred_top_10_k8, col_taxa = 3:10, id = FALSE)), 
                       y = blastn_host, by.x = "contig", by.y = "contig")
lca_8_pred_k8$is_correct <- is_ancestor(lca_8_pred_k8, col_ancestor = 3, col_target = 5)

#lca of first 9
lca_9_pred_k8 <- merge(x = data.frame(pred_top_10_k8[,1:2], merge_into_lca(pred_top_10_k8, col_taxa = 3:11, id = FALSE)), 
                       y = blastn_host, by.x = "contig", by.y = "contig")
lca_9_pred_k8$is_correct <- is_ancestor(lca_9_pred_k8, col_ancestor = 3, col_target = 5)

#lca of first 10
lca_10_pred_k8 <- merge(x = data.frame(pred_top_10_k8[,1:2], merge_into_lca(pred_top_10_k8, col_taxa = 3:12, id = FALSE)), 
                        y = blastn_host, by.x = "contig", by.y = "contig")
lca_10_pred_k8$is_correct <- is_ancestor(lca_10_pred_k8, col_ancestor = 3, col_target = 5)

################
# Save results #
################

save(list = c("lca_2_pred_k8", "lca_3_pred_k8", "lca_4_pred_k8",
              "lca_5_pred_k8", "lca_6_pred_k8", "lca_7_pred_k8",
              "lca_8_pred_k8", "lca_9_pred_k8", "lca_10_pred_k8"), 
     file = paste(location, "data/intermediate/lca/lca_k8.RData", sep = ""))


#################
# End of script #
#################

print("\t lca_k8.R is finished.")