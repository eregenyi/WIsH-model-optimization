
print("co_occurence_viral.R is running...")

###########################
# Determine the threshold #
###########################

# get similarity amount (number of co-occurences) with respect to the cutoff value
gswc <- get_sim_wrp_cutoff(binary_otu, 295:305)

# Plot number of co-occurences w.r.p. cutoff value from 295-305
plot(y = gswc, x = names(gswc), ty = 'l', xlab = "Similarity threshold", ylab = "Number of co-occurences")

#####################
# Check predictions #
#####################

# get contigs that are more identical than 302
co_occurence_list <- get_all_similar(binary_otu, 302)
co_occuring_seqs_idx <- unlist(lapply(1:length(co_occurence_list), function(x) length(co_occurence_list[[x]]) != 0L))
co_occuring_seqs <- co_occurence_list[co_occuring_seqs_idx]

# maybe duplicate removal is needed at this stage!
# get host list & remove the entries from the list, where no predictions were made
co_occuring_seq_hosts_list <- host_list_for_cname_list(co_occuring_seqs)
co_occuring_seq_host_idx <- unlist(lapply(1:length(co_occuring_seq_hosts_list), function(x) length(co_occuring_seq_hosts_list[[x]]) != 0L))
co_occuring_seq_hosts <- co_occuring_seq_hosts_list[co_occuring_seq_host_idx]

# get the lowest common ancestor of the predictions, and check if the prediction list contained it. If yes, the prediction is deemed correct.
lca_co_occuring_seq_hosts <- sapply(co_occuring_seq_hosts, function(x) get_lca(unlist(x)))
lca_co_occuring_seq_hosts_df <- as.data.frame(t(lca_co_occuring_seq_hosts))
lca_co_occuring_seq_hosts_df$is_pred_host_the_same <- as.numeric(are_in_lists(lca_co_occuring_seq_hosts_df, co_occuring_seq_hosts))

# how many of the predictions were "correct" ?
View(lca_co_occuring_seq_hosts_df)
lca_co_occuring_seq_hosts_df_no_na <- na.omit(lca_co_occuring_seq_hosts_df)
View(lca_co_occuring_seq_hosts_df_no_na)
sum(lca_co_occuring_seq_hosts_df_no_na$is_pred_host_the_same)/nrow(lca_co_occuring_seq_hosts_df_no_na)
#84.2% of the cooccuring sequences have the same predicted hosts. PROBLEM: many are superkingdom level

#############################
# Require 100% co-occurence #
#############################

# get contigs names into a list that are perfectly identical
co_occurence_list_100 <- get_all_similar_100(binary_otu)
co_occuring_seqs_idx_100 <- unlist(lapply(1:length(co_occurence_list_100), function(x) length(co_occurence_list_100[[x]]) != 0L))
co_occuring_seqs_100 <- co_occurence_list_100[co_occuring_seqs_idx_100]

# remove duplicates & elements which only contain one node
is_co_occuring_seqs_100 <- unlist(lapply(1:length(co_occuring_seqs_100), function(x) length(co_occuring_seqs_100[[x]]) != 1))
co_occuring_seqs_100 <- co_occuring_seqs_100[is_co_occuring_seqs_100]
co_occuring_seqs_100 <- co_occuring_seqs_100[!duplicated(lapply(co_occuring_seqs_100, sort))]

# get host list & remove the entries from the list, where no predictions were made
co_occuring_seq_hosts_list_100 <- host_list_for_cname_list(co_occuring_seqs_100)
co_occuring_seq_host_idx_100 <- unlist(lapply(1:length(co_occuring_seq_hosts_list_100), function(x) length(co_occuring_seq_hosts_list_100[[x]]) != 0L))
co_occuring_seq_hosts_100 <- co_occuring_seq_hosts_list_100[co_occuring_seq_host_idx_100]
co_occuring_seq_host_idx2_100 <- unlist(lapply(1:length(co_occuring_seq_hosts_100), function(x) length(co_occuring_seq_hosts_100[[x]]) != 1))
co_occuring_seq_hosts_100 <- co_occuring_seq_hosts_100 [co_occuring_seq_host_idx2_100]

# get the lowest common ancestor of the predictions, and check if the prediction list contained it. If yes, the prediction is deemed correct.
lca_co_occuring_seq_hosts_100 <- sapply(co_occuring_seq_hosts_100, function(x) get_lca(unlist(x)))
lca_co_occuring_seq_hosts_df_100 <- as.data.frame(t(lca_co_occuring_seq_hosts_100))
lca_co_occuring_seq_hosts_df_100$is_pred_host_the_same <- as.numeric(are_in_lists(lca_co_occuring_seq_hosts_df_100, co_occuring_seq_hosts_100))

# how many of the predictions were "correct" ?
View(lca_co_occuring_seq_hosts_df_100)
lca_co_occuring_seq_hosts_df_no_na_100 <- na.omit(lca_co_occuring_seq_hosts_df_100)
View(lca_co_occuring_seq_hosts_df_no_na_100)
sum(lca_co_occuring_seq_hosts_df_no_na_100$is_pred_host_the_same)/nrow(lca_co_occuring_seq_hosts_df_no_na_100)
#76%

# save data 
node_list_100 <- co_occuring_seqs_100[co_occuring_seq_host_idx_100][co_occuring_seq_host_idx2_100]
host_list_100 <- co_occuring_seq_hosts_100
correctness_100 <- lca_co_occuring_seq_hosts_df_no_na_100

write.xlsx(x = correctness_100, file = paste(location, "data/output/correctness_100.xlsx", sep = ""))
lapply(lapply(host_list_100, as.character), write, file = paste(location, "data/output/host_list_100.txt", sep = ""), append=TRUE, ncolumns=1000)
lapply(lapply(node_list_100, as.character), write, file = paste(location, "data/output/node_list_100.txt", sep = ""), append=TRUE, ncolumns=1000)
lapply(lapply(co_occuring_seqs, as.character), write, file = paste(location, "data/output/co_occ_list_gt_302.txt", sep = ""), append=TRUE, ncolumns=1000)

########
# Save #
########

save(list = c("gswc", "co_occurence_list", "lca_co_occuring_seq_hosts", "lca_co_occuring_seq_hosts_df",
              "co_occurence_list_100", "lca_co_occuring_seq_hosts_100"), 
     file = paste(location, "data/intermediate/co_vir/co_vir.RData", sep = ""))


#################
# End of script #
#################

print("co_occurence_viral.R is done.")
