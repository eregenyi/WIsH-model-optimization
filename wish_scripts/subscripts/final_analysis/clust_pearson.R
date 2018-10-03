
cat('\014')

##############
# Clustering #
##############

# script variables
k_range <- 2:26
dist_mx <- 1 - cor(ll_p_k7) #here cor computes distances differently. we should not transpose the data matrix
dist_ll <- as.dist(dist_mx)

################
# 1. Dendogram #      (Via hierarchical clustering)
################

d_hc <- hclust(dist_ll, method = "complete")
plot(d_hc) #4 seems to be a sensible choice

ggdendrogram(d_hc, rotate = FALSE, size = 2, labels = FALSE) +
  labs(title = 
         "Hierarchical clustering", 
       x = "Phage contigs", y = "") +
  theme_gdocs() +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank())

#################
# 2. Elbow plot #
#################

d_kmc <- lapply(k_range, function(x) kmeans(dist_ll, x, 100))

mean_wws <- lapply(1:length(d_kmc), function(x) mean(d_kmc[[x]]$withinss))
tot_wws <- lapply(1:length(d_kmc), function(x) d_kmc[[x]]$tot.withinss)

m_wws <- data.frame(k_range, mean_wws = unlist(mean_wws))
t_wws <- data.frame(k_range, tot_wws = unlist(tot_wws))

plot(k_range, mean_wws, type = 'o')
plot(k_range, tot_wws, type = 'o')

col_k7 <- viridis(3)[2]

# mean withinss - withinss just means the sum of squared error within classes
ggplot(data = m_wws, aes(x = k_range, y = mean_wws)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks=seq(2,38,2)) + 
  labs(title = "Mean within SS", x = "Number of clusters", y = "Mean withinness") +
  theme_gdocs() +
  geom_vline(xintercept = 38, col = col_k7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
      axis.text.x=element_blank(), axis.text.y=element_blank())


# total withinness
ggplot(data = t_wws, aes(x = k_range, y = tot_wws)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks=seq(2,38,2)) + 
  labs(title = "Total within SS", x = "Number of clusters", y = "Total withinness") +
  theme_gdocs() +
  geom_vline(xintercept = 38, col = col_k7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x=element_blank())

# suggests 3-4 classes?
# change k_range

k_range <- 2:38

#################
# 3. Silhouette #
#################

d_pam <- lapply(k_range, function(x) pam(dist_ll, diss = TRUE,  x))
d_pam_sil <- unlist(lapply(1:length(k_range), function(x) d_pam[[x]]$silinfo$avg.width))
plot(y = c(d_pam_sil_2, d_pam_sil), x = k_range, type = 'ol') # optimal nr of clusters = 4

sil_crit <- data.frame(k_range = k_range, sil = d_pam_sil)

ggplot(data = sil_crit, aes(x = k_range, y = sil)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks=seq(2,38,2)) + 
  labs(title = "Silhouette criterion", x = "Number of clusters", y = "Average silhouette width") +
  theme_gdocs() 

which.max(d_pam_sil) 
#so 2 clusters? does not give much

##################################
# 4. Calinski-Harabasz criterion #
##################################
# For k=1:10, with 100 iterations, this one runs for appr. 2 days for data of size of 5000x2000

d_cascade_km <- cascadeKM(scale(dist_mx, center = TRUE,  scale = FALSE), 2, 38, iter = 100)

plot(d_cascade_km, sortg = TRUE, grpmts.plot = TRUE)
d_cascade_km
d_cascade_km$results
calinski.best <- as.numeric(which.max(d_cascade_km$results[2,]))
cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")
plot(c(d_cascade_km$results[2,],d_cascade_km_2$results[2,]), type = "ol")

ch_crit <- data.frame(k_range = k_range, 
                      ch_idx = d_cascade_km$results[2,])

ggplot(data = ch_crit, aes(x = k_range, y = ch_idx)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks=seq(2,38,2)) + 
  labs(title = "Celinski-Harabasz criterion", x = "Number of clusters", y = "Celinski index") +
  theme_gdocs() 

#########################
# 5. Gap Statistic (GS) #
#########################

#clusGap(dist_mx, kmeans, 10, B = 100, verbose = interactive())
#graphical representation

g_gap_1 <- gap_statistic((ll_p_k7), min_num_clusters = 1, max_num_clusters = 10, num_reference_bootstraps = 10)
save(list = c("g_gap_1"), file = paste(location, "data/intermediate/final_analysis/gap_1.RData", sep = ""))
g_gap_2 <- gap_statistic(ll_p_k7, min_num_clusters = 11, max_num_clusters = 20, num_reference_bootstraps = 10)
save(list = c("g_gap_2"), file = paste(location, "data/intermediate/final_analysis/gap_2.RData", sep = ""))
g_gap_3 <- gap_statistic(ll_p_k7, min_num_clusters = 21, max_num_clusters = 30, num_reference_bootstraps = 10)
save(list = c("g_gap_3"), file = paste(location, "data/intermediate/final_analysis/gap_3.RData", sep = ""))
g_gap_4 <- gap_statistic(ll_p_k7, min_num_clusters = 31, max_num_clusters = 38, num_reference_bootstraps = 10)
save(list = c("g_gap_4"), file = paste(location, "data/intermediate/final_analysis/gap_4.RData", sep = ""))

gaps <- c(g_gap_1[[1]],g_gap_2[[1]],g_gap_3[[1]], g_gap_4[[1]])
stddevs <- c(g_gap_1[[2]],g_gap_2[[2]],g_gap_3[[2]], g_gap_4[[2]])
num_clusters <- c(g_gap_1[[3]],g_gap_2[[3]],g_gap_3[[3]], g_gap_4[[3]])

g_gap_criterion <- plot_gap_statistic(gaps[8:38],stddevs[8:38],num_clusters[8:38])

g_gap_criterion + 
  scale_x_continuous(breaks=seq(2, 38, 2)) + 
  labs(title = "Gap statistic", x = "Number of clusters", y = "Gap") +
  theme_gdocs() 

plot(c(ggplot_build(g_gap_crit)$data[[2]]$y, ggplot_build(g_gap_crit_2)$data[[2]]$y))

###############
# Final plots #
###############

pcoa_cor <- pcoa(dist_mx, correction="none", rn=NULL)
str(d_pam[[19]]$clustering) #the range was 2:269, so the 16th element of the list will correspond to k=18
str(pcoa_cor)
biplot(pcoa_cor)

p_col <- plot(pcoa_cor$vectors[,1:2], col = d_pam[[16]]$clustering)
p <- plot(pcoa_cor$vectors[,1:2])
p3_col <- plot3d(pcoa_cor$vectors[,1:3], col = d_pam[[16]]$clustering)
p3 <- plot3d(pcoa_cor$vectors[,1:3])

plot(pcoa_cor$values$Eigenvalues[1:10])

rgl.postscript( filename = paste(location, "figures/final_analysis/pcoa_3d_col_k=18.pdf", sep = ""),"pdf")

##################
# Check clusters #
##################

#load final_analysis.RData, as it contains the reference dataset, test_phage_clustering
phage_clusters <- data.frame(contig = names(d_pam[[19]]$clustering), cluster = d_pam[[19]]$clustering)
phage_clusters_and_true_taxa <- merge(x = test_phage_clustering, y = phage_clusters, by.x = "contig", by.y = "contig")
ph_ctr <- merge(x = phage_clusters_and_true_taxa, y = pred_lca_5_p_full_k7[,c(1,4)], by.x = "contig", by.y = "contig")
View(ph_ctr[ph_ctr$lca_rank == "genus",])
View(ph_ctr)
# for my latex file:

contig_taxa_host <- merge(x = test_phage_clustering, y = blastn_host, by.x = "contig", by.y = "contig")
contig_taxa_host_cluster <- merge(x = contig_taxa_host, y = phage_clusters, by.x = "contig", by.y = "contig", all.x = TRUE)

########
# Save #
########

save(list = c("ll_p_k7", "k_range", "dist_mx", "dist_ll", "d_hc",
              "d_kmc", "mean_wws", "tot_wws", "m_wws", "t_wws",
              "d_pam", "d_pam_sil", 
              "d_cascade_km", 
              "g_gap_crit", "g_gap_crit_2",
              "d_pam_2", "d_cascade_km_2", "pcoa_cor",
              "g_gap_criterion", "num_clusters", "stddevs", "gaps"),
     file = paste(location, "data/intermediate/final_analysis/clust_pearson.RData", sep = ""))
