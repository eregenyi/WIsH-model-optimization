#example data
n = 100
g = 6 
set.seed(g)
ex_dat <- data.frame(x = unlist(lapply(1:g, function(i) rnorm(n/g, runif(1)*i^2))), 
                y = unlist(lapply(1:g, function(i) rnorm(n/g, runif(1)*i^2))))
plot(ex_dat)
gap_statistic(ex_dat)

# load the real data
load(file = file.choose()) # load in.RData or final_analysis.RData

# gap statistic on real data 
g_gap_1 <- gap_statistic((ll_p_k7), min_num_clusters = 1, max_num_clusters = 10, num_reference_bootstraps = 10)
save(list = c("g_gap_1"), file = paste(location, "data/intermediate/final_analysis/gap_1.RData", sep = ""))
g_gap_2 <- gap_statistic(ll_p_k7, min_num_clusters = 11, max_num_clusters = 20, num_reference_bootstraps = 10)
save(list = c("g_gap_2"), file = paste(location, "data/intermediate/final_analysis/gap_2.RData", sep = ""))
g_gap_3 <- gap_statistic(ll_p_k7, min_num_clusters = 21, max_num_clusters = 30, num_reference_bootstraps = 10)
save(list = c("g_gap_3"), file = paste(location, "data/intermediate/final_analysis/gap_3.RData", sep = ""))
g_gap_4 <- gap_statistic(ll_p_k7, min_num_clusters = 31, max_num_clusters = 38, num_reference_bootstraps = 10)
save(list = c("g_gap_4"), file = paste(location, "data/intermediate/final_analysis/gap_4.RData", sep = ""))

df <- ggplot_build(g_gap_4)

plot_gap_statistic(g_gap_1[[1]],g_gap_1[[2]],g_gap_1[[3]])

gaps <- c(g_gap_1[[1]],g_gap_2[[1]],g_gap_3[[1]], ggplot_build(g_gap_4)$data[[2]]$y)
stddevs <- c(g_gap_1[[2]],g_gap_2[[2]],g_gap_3[[2]], (ggplot_build(g_gap_4)$data[[2]]$ymax - ggplot_build(g_gap_4)$data[[2]]$ymin))
num_clusters <- c(g_gap_1[[3]],g_gap_2[[3]],g_gap_3[[3]], ggplot_build(g_gap_4)$data[[2]]$x)

g_gap_criterion <- plot_gap_statistic(gaps[8:38],stddevs[8:38],num_clusters[8:38])

g_gap_criterion + 
  scale_x_continuous(breaks=seq(2, 38, 2)) + 
  labs(title = "Gap statistic", x = "Number of clusters", y = "Gap") +
  theme_gdocs() 
dim(ll_p_k7)

# save
save(list = c("gaps", "stddevs","num_clusters"),
     file = paste(location, "data/intermediate/final_analysis/gap.RData", sep = ""))
