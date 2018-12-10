

pcoa_dif_dist <- 
  lapply(getDistMethods()[-c(8,11,12,16,18,19:22,25:28,32,37:38,40:41)], 
         function(x) tryCatch(pcoa(distance(t(ll_p_k7), method = x, p = 1), correction="none", rn=NULL),error=function(e) NA))
  
class(pcoa_dif_dist)
length(pcoa_dif_dist)

str(pcoa_dif_dist[1])
pcoa_dif_dist[[1]]$values

################
# Plot results #
################

## 2D plots 

plot(y = pcoa_dif_dist[[1]]$vectors[,2], x = pcoa_dif_dist[[1]]$vectors[,1])
plot(y = pcoa_dif_dist[[2]]$vectors[,2], x = pcoa_dif_dist[[2]]$vectors[,1])
plot(y = pcoa_dif_dist[[3]]$vectors[,2], x = pcoa_dif_dist[[3]]$vectors[,1])

## 3D plots

plot3d(pcoa_dif_dist[[1]]$vectors[,1], pcoa_dif_dist[[1]]$vectors[,2],pcoa_dif_dist[[1]]$vectors[,3])



########
# Save #
########

save(list = c("pcoa_dif_dist"),
     file = paste(location, "data/intermediate/final_analysis/pcoa.RData", sep = ""))
