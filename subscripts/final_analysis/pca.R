# PCA with scaling and centering

pca_ll_k7 <- prcomp(ll_k7,
                    center = TRUE,
                    scale. = TRUE)
pca_ll_p_k7 <- prcomp(ll_p_k7,
                    center = TRUE,
                    scale. = TRUE)
pca_ll_p_not_scaled_k7 <- prcomp(t(ll_p_k7),
                      center = TRUE,
                      scale. = FALSE)

# scree plot
plot(pca_ll_k7, type = "l")
plot(pca_ll_p_k7, type = "l")
# scores wrp the first 2 PCs
plot(x = pca_ll_k7$x[,1], y = pca_ll_k7$x[,2])
plot(x = pca_ll_p_k7$x[,1], y = pca_ll_p_k7$x[,2])
# importance of PCs
summary(pca_ll_k7)$importance[,0:10]
summary(pca_ll_p_not_scaled_k7)$importance[,0:10]
# 3D plots
plot3d(pca_ll_k7$x[,1], pca_ll_k7$x[,2], pca_ll_k7$x[,3])
plot3d(pca_ll_p_k7$x[,1], pca_ll_p_k7$x[,2], pca_ll_p_k7$x[,3])
plot3d(pca_ll_p_not_scaled_k7$x[,1], pca_ll_p_not_scaled_k7$x[,2], pca_ll_p_not_scaled_k7$x[,3])

# biplots - based on these, scalng may be cool actually. 
j <- ggbiplot(pca_ll_p_k7, obs.scale = 1, var.scale = 1, circle = TRUE) + 
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
print(j)

