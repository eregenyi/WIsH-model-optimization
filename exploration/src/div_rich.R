#########
# Viral #
#########

#append the meta table with diversity
sample_data(vir_ps)$diversity_contig <- diversity(otu_table(vir_ps), index = "shannon")

#append the meta table with richness
sample_data(vir_ps)$richness_contig <- specnumber(otu_table(vir_ps))

#plot diversity
ggplot(sample_data(vir_ps), aes(x = X.days, y = diversity_contig)) + 
  geom_point() + facet_wrap(~baby, ncol = 4) + 
  geom_smooth(method = "loess", span = 0.3) + 
  theme_gdocs()+
  xlab("Time [days after birth]") + ylab("Viral diversity (Shannon)") + ggtitle("Diversity (contig)")

#plot richness
ggplot(sample_data(vir_ps), aes(x = X.days, y = richness_contig)) + 
  geom_point() + facet_wrap(~baby, ncol = 4) + 
  geom_smooth(method = "loess", span = 0.3, colour = "darkblue") + 
  theme_gdocs()+
  ylab("Viral richness") + xlab("Time [days after birth]") + ggtitle("Richness (contig)")

###############
# Just phages #
###############

vir_ps_p <- subset_taxa(vir_ps, final_virus_type == "Phage")
vir_ps_p <- prune_taxa(taxa_sums(vir_ps_p) > 0, vir_ps_p)

#append the meta table with diversity
sample_data(vir_ps_p)$diversity_contig <- diversity(otu_table(vir_ps_p), index = "shannon")

#append the meta table with richness
sample_data(vir_ps_p)$richness_contig <- specnumber(otu_table(vir_ps_p))

#plot diversity
ggplot(sample_data(vir_ps_p), aes(x = X.days, y = diversity_contig)) + 
  geom_point() + facet_wrap(~baby, ncol = 4) + 
  geom_smooth(method = "loess", span = 0.3, colour = "pink3") + 
  theme_gdocs()+
  xlab("Time [days after birth]") + ylab("Phage diversity (Shannon)") + ggtitle("Diversity (contig), phages")

#plot richness
ggplot(sample_data(vir_ps_p), aes(x = X.days, y = richness_contig)) + 
  geom_point() + facet_wrap(~baby, ncol = 4) + 
  geom_smooth(method = "loess", span = 0.3, colour = "darkred") + 
  theme_gdocs()+
  ylab("Phage richness") + xlab("Time [days after birth]") + ggtitle("Richness (contig), phages")

###########################
# Just eukaryotic viruses #
###########################

vir_ps_e <- subset_taxa(vir_ps, final_virus_type == "Virus_eukaryotic")
vir_ps_e <- prune_taxa(taxa_sums(vir_ps_e) > 0, vir_ps_e)

#append the meta table with diversity
sample_data(vir_ps_e)$diversity_contig <- diversity(otu_table(vir_ps_e), index = "shannon")

#append the meta table with richness
sample_data(vir_ps_e)$richness_contig <- specnumber(otu_table(vir_ps_e))

#plot diversity
ggplot(sample_data(vir_ps_e), aes(x = X.days, y = diversity_contig)) + 
  geom_point() + facet_wrap(~baby, ncol = 4) + 
  geom_smooth(method = "loess", span = 0.3, colour = "orange3") + 
  theme_gdocs()+
  xlab("Time [days after birth]") + ylab("Eukaryotic viral diversity (Shannon)") + ggtitle("Diversity (contig), eukaryotic virus")

#plot richness
ggplot(sample_data(vir_ps_e), aes(x = X.days, y = richness_contig)) + 
  geom_point() + facet_wrap(~baby, ncol = 4) + 
  geom_smooth(method = "loess", span = 0.3, colour = "orange4") + 
  theme_gdocs()+
  ylab("Eukaryotic viral richness") + xlab("Time [days after birth]") + ggtitle("Richness (contig), eukaryotic virus")



###############
# Prokaryotic #
###############

sample_data(bac_ps)$unique_ID <- gsub("_", ".", sample_data(bac_ps)$unique_ID, fixed=TRUE)
merged <- merge(y = as.data.frame(as(sample_data(vir_ps), "matrix")), 
                             x = as.data.frame(as(sample_data(bac_ps), "matrix")), 
                             by.y = "Sample.ID", by.x = "unique_ID", all = TRUE)
row.names(merged) <- merged$ID
sample_data(bac_ps) <- merged
bac_ps <- subset_samples(bac_ps, !is.na(baby))

#append the meta table with diversity
sample_data(bac_ps)$diversity_contig <- diversity(otu_table(bac_ps), index = "shannon")

#append the meta table with richness
sample_data(bac_ps)$richness_contig <- specnumber(otu_table(bac_ps))

bac_plot_data <- na.omit(data.frame(X.days = sample_data(bac_ps)$X.days,
                                    diversity_contig = sample_data(bac_ps)$diversity_contig,
                                    richness_contig = sample_data(bac_ps)$richness_contig,
                                    baby = sample_data(bac_ps)$baby,
                                    nr = sample_data(bac_ps)$nr,
                             ID = sample_data(bac_ps)$ID))

bac_plot_data <- bac_plot_data[order(bac_plot_data$nr),]
str(bac_plot_data)
bac_plot_data$X.days <- as.numeric(bac_plot_data$X.days)
bac_plot_data$richness_contig <- as.numeric(bac_plot_data$richness_contig)

#plot diversity
ggplot(bac_plot_data , aes(x = X.days, y = diversity_contig)) + 
  geom_point() + facet_wrap(~baby, ncol = 4) + 
  geom_smooth(method = "loess", span = 0.3, colour = "green3") + 
  theme_gdocs()+
  xlab("Time [days after birth]") + ylab("Prokaryotic diversity (Shannon)") + ggtitle("Diversity (contig)")

#plot richness
ggplot(bac_plot_data, aes(x = X.days, y = richness_contig)) + 
  geom_point() + facet_wrap(~baby, ncol = 4) + 
  geom_smooth(method = "loess", span = 0.3, colour = "darkgreen") + 
  theme_gdocs()+
  ylab("Prokaryotic richness") + xlab("Time [days after birth]") + ggtitle("Richness (contig)")

View(unique(tax_table(vir_ps)))
