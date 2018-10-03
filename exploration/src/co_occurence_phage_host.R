"
Get genus level predictions, and the corresponding viral abundance vector. 
Get the corresponding prokaryotic genera abundance vectors.
Plot them for each baby onto the same plot.
"


# Try things for Salmonella for instance
View(pred_lca_5_p_full_k7)
View(pred_lca_5_p_full_k7[pred_lca_5_p_full_k7$lca_rank == "genus",])
unique(pred_lca_5_p_full_k7[pred_lca_5_p_full_k7$lca_rank == "genus",3])
nrow(pred_lca_5_p_full_k7[pred_lca_5_p_full_k7$lca_name == "Salmonella" ,])

otu_table(bac_ps)$row.names(tax_table(bac_ps)$Genus == "Salmonella")

View(tax_table(bac_ps))
View(sample_data(bac_ps))
View(as(bac_ps, "matrix"))

salmonella <-dimnames(which((as(tax_table(bac_ps), "matrix")) == "Clostridium", arr.ind = TRUE))[[1]]
salmonella2 <- names(which((as(tax_table(bac_ps), "matrix"))[,6] == "Salmonella"))
taxa_names(bac_ps)
sample_names(bac_ps)
which((get_sample(bac_ps, salmonella)) !=0)

which((get_otus_by_taxa(bac_ps, c("Salmonella", "Pseudomonas"))) != 0 )
