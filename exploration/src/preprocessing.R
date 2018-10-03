
print("preprocessing.R is running...")

#########################
# Parse into binary OTU #
#########################

vir_ps_co <- vir_ps
binary_otu <- ifelse(as.data.frame(as(otu_table(vir_ps_co), "matrix")) == 0, 0, 1)
otu_table(vir_ps_co) <- otu_table(binary_otu, taxa_are_rows = FALSE)

#########################
# Co-occurence analysis #
#########################

similarity_list <- get_all_similar(binary_otu, threshold = 300)

get_all_similar(binary_otu[,1:600], threshold = 300)
# how many similarities are detected?
indices <- c()
for (i in 1:length(similarity_list)){  if (length(similarity_list[[i]])>0){    indices <- c(indices, i)  }}

# Check out the contigs
length(similarity_list[[indices[6]]])
similarity_list[[indices[6]]]
colnames(binary_otu)[566]

tax <- data.frame(as(tax_table(vir_ps_co), "matrix"))
tax$contig <- row.names(tax)
View(tax)

preds <- na.omit(pred_lca_5_full_k7)
tax <- merge(tax, preds, by.x = "contig", by.y = "contig")

#############################
# Save intermediate results #
#############################

save(list = c("vir_ps_co", "binary_otu", "similarity_list", "tax", "similarity_list"), 
     file = paste(location, "data/intermediate/preprocessing/preprocessing.RData", sep = ""))

#################
# End of script #
#################

print("preprocessing.R is done.")

