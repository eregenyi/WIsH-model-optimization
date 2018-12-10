
print("import_data.R is running...")

###############
# import data #
###############

#wish output for different order models
ll_k6 <- read.table(file = paste(location, "data/input/output_kegg_k6/llikelihood.matrix", sep = "" ), header = TRUE) #llikelihood.matrix for order6
predictions_k6 <- read.table(file = paste(location, "data/input/output_kegg_k6/prediction.list", sep = ""), header = TRUE) #prediction.list order6

ll_k7 <- read.table(file = paste(location, "data/input/output_kegg_k7/llikelihood.matrix", sep = ""), header = TRUE) #llikelihood.matrix for order7
predictions_k7 <- read.table(file = paste(location, "data/input/output_kegg_k7/prediction.list", sep = ""), header = TRUE) #prediction.list order7

ll_k8 <- read.table(file = paste(location, "data/input/output_kegg_k8/llikelihood.matrix", sep = ""), header = TRUE) #llikelihood.matrix for order8
predictions_k8 <- read.table(file = paste(location, "data/input/output_kegg_k8/prediction.list", sep = ""), header = TRUE) #prediction.list order8

#metadata: wish identifies contigs with the filenames. These files contain the filenames and the corresponding contig names in case of phages, and the filenames, accession numbers (ncbi), and taxonomy in case of prokaryotes
phages_list <- read.table(file = paste(location, "data/input/meta/phage_genome_meta_data", sep = ""), header = FALSE) #phage_genomes_meta_data, contains the file names and the phage names (NODE)
prok_list <-  read_delim(file = paste(location, "data/input/meta/prok_genomes_meta_data", sep = ""), "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE) #prok_genomes_meta_data, contains the prokaryote genome file names and the prokaryote names

#other input data files: for model validation and final assession
test_phage_clustering <- read.table(file = paste(location, "data/input/meta/blastn_annotated_phage.txt", sep = ""), header = TRUE)
blastn_host <- read.table(file = paste(location, "data/input/meta/blastn_host.txt", sep = ""), header = TRUE)
kt_family <- read.table(file = paste(location, "data/input/meta/kt_blastn.txt", sep = ""), header = TRUE)

#16S data corresponding to the samples
load(file = paste(location, "data/input/meta/bac_ps.RData", sep = ""))

#################
# End of script #
#################

print("import_data.R is done.")