# sim_wrp_cutoff() returns a named array, where the names are the cutoff values, the values are the number of cooccurences
# @param dat: data matrix with the columns being the binary vectors to compare to each other
# @param cutoff_range: numeric vector of cutoffs to try
# @return: named array, where the names are the cutoff values, the values are the number of cooccurences
get_sim_wrp_cutoff <- function(dat, cutoff_range){
  sim_wrp_cutoff <- sapply(cutoff_range, function(x) get_similarity_amount(dat, x))
  names(sim_wrp_cutoff) <- cutoff_range
  return(sim_wrp_cutoff)
}

# get_amount_of_similarity()
# @param dat:
# @param threshold:
# @return:
get_similarity_amount <- function(dat, threshold){
  similarity_amount <- sum(unlist(lapply(1:ncol(dat), function(y) length(unlist(get_similar(dat, y, threshold))))))/2
  return(similarity_amount)
}

# get_all_similar() given a dataframe and a threshold, returns a list, 
# where each element corresponds to the columns of the same index, and contains 
# the column names of other columns that are more identical to it than the threshold
# @param dat:
# @param threshold:
# @return:
get_all_similar <- function(dat, threshold){
  all_similar <- lapply(1:ncol(dat), function(y) unlist(get_similar(dat, y, threshold)))
  return(all_similar)
}

# get_all_similar() given a dataframe and a threshold, returns a list, 
# where each element corresponds to the columns of the same index, and contains 
# the column names of other columns that are more identical to it than the threshold
# @param dat:
# @param threshold:
# @return:
get_all_similar_100 <- function(dat){
  all_similar <- lapply(1:ncol(dat), function(y) unlist(get_similar_100(dat, y)))
  return(all_similar)
}

# get_similar() takes a data-frame, a column index and threshold, and returns the column name of those columns, that are more identical to col than the threshold
# @param dat:
# @param col:
# @param threshold:
# @return:
get_similar <- function(dat, col, threshold){
  cat("Column being analysed:", col, ", with threshold:", threshold, '\n')
  identity <- apply(dat, 2, function(x) sum(dat[,col] == x))
  take <- (identity  > threshold) & (names(identity ) != names(identity )[col])
  names_of_similar <- names(identity )[take]
  return(names_of_similar)
}


# get_similar() takes a data-frame, a column index and threshold, and returns the column name of those columns, that are more identical to col than the threshold
# @param dat:
# @param col:
# @param threshold:
# @return:
get_similar_100 <- function(dat, col){
  cat("Column being analysed:", col, ", with threshold: 100%", '\n')
  #same <- apply(dat[,-(col)], 2, function(x) (dat[,col] == x))
  same <- apply(dat[,], 2, function(x) (dat[,col] == x))
  identical <- apply(same, 2, function(x) (sum(x)==length(x)))
  names_of_similar_100 <- colnames(same)[identical]
  return(names_of_similar_100)
}



# return the prediction instead of contig name. WARNING: not customizable in terms of information containing columns! No error handling!
# @param node_names_list: a list of character vectors or lists containing node names
# @param predictions: dataframe with the predictions. Direct output of the final_analysis script of the wish_scripts module.
# @return: a list of character vectors of the predicted hosts, containing only observations that HAD a prediction. (for each element of the input list, will return an element lesser of equal to the input in size.)
host_list_for_cname_list <- function(node_names_list, predictions = pred_lca_5_p_full_k7){
  hosts_list <- (lapply(node_names_list, function(x) host_for_cname(x, predictions)))
  return(hosts_list)
}
# takes a column of a dataframe and checks its presence in a list of lists. returns a logical vector.
# @param df:
# @param list_of_lists:
# @param column:
are_in_lists <- function(df, list_of_lists, column = 1){
  query <- df[,column]
  results <- sapply(1:length(query), function(x) is_in_list(query[x], list_of_lists[[x]]))
  return(results)
}

jaccard_idx <- function(ls){
  dist <- unlist(lapply(combn(ls, 2, simplify = FALSE), function(x) {
    length(intersect(x[[1]], x[[2]]))/length(union(x[[1]], x[[2]])) }))
  return(dist)
}


##################################
# functions used by my functions #
##################################

# return the prediction instead of contig name. WARNING: not customizable in terms of information containing columns! No error handling!
# @param node_name: a character vector or list containing node names
# @param predictions: dataframe with the predictions. Direct output of the final_analysis script of the wish_scripts module.
# @return: a character vector of the predicted hosts, containing only observations that HAD a prediction. (will return a size lesser of equal to the input size)
host_for_cname <- function(node_names, predictions = pred_lca_5_p_full_k7){
  node_names2 <- unlist(node_names)
  hosts <- unlist(lapply(node_names2, function(x) predictions[predictions$contig == x,3]))
  return(hosts)
}

# get_lca() returns the lowest common ancestor of a given list of taxa
# @param taxa_list:
# @return: the lowest common ancestor name, rank and id, as present in ncbi's database
get_lca <- function (taxa_list) {
  if(length(taxa_list) == 1){
    return(c(name = as.character(taxa_list), rank = NA, id = NA))
  } else {
    return(tryCatch(unlist(lowest_common(taxa_list, db = "ncbi", rows = 2)), error=function(e) c(name = NA, rank = NA, id = NA)))
  }
}


#returns true or false depending on the queried element being in the queried list or not.
is_in_list <- function(element, list_of_elements){
  return(element %in% list_of_elements)
}

###################################
# phyloseq manipulating functions #
###################################

# returns the otu count vectors as a dataframe for each otu matching the given taxa in the list 
get_otus_by_taxa <- function(pseq, tax_list){
  # get rownames from the tax table
  tax_otu_name <- unlist(lapply(tax_list, function(x) dimnames(which((as(tax_table(pseq), "matrix")) == x, arr.ind = TRUE))[[1]]))
  # get the otu's accross samples
  otu_vectors <- data.frame(t(sapply(tax_otu_name, function(y) get_sample(pseq, y))))
  return(otu_vectors)
}




