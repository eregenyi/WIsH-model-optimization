######################
# Plotting functions #
######################

# plot_heatmap_and_sample()
# @param ll:
# @param sample_row:
# @param sample_col:
# @param seed: a seed for the random number generator. For seed = 1, the same set of random numbers will be generated always.
# @param sample_size:
# @return: a graphical, heatmap type object
plot_hmp<- function(likelihood_matrix, sample_row = TRUE, sample_col = TRUE, seed = 1, 
                                    sample_size_row = 80, sample_size_col = 1500, prok_l = prok_list, col_fname = 1, col_gen = 3){
  
  ll <- sample_matrix(likelihood_matrix, sample_row, sample_col, seed, sample_size_row, sample_size_col)
  genera <- get_host_names(vec = row.names(ll), prok_l, col_fname, col_gen)
  ll <- as.matrix(ll)
  return(heatmap(ll, labRow=genera, labCol = FALSE))
}

# plot_pie_of_preds() plots a pie chart creates a pie chart representing the 
# @param df: a dataframe, whose one specific column is to be depicted
# @param column: the specific column to analyse
# @param lab: label of the plot
# @param trsh: treshold value. Everything with a lower percentage than the threshold will be summed up under "Other"
# @post: draws a piecart with the percentage of occurences considered as proportions
# @return: void
plot_pie_of_preds <- function(df, column = 3, lab = "Pie Chart", trsh = 2){
  pie_data <- sum_perc_of_unique(df[,column])
  pie_data <- rbind(pie_data[pie_data[,2] > trsh,], data.frame(unique_elements = "Other", perc_of_occ = sum(pie_data[pie_data[,2] <= trsh,2]), sum_of_occ = sum(pie_data[pie_data[,2] <= trsh,3])))
  slices <- pie_data[,2] 
  lbls <- pie_data[,1]
  lbls <- paste(lbls,pie_data[,2],sep=" ") # add percents to labels 
  lbls <- paste(lbls,"%",sep="") # ad % to labels 
  
  pie(slices,labels = lbls, col=rainbow(length(lbls)),
      main=lab)
}

# Plots a piechart from a dataframe already containing percentages
plot_pie_perc <- function(df, c1, c2, lab = "Pie Chart", trsh = 2) {
  pie_data <- df[,c(c1,c2)]
  names(pie_data) <- c("unique_elements", "perc_of_occ")
  pie_data <- pie_data[order(-pie_data[,2]),]
  pie_data <- rbind(pie_data[pie_data[,2] > trsh,], data.frame(unique_elements = "Other", 
                                                               perc_of_occ = sum(pie_data[pie_data[,2] <= trsh,2])))
  slices <- pie_data[,2] 
  lbls <- pie_data[,1]
  lbls <- paste(lbls,pie_data[,2],sep=" ") # add percents to labels 
  lbls <- paste(lbls,"%",sep="") # ad % to labels 
  
  pie(slices,labels = lbls, col=rainbow(length(lbls)),
      main=lab)
}

##################################
# Functions used by my functions #
##################################

# perc_of_unique() takes a vector, vec, and returns a data frame, df, where df[,1] = unique(vec) and df[,2] = the percentage of occurence of df[,1] in vec.
# @param vec: a vector of unique or non-unique elements
# @return: a dataframe where the first row contains only unique elements, and the second the percentage of occurence of that element in the vector
sum_perc_of_unique <- function(vec){
  unique_elements <- unique(vec)
  sum_of_occ <- round(sapply(unique_elements, function(x) sum(x == vec)), 1)
  perc_of_occ <- round((sum_of_occ*100/length(vec)), 1)
  df <- data.frame(unique_elements, perc_of_occ, sum_of_occ)
  df_ordered <- df[order(-df$perc_of_occ),]
  return(df_ordered)
}


# get_host_names() takes a character vector of the names of the files containing prokaryotic host genomes, 
# and reutrns a vector of the coresponding genera or family
# @param vec: character vector of file-names
# @param prok_list: dataframe with the information on the putative host. Found in the input files
# @param rank: not implemented yet. taxonomic rank of the vector to return ("g" - genus, or "f" - family)
# @param col_fname: column index of the column containing the filenames
# @param col_gen: column index of the column containing the genus names
# @return: vector of genera corresponding to the bacteria
get_host_names<- function(vec, prok_l = prok_list, col_fname = 1, col_gen = 3){
  genera <- unlist(lapply(vec, function(x) as.character(prok_l[x == prok_l[,col_fname], col_gen])))
  return(genera)
}

# sample_matrix() samples a matrix randomly
# @param sample_row
# @param sample_col
# @param seed
# @param sample_size_row 
# @param sample_size_col 
sample_matrix <- function(matrix, sample_row = TRUE, sample_col = FALSE, seed = 1, sample_size_row = 80, sample_size_col = 2000){
  mx <- matrix
  set.seed(seed)
  rowsample <- sample.int(nrow(mx), sample_size_row)
  colsample = sample.int(ncol(mx), sample_size_col)
  if (sample_row == TRUE){
    mx2 <- mx[rowsample,]
  } else {
    mx2 <- matrix
  }
  if (sample_col == TRUE){
    mx3 <- mx2[,colsample]
  } else {
    mx3 <- mx2
  }
  return(mx3)
}


# get_tax() takes a vector of string elements (that correspond to a family or higher level taxonomic rank) and returns a dataframe with the unique elements in the first column, and the corresponding family name in the second column
# @param vec: string vector of taxonomic information
# @param tax: taxonomic level
# @return:
get_tax <- function(vec, tax = c("family", "order")){
  unique_elements <- unique(vec)
  tax <- t(sapply(unique_elements, function(x) get_fop(x)))
  return(data.frame(unique_elements,tax))
}

# collapse_and_sum() takes a dataframe and 2 columns, and returns a dataframe where the first column only has unique elements, and the second columns's cells are summed up for duplicates
# @param df: dataframe of interest
# @param c1: column intdex (integer) of the column with duplicate elements in df
# @param c2: column index (integer) of numeric (to sum) elements in df
collapse_and_sum <- function(df, c1, c2){
  unique_c1 <- unique(df[,c1])
  sum_c2 <- sapply(unique_c1, function(x) sum(df[df[,c1] == x,c2]))
  df_out <- data.frame(unique_c1, sum_c2)
  return(df_out)
}


# merge_and_sum() takes 2 dataframes and 2 specified columns of each. The first specified columns are compared, and if they match, the second specified columns are summed
merge_and_sum <- function(dfa, dfb, c1a, c2a, c1b, c2b){
  merged_df <- merge(dfa[, c(c1a, c2a)], dfb[, c(c1b, c2b)], by.x = names(dfa)[c1a], by.y = names(dfb)[c1b], all = TRUE)
  merged_df[is.na(merged_df)] <- 0
  merged_df[,2] <- merged_df[,2] + merged_df[,3]
  merged_df[,3] <-  round(merged_df[,2]*100/sum(merged_df[,2]), 1)
  names(merged_df) <- c("family", "count", "perc")
  return(merged_df)
}
