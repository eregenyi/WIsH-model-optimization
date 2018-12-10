"
This script contains functions that are essential for other scripts. 

Functions:

"


print("functions.R is now running...")


##################
# Main functions #
##################

#takes 3 tables containing the predicted and true host genus, family, order and phylum, and creates an accuracy table.
# @param p: p-value cutoff. default: 1
# @param t1_in, t2_in, t3_in: outputs of get_host_fop_pred_and_true(), equal number of row required!
get_bp_acc_table <- function(t1_in, t2_in, t3_in, p = 1){
  #subset according to the p-value
  t1 <- t1_in[t1_in[,2]<=p,]
  t2 <- t2_in[t2_in[,2]<=p,]
  t3 <- t3_in[t3_in[,2]<=p,]
  #get accuracies
  nrow_1 <- nrow(t1) 
  nrow_2 <- nrow(t2) 
  nrow_3 <- nrow(t3) 
  ac_1 <- c(sum(t1[,11])/nrow_1, sum(t1[,12])/nrow_1, sum(t1[,13])/nrow_1, sum(t1[,14])/nrow_1)
  ac_2 <- c(sum(t2[,11])/nrow_2, sum(t2[,12])/nrow_2, sum(t2[,13])/nrow_2, sum(t2[,14])/nrow_2)
  ac_3 <- c(sum(t3[,11])/nrow_3, sum(t3[,12])/nrow_3, sum(t3[,13])/nrow_3, sum(t3[,14])/nrow_3)
  
  ac_tab <- data.frame(ac_1, ac_2, ac_3)
  return(ac_tab)
}


# file_to_contig_names() is a function for finding the contig names
# correspoding to the given input file names. (which are the identifiers of sequences in WIsH's output) 
# @param file_names: a list of filenames, for which we want to find the corresponding contig names
# @param phages_list: a data-frame containing the file names in the first, and the corresponding contig names in the second column
# @return: corresponding contig names
file_to_contig_names <- function(file_names, phages_list){
  
  node_names_list <- lapply(file_names, function(x) 
    as.character(phages_list[(phages_list[,1] == x), 2]))
  node_names <- unlist(node_names_list)
  
  if (length(node_names) != 0 ){
    return(node_names)
  } else {
    stop("Returned value is of length 0, check your input.")
  }
}


# parse_into_top_10_predictions() parses the provided input into a table with 
# the contig name, the p-value, and the genus names of the 10 best predicted hosts (in separate columns)
# @param contig_list:
# @param ll:
# @param predictions:
# @param prok_list:
# @return:
parse_into_top_10_predictions <- function(contig_list, ll, predictions, prok_list){
  #top_10 with filenames
  top_10 <- t(sapply(contig_list, function(x) get_top_10(ll,x))) 
  #get genus names instead of filenames
  top_10 <- t((apply(top_10, 1, function(x) # apply to every row (1),
    (sapply(x, function(y) prok_list$genus[prok_list$file_name == y]))))) # lapply every element of that row
  #create final dataframe
  preds <- merge(x = data.frame(contig = contig_list), y = predictions, 
                 by.x = "contig", by.y = "contig")
  contig_pvalue_top_10 <- data.frame(preds[,c("contig", "p-value")], top_10)
  return(contig_pvalue_top_10)
}

# get_host_fop_pred_and_true() given a 
# @param
# @param
# @param
# @return
get_host_fop_pred_and_true <- function(predictions, col_nr_pred, col_nr_true){
  #taxonomic level of info to retrieve
  tax_levels <- c("family", "order", "phylum")
  #retrieve fop for predicted hosts
  fop_pred <- tax_name(query = predictions[,col_nr_pred], get = tax_levels, db = "ncbi", rows = 2)[,3:5]
  names(fop_pred) <- c("pred_host_fam", "pred_host_ord",  "pred_host_phy")
  #retrieve fop for true hosts
  fop_true <- tax_name(query = predictions[,col_nr_true], get = tax_levels, db = "ncbi", rows = 2)[,3:5]
  names(fop_true) <- c("true_host_fam", "true_host_ord",  "true_host_phy")
  #create output
  predictions_fop <- data.frame(predictions, fop_pred, fop_true)
  
  #get indices for is_correct
  nr_col <- ncol(predictions_fop)
  cols_fop_pred <- c(col_nr_pred, (nr_col-5):(nr_col-3))
  cols_fop_true <- c(col_nr_true, (nr_col-2):nr_col)
  
  #get is_correct
  predictions_fop_is_correct <- get_fop_is_correct(predictions_fop, cols_fop_pred, cols_fop_true)
  return(predictions_fop_is_correct)
}

# darken a colour
# @param color:
# @param factor:
darken <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  return(col)
}


# lighten a colour
# @param color:
# @param factor:
lighten <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(col), maxColorValue=255)
  return(col)
}


######################
# plotting functions #
######################

# plot_al_gfop() plots the accuracy (of the prediction) against the length of the sequence, for genus, family, order and phylum level
# @param predictions_fop_is_correct:
# @param col_contig_name:
# @param cols_is_correct:
# @param color:
# @return: void (no return value)
plot_al_gfop <- function(predictions_fop_is_correct, col_contig_name = 1, cols_is_correct = 11:14, color ="black"){
  
  # establish variables
  pfic <- predictions_fop_is_correct
  cols <- cols_is_correct
  col_cname <- col_contig_name
  pfic$length <- unlist(lapply(pfic[,col_cname], function(x) get_contig_length(x))) #contig length
  
  #order data by contig length
  pfic_o <- pfic[order(pfic$length),]
  pfic_o$running_length <- c(rollapply(pfic_o$length, width = 10, by = 1, FUN = mean, align = "left"),
                             rep(33000,9))
  
  #add the running means
  pfic_o$running_ic_gen <- running_mean(pfic_o[,cols[1]])
  pfic_o$running_ic_fam <- running_mean(pfic_o[,cols[2]])
  pfic_o$running_ic_ord <- running_mean(pfic_o[,cols[3]])
  pfic_o$running_ic_phy <- running_mean(pfic_o[,cols[4]])
  
  # plot 
  ggplot(data = pfic_o, aes(x = running_length)) +
    #genus
    geom_smooth(aes(y = running_ic_gen), 
              se = FALSE, span = 0.14, method = "loess", 
              colour = color, linetype = "solid") + 
    #family
    geom_smooth(aes(y = running_ic_fam), 
                se = FALSE, span = 0.14, method = "loess", 
                colour = color, linetype = "longdash") + 
    #order
    geom_smooth(aes(y = running_ic_ord), 
                se = FALSE, span = 0.14,method = "loess",
                colour = color, linetype = "twodash") + 
    #phylum
    geom_smooth(aes(y = running_ic_phy), 
                se = FALSE, span = 0.14, method = "loess",
                colour = color, linetype = "dotted") + 
#try plotting with the lines bellow uncommented. That will reveal if the smoothing splines distort the trend
    #geom_line(aes(y = running_ic_gen)) + #genus
    #geom_line(aes(y = running_ic_fam)) + #family
    #geom_line(aes(y = running_ic_ord)) + #order
    #geom_line(aes(y = running_ic_phy)) + #phylum
    
    # add a theme, legend and axis labels
    theme_gdocs() +
    xlab("contig length") +
    ylab("prediction accuracy") +
    ylim(0,1)
}

# plot_ap_gfop() plots the average accuracy (of the prediction) against the p-value cutoff, for genus, family, order and phylum level
# @param predictions_fop_is_correct:
# @param col_pvalue: numeric, column number of the column containing the p-values
# @param cols_is_correct:
# @param color: color of the line graphs
# @return: void (no return value)
plot_ap_gfop <- function(predictions_fop_is_correct, col_pvalue = 2, cols_is_correct = 11:14, color ="black") {
  
  #establish variables
  pfic <- predictions_fop_is_correct
  cols <- cols_is_correct
  p <- col_pvalue
  
  # establish data to plot from
  cutoffs <- rev(seq(0.01,0.3,0.01))
  is_cor_gen <- sapply(cutoffs, function(x) mean(pfic[pfic[,p] < x, cols[1]]))
  is_cor_fam <- sapply(cutoffs, function(x) mean(pfic[pfic[,p] < x, cols[2]]))
  is_cor_ord <- sapply(cutoffs, function(x) mean(pfic[pfic[,p] < x, cols[3]]))
  is_cor_phy <- sapply(cutoffs, function(x) mean(pfic[pfic[,p] < x, cols[4]]))
  p_acc <- data.frame(cutoff = cutoffs, is_cor_gen, is_cor_fam, is_cor_ord, is_cor_phy)
  
  # plot
  ggplot(data = p_acc, aes(x = cutoff)) + #base
    #genus
    geom_smooth(aes(y = is_cor_gen), 
                se = FALSE, span = 0.4, method = "loess", 
                colour = color, linetype = "solid") +
    geom_smooth(aes(y = is_cor_fam), 
                se = FALSE, span = 0.4, method = "loess", 
                colour = color, linetype = "longdash") +
    geom_smooth(aes(y = is_cor_ord), 
                se = FALSE, span = 0.4, method = "loess", 
                colour = color, linetype = "twodash") +
    geom_smooth(aes(y = is_cor_phy), 
                se = FALSE, span = 0.4, method = "loess", 
                colour = color, linetype = "dotted") +
    #geom_line(aes(y = is_cor_gen))+ #genus
    #geom_line(aes(y = is_cor_fam))+ #family
    #geom_line(aes(y = is_cor_ord))+ #order
    #geom_line(aes(y = is_cor_phy)) + #phylum
    theme_gdocs() +
    xlab("p-value cutoff") +
    ylab("prediction accuracy") +
    xlim(0.3, 0.01) +
    ylim(0,1)
}

#plot_nobs_vs_cutoff() plot p-value cutoffs vs #observarions
# @param predictions_fop_is_correct:
# @param col_pvalue: numeric, column number of the column containing the p-value
# @param color: color of the line graphs
# @return: void (no return value)
plot_nobs_vs_cutoff <- function(predictions_fop_is_correct, col_pvalue = 2, color ="black") {
  
  #establish variables
  pfic <- predictions_fop_is_correct
  p <- col_pvalue
  
  # establish data to plot from
  cutoffs <- rev(seq(0.01,0.3,0.01))
  nr_obs <- sapply(cutoffs, function(x) sum(pfic[,p] < x))
  cutoff_nr_obs <- data.frame(cutoff = cutoffs, nobs = nr_obs)

  # plot
  ggplot(data = cutoff_nr_obs, aes(x = cutoff)) + #base
    geom_smooth(aes(y = nobs), 
                se = FALSE, span = 0.4, method = "loess", 
                colour = color, linetype = "solid") +
    #geom_line(aes(y = nobs)) +
    xlim(0.3, 0.01) +
    theme_gdocs() +
    xlab("p-value cutoff") +
    ylab("number of observations remaining in the test set")
}

##################################
# Functions used by my functions #
##################################

# get_top_10() is a function returning the best 10 predictions for a given contg name
# @param file_names: a list of filenames, for which we want to find the corresponding contig names
# @param phages_list: a data-frame containing the file names in the first, and the corresponding contig names in the second column
# @return: corresponding contig names
get_top_10 <- function(ll, contig){
  
  likelihoods <- ll[,colnames(ll) == contig]
  best_likelihoods <- head(sort(likelihoods, decreasing = TRUE), 10)
  hosts <- lapply(best_likelihoods, function(x) row.names(ll[likelihoods == x,]))
  hosts <- head(unlist(hosts), 10) #sometimes here are multiple predictions with the same likelihood, in that case we take only the top 10
  names(hosts) <- c("pred_host_1", "pred_host_2", "pred_host_3", 
                    "pred_host_4", "pred_host_5", "pred_host_6", "pred_host_7", 
                    "pred_host_8", "pred_host_9", "pred_host_10")
  return(hosts)
}

# get_fop() retrieves the family, the order and he phylum of a query sequence. (fop)
# ! rows = 2, bacuse of Bacillus! This may result in unexpected results in other cases!
# @param query: the tax name to be queried (in this code it is used n genus level queries)
# @return: qfop is a named list with the query, family, order and phylum.
get_fop <- function(query){
  tax_query <- query
  tax_levels <- c("family", "order", "phylum")
  qfop <- tryCatch(unlist(tax_name(query = tax_query, get = tax_levels, db = "ncbi", rows = 2, message = FALSE)[,2:5]), error=function(e) c(name = NA, rank = NA, id = NA))
  return(qfop)
}

#get_fop_is_correct
# @param
# @param
# @param
# @return
get_fop_is_correct <- function(predictions_fop, cols_fop_pred, cols_fop_true){
  #create copy
  predictions_fop_is_correct <- predictions_fop
  
  #compare genus
  pred_gen <- predictions_fop_is_correct[,cols_fop_pred[1]]
  true_gen <- predictions_fop_is_correct[,cols_fop_true[1]]
  predictions_fop_is_correct$is_corr_gen <- ifelse(pred_gen == true_gen, 1, 0)
  
  #compare family
  pred_fam <- predictions_fop_is_correct[,cols_fop_pred[2]]
  true_fam <- predictions_fop_is_correct[,cols_fop_true[2]]
  predictions_fop_is_correct$is_corr_fam <- ifelse(pred_fam == true_fam, 1, 0)
  
  #compare order
  pred_ord <- predictions_fop_is_correct[,cols_fop_pred[3]]
  true_ord <- predictions_fop_is_correct[,cols_fop_true[3]]
  predictions_fop_is_correct$is_corr_ord <- ifelse(pred_ord == true_ord, 1, 0)
  
  #compare phylum
  pred_phy <- predictions_fop_is_correct[,cols_fop_pred[4]]
  true_phy <- predictions_fop_is_correct[,cols_fop_true[4]]
  predictions_fop_is_correct$is_corr_phy <- ifelse(pred_phy == true_phy, 1, 0)
  
  #return
  return(predictions_fop_is_correct)
}


# get_mode() gets the most common element of an array. If multiple elements
# would be chosen, returns the one that occurs first (has the lowest index)
# @param vec:
# @return: most common element, in case of multiple possibilities, the one occuring first
get_mode <- function(vec){
  uniq_vec <- unique(vec)
  mode <- uniq_vec[which.max(tabulate(match(vec, uniq_vec)))]
  return(mode)
}

# get_contig_length extracts the length of the contig give the name. The format should be standard.
# ! contig length from the contig name. If the contig names follow a different standard, it may not work!
# @param contig_name: a string, where distinct elements are separated by "_" (underscore), 
# and the 4th such segment contains the length.
# @return: numeric value of the contig length
get_contig_length <- function(contig_name){ 
  
  length_str <- strsplit(as.character(contig_name), "_", fixed= TRUE)[[1]][4]
  length_nrc <- as.numeric(length_str)
  
  return(length_nrc)
}

# rolling_mean computes the rolling mean of a vector. The last nr-1 elements are going to be ones.
# @param vec:
# @return:
running_mean <- function(vec, windowsize = 10){
  r_mean <- c(rollapply(vec, width = windowsize, by = 1, FUN = mean, align = "left"),
       rep(1,windowsize-1))
  return(r_mean)
}

#########################
# Source more functions #
#########################

source(paste(location, "subscripts/functions/functions_lca.R", sep = ""))
source(paste(location, "subscripts/functions/functions_final_analysis.R", sep = ""))

#################
# End of script #
#################

print("functions.R is done.")