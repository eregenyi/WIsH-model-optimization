# merge_into_lca() gets the by row of for the given columns, 
# and replaces them with 3 columns by default: name, rank and id
# @param df: dataframe containing taxonomic information, including the taxa for which lca is to be retrieved. (lca is retrieved by row.)
# @param col_taxa: column indices containing the taxa to merge
# @param name: logical, if TRUE, the lca name is returned
# @param rank: logical, it TRUE, the rank of the lca is returned
# @param id: logical, if TRUE, the id of the lca is returned
# @return: information on the lca (name rank id) as specified by parameters name rank and id.
merge_into_lca <- function(df, col_taxa, name = TRUE, rank = TRUE, id = TRUE){
  #do some error handling:
  to_return <- c(name, rank, id)
  is_gt_one(vec = to_return, message = "At least one of the parameters \"name\", \"rank\" or \"id\" has to be TRUE.")
  #establish input, process data and return
  df_tax <- df[,col_taxa]
  print(df_tax)
  df_lca <- t((apply(df_tax, 1, function(x) get_lca(x))))
  colnames(df_lca) <- c("lca_name", "lca_rank", "lca_id")
  return(df_lca[,to_return])
}

# is_ancestor() is a function that takes a dataframe with 2 speciied columns (ancestor and target),
# compares this ancestor and target by row, and returns a binary vecor, 
# with 1 for cases where ancestor is indeed the ancestor of the target, 0 otherwise.
# @param df:
# @param: col_ancestor
# @param: col_target
is_ancestor <- function(df, col_ancestor = 3, col_target = 5){
  #establish input
  ancestors <- df[,col_ancestor]
  targets <- df[,col_target]
  ta <- data.frame(ancestors, targets)
  #compute and return
  is_correct <- as.numeric(t((apply(ta, 1, function(x) get_lca(x))))[,1] == ancestors)
  return(is_correct)
}

# get_gen_rate() returns the rate of "genus" given a vecotr
# @param vec: 
# @return:
get_gen_rate <- function(vec){
  rate <- round(sum(vec == "genus")/length(vec), 2)
  return(rate)
}

######################
# Plotting functions #
######################

# plot_acc_vs_nr_lca() plots the accuracy versus the 
# number of best predictions of which the lca was taken for the final prediction
# @param df_acc_bp: a dataframe with the first column contaning the number of best predictions, the second the accuracy estimate
# @param color: color of the line
# @return: ggplot object (an be put into variable, othervise it just plots)
plot_acc_vs_nr_bp <- function(df_acc_bp, color = "black"){
  df <- df_acc_bp
  names(df)[1:2] <- c("bp", "acc")
  
  ggplot(data = df, aes(x = bp)) +
    geom_line(aes(y = acc), color = color) +
    theme_gdocs() +
    xlab("Number of best predictions used") +
    ylab("Prediction accuracy") +
    ylim(0,1) +
    scale_x_discrete(limits=1:10)
}

plot_acc_and_genperc_vs_nr_bp <- function(df_acc_bp, color_acc = "black", color_perc = "darkgrey"){
  df <- df_acc_bp
  names(df)[1:3] <- c("bp", "acc", "genperc")
  
  ggplot(data = df, aes(x = bp)) +
    geom_line(aes(y = acc), color = color_acc, size = 1.3) +
    geom_line(aes(y = genperc), color = color_perc, size = 1.3, linetype="dashed") +
    theme_gdocs() +
    xlab("Number of best predictions used") +
    ylab("Prediction accuracy & percentage of predictions being at genus level") +
    ylim(0,1) +
    scale_x_discrete(limits=1:10)
}

# plot_bar_tax_ranks() creates a barplot depicting the percentage of predictions of each taxonomic rank
# @param df_with_rank: a dataframe of the predictions containing the ranks in one column
# @param col_rank: the indedx of the column containing the ranks
# @return: ggplot
# could improve: implement colorchoice
plot_bar_tax_ranks <- function(df_with_rank, col_rank = 4, colour = "steelblue"){
  df <- df_with_rank
  #possible ranks
  df_rank_perc <- data.frame(rank = c("genus", "family", "order", 
                  "class", "phylum", "below-superkingdom",  "superkingdom", 
                  "below-no rank"))
  #get percentage of predictions of each rank
  for (i in 1:nrow(df_rank_perc)){
    df_rank_perc$perc[i] <- round(sum(as.character(df_rank_perc[i,1]) ==  df[,col_rank])*100/nrow(df), 1)
  }
  #plot
  df_rank_perc$rank <- c("genus", "family", "order", 
                         "class", "phylum", "bellow s.kingdom", 
                         "s. kingdom", "above s.kingdom")
  df_rank_perc$rank <- factor(df_rank_perc$rank, levels = df_rank_perc$rank)
  ggplot(data=df_rank_perc, aes(x = rank, y = perc)) +
    geom_bar(stat = "identity", color = "black", fill = colour) +
    geom_text(aes(label = perc), vjust = -1.6, color = colour, size = 3.5) +
    ylab(label = "Percentage of predictions [%]") +
    xlab(label = "Taxonomic rank") +
    theme_gdocs() +
    ylim(0,60) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}


#plot_overall_accuracy creates a scatterplot 
plot_overall_accuracy <- function(df_plot, col1 = "purple", col2 = "steelblue", col3 = "orange"){
  df <- df_plot
  names(df) <- c("strategy","order6", "order7", "order8")
  df$strategy <- factor(df$strategy, levels = df$strategy)
  
  ggplot(data = df, aes(x = strategy)) +
    geom_point(aes(y = order6), shape = "___", size = 6, color = col1) +
    geom_point(aes(y = order7), shape = "___", size = 6, color = col2) +
    geom_point(aes(y = order8), shape = "___", size = 6, color = col3) +
    ylab(label = "Accuracy") +
    ylim(0,1) +
    xlab(label = "Accary boosting strategy") +
    ggtitle("Oveall accuracy vs boosting strategy") +
    #guides(colour = c(col1,col2,col3), shape = "___", label = c("7", "6", "8"))+
    theme_gdocs() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
} 

################################
#function used by my functions #
################################


# get_lca() returns the lowest common ancestor of a given list of taxa
# @param taxa_list:
# @return: the lowest common ancestor name, rank and id, as present in ncbi's database
get_lca <- function (taxa_list) {
  return(tryCatch(unlist(lowest_common(taxa_list, db = "ncbi", rows = 2)), error=function(e) c(name = NA, rank = NA, id = NA)))
}


# is_gt_one() throws an error if the sum of the provided vector, and prints the message
# @param vec:
# @param message:
# @return: void
is_gt_one <- function(vec, message = "Oops, something went wrong. Check your input parameters (esp the logical ones!)!"){
  if (sum(vec) < 1) {
    stop(message)
  }
}
