# take contig names and length
#print a barplot nr_of_sequences vs. sequence_length
#boxplot of the sequence length
#how many/what percentage of the sequences are below this and that? (WIsH paper)


## Actually, how long are my contigs in general?
exp_predictions_k7 <- predictions_k7
exp_predictions_k7$contig_length <- as.numeric(apply(data.frame(exp_predictions_k7$contig), 1, function(x) strsplit(as.character(x), "_", fixed= TRUE)[[1]][4]))

mean(exp_predictions_k7$contig_length) #2997
median(exp_predictions_k7$contig_length) #1065

boxplot(exp_predictions_k7$contig_length, ylim = c(0,4500)) #the IQR is around 1000-

ggplot(data = exp_predictions_k7, aes(y = contig_length, x = "")) +
  geom_boxplot(outlier.colour="black", outlier.shape='_', outlier.size=2, notch=FALSE, fill = "steelblue") +
  theme_gdocs() +
  xlab("Contig length") +
  ylab("Base number") +
  ylim(0,3000)

quantile(exp_predictions_k7$contig_length) # 500    693   1065   2086  165703

sum(exp_predictions_k7$contig_length >= 3000) # 1682 
sum(exp_predictions_k7$contig_length >= 3000)/nrow(exp_predictions_k7) # 17%

sum(exp_predictions_k7$contig_length >= 2000) #2522
sum(exp_predictions_k7$contig_length >= 2000)/nrow(exp_predictions_k7) # 26%

nrow(exp_predictions_k7) #9626 contigs


############ barplot

# frequency bar-plot
head(exp_predictions_k7)

bar_dat <- data.frame(
  contig_lengths = sort(unique(exp_predictions_k7$contig_length)))

bar_dat$freq <- sapply(bar_dat$contig_lengths, function(x) 
  sum(x == exp_predictions_k7$contig_length))

max(bar_dat$contig_length)
min(bar_dat$contig_length)
max(bar_dat$freq)

p<-ggplot(data=bar_dat, aes(x = contig_lengths, y = freq)) +
  geom_bar(stat = "identity", fill = 'steelblue', width = 1000) +
  #geom_point(colour = "steelblue")+
  labs(title = "Contig length frequency", x = "Contig length", y = "Frequency of occurence") +
  theme_gdocs() 
p

plot(x = bar_dat$contig_lengths, y = bar_dat$freq)
