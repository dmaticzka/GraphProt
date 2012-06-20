# first part: this is how RNAcontext does it

# best for my own computations:
d <- read.table("results.csv", col.names=c('pred','protein','type','set','length','apr'));

# get best parameters from training data tests
selection <- subset(ddply(subset(d, pred=='train'), .(protein,type,set), transform, max_apr=max(apr)), apr==max_apr)[,c(2,3,4,5)]

# select results from "test on test"
selected_values <- join(selection, subset(d, pred=='test'))

# compute mean and decide between full and weak
selected_mean_apr <- ddply(selected_values, .(protein,type), summarize, mean_apr=mean(apr))

# select type with maximum apr
max_apr <- subset(ddply(selected_mean_apr, .(protein), transform, max_apr=max(mean_apr)), mean_apr==max_apr)

# get parameters for best values
result <- join(max_apr,selected_values)

write.table(result[,c(1,2,5,6,8,3)], "best.csv", quote=F, sep="\t", row.names=F)


# second part: we just use the full sets
fullselection <- subset(ddply(subset(d, pred=='train' & type=='full'), .(protein,type,set), transform, max_apr=max(apr)), apr==max_apr)[,c(2,3,4,5)]
full_selected_values <- join(fullselection, subset(d, pred=='test'))
mean_selected_values <- ddply(full_selected_values, .(protein,type), transform, mean_apr=mean(apr))
write.table(mean_selected_values[,c(1,2,3,4,6,7)], "best_full_only.csv", quote=F, sep="\t", row.names=F)

# everything
best_separate_types <- join(selected_mean_apr, selected_values)
write.table(best_separate_types[,c(1,2,4,5,7,3)], "best_all.csv", quote=F, sep="\t", row.names=F)
