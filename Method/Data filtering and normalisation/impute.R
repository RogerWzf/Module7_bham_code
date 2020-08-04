library("impute")

# Input the data
peaks_list <- read.table("Peak_intensity_matrix_norm.tsv", header=F)

# Change the column and row name
colnames(peaks_list) <- peaks_list[1,]
rownames(peaks_list) <- peaks_list[,1]

# Remove the redundant rown and column
peaks_list <- peaks_list[,-1]
peaks_list <- peaks_list[-1,]

# Impute the data using KNN
imputed <- impute.knn(as.matrix(peaks_list), k = 10, rowmax = 0.8, colmax = 0.7, maxp = 1500)

write.csv(imputed$data,file='imputed_data.csv')
