# Load the packages
if(!require("glmpca", quietly = TRUE)){
  install.packages("glmpca",dependencies=TRUE)
  library("glmpca")
}

if(!require("ggplot2", quietly = TRUE)){
  install.packages("ggplot2",dependencies=TRUE)
  library("ggplot2")
}

# Set the working directory.
path <- "/Users/zhifan/Downloads/tmo/processed_data/Treatment/car_con"
setwd(path)

# Load the metadata
sample_sheet <- readxl::read_excel("sampleMetadata.xlsx")
sample_sheet <- as.data.frame(sample_sheet)

# load peaks
peaks_matrix_file <- 'PQN_KNN_Glog_peaks.csv'
peaks_matrix <- read.csv(peaks_matrix_file, header = T, row.names = 1)
peaks_matrix <- t(peaks_matrix)

# Run the glmPCA
gpca <- glmpca(peaks_matrix, L=2, fam="nb")
gpca.dat <- gpca$factors
gpca.dat$Genotype <- sample_sheet$Genotype
gpca.dat$time <- sample_sheet$Time_point
gpca.dat$treatment <- sample_sheet$Treatment

# Plot
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = time, shape = Genotype)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA") +
  theme_bw()+
  geom_text(aes(label=rownames(gpca.dat)),hjust=0, vjust=0)
