# The code was based on the WGCNA official tutorial. 
# cf.https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html
# Load the library
library(WGCNA)
library(formattable)

# Allow multiple threading
allowWGCNAThreads(nThreads = 2)

# Set the working directory.
path <- "/Users/zhifan/Downloads/tmo/processed_data/WGCNA/con_car"
setwd(path)

#Read in the peak table
peaks_list <- read.csv("PQN_KNN_Glog_peaks.csv", header = TRUE, row.names=1, quote="\"", 
                    dec = '.', fill = TRUE, comment.char="",na.strings = "NA", check.names=FALSE)

# Cluster the samples to check for outliers
sampleTree <- hclust(dist(peaks_list), method = "average") 

# Plot the sample tree
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 50, col = "blue")
# Determine cluster under the line
clust <- cutreeStatic(sampleTree, cutHeight = 55, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples <- (clust==1)
datExpr <- peaks_list[keepSamples, ]

# Loading the clinical trait data
traitData <- readxl::read_excel("sampleMetadata.xlsx")
traitData <- as.data.frame(traitData)

# remove columns that hold information we do not need.
names(traitData)
allTraits <- traitData[, -c(2, 3, 6, 8)]

# Form a data frame analogous to expression data that will hold the clinical traits.
selected_Samples <- rownames(datExpr)
traitRows <- match(selected_Samples, allTraits$sample)
datTraits <- allTraits[traitRows, -1]
row.names(datTraits) <- allTraits[traitRows, 1]
collectGarbage()

# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector <- powers, verbose = 5, indent = 1)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 <-- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Automatic, one-step network construction and module detection
net <- blockwiseModules(datExpr, power = 14, networkType = "unsigned",
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "DaphniaTOM", 
                       verbose = 3)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors <- labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Relating modules to external clinical traits
# 3.a Quantifying module–trait associations
# Define numbers of genes and samples
nSamples <- nrow(datExpr)
# Recalculate MEs with color labels
moduleColors <- labels2colors(net$colors)
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigenmetabolites
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# 3.b Gene relationship to trait and important modules: Gene Significance and Module Membership

# Define variable Treatment containing the Treatment column of datTrait
Treatment <- as.data.frame(datTraits$Treatment)
names(Treatment) <- "Treatment"

# names (colors) of the modules
modNames <- substring(names(MEs), 3)

# Calculate the module membership
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep="")
names(MMPvalue) <- paste("p.MM", modNames, sep="")

# Calculate the gene significant
geneTraitSignificance <- as.data.frame(cor(datExpr, Treatment, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(Treatment), sep="")
names(GSPvalue) <- paste("p.GS.", names(Treatment), sep="")

# 3.c Intramodular analysis: identifying genes with high GS and MM
module <- "brown"
column <- match(module, modNames)
moduleGenes <- moduleColors==module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Peak significance for Carbaryl treatment",
                   main = paste("Module membership vs. peaks significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# Summary the output
names(datExpr)[moduleColors=="brown"]
annot <- read.csv(file = "../KEGG_beamspy_annot.csv")
annot$mz <- digits(annot$mz,7)
#annot <- annot[rownames(datExpr),]
probes <- names(datExpr)
probes2annot <- match(probes, annot$mz)
# The following is the number or probes without annotation
sum(is.na(probes2annot))
# Should return 0.
# Create the starting data frame
geneInfo0 <- data.frame(mz = probes,
                       intensity = annot$intensity[probes2annot],
                       molecular_formula = annot$molecular_formula[probes2annot],
                       compound_name = annot$compound_name[probes2annot],
                       compound_id = annot$compound_id[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder <- order(-abs(cor(MEs, Treatment, use = "p")))
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)){
        oldNames = names(geneInfo0)
        geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                               MMPvalue[, modOrder[mod]]);
        names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                             paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Treatment))
geneInfo <- geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo.csv")

# 5. Visualization of networks within R

# 5.a Visualizing the gene network

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM <- 1-TOMsimilarityFromExpr(datExpr, power = 12)
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM <- dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) <- NA
# Call the plot function
sizeGrWindow(9,9)
geneTree <- net$dendrograms[[1]]
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all metabolites")