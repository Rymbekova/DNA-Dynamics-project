rm(list=ls())
setwd("~/Downloads/DRD/Input_data")
suppressMessages(library(minfi)) 

# Step 8. Perform a Principal Component Analysis on the beta matrix obtained previously
load("/home/aigerim/Downloads/DRD/Input_data/RGSET.RData")
preprocessQuantile_results <- preprocessQuantile(RGSET)
beta_preprocessQuantile <- getBeta(preprocessQuantile_results)
pheno <- read.csv("~/Downloads/DRD/Input_data/Samplesheet_report_2020.csv",header=T, stringsAsFactors=T)
pca_results <- prcomp(t(beta_preprocessQuantile),scale=T)
print(summary(pca_results))
plot(pca_results)
palette(c("orange","purple"))
plot(pca_results$x[,1], pca_results$x[,2],cex=2,pch=2,col=pheno$Group,xlab="PC1",ylab="PC2",xlim=c(-700,500),ylim=c(-700,700))
text(pca_results$x[,1], pca_results$x[,2],labels=rownames(pca_results$x),cex=0.5,pos=1)
legend("bottomright",legend=levels(pheno$Group),col=c(1:nlevels(pheno$Group)),pch=2)
# The variability seems to be quite low and at the plot there is no clear pattern in forming clusters, probably the dimensionality should be higher than 2


# Step 9. Identify differentially methylated probes between DS and WT samples groups using t-test
t_test <- t.test(beta_preprocessQuantile[1,] ~ pheno$Group)
My_ttest_function <- function(x) {
  t_test <- t.test(x~ pheno$Group)
  return(t_test$p.value)
} 
pValues_ttest_entire <- apply(beta_preprocessQuantile,1, My_ttest_function)
final_ttest_entire <- data.frame(beta_preprocessQuantile, pValues_ttest_entire)
head(final_ttest_entire)
final_ttest_entire <- final_ttest_entire[order(final_ttest_entire$pValues_ttest_entire),]
head(final_ttest_entire)

# Step 10. Apply multiple test correction and set a threshold of 0.05. How many probes are differentially methylated considering nominal pValues? How many after Bonferroni correction, BH correction?
raw_pValues <- final_ttest_entire[,9]
corrected_pValues_BH <- p.adjust(raw_pValues,"BH")
corrected_pValues_Bonf <- p.adjust(raw_pValues,"bonferroni")
final_ttest_entire_corrected <- data.frame(final_ttest_entire, corrected_pValues_BH, corrected_pValues_Bonf)
head(final_ttest_entire_corrected)
dim(final_ttest_entire[final_ttest_entire$pValues_ttest_entire<=0.05,]) # how many probes identified as differentially methylated considering nominal pValues
dim(final_ttest_entire[final_ttest_entire$corrected_pValues_BH<=0.05,]) # how many probes identified as differentially methylated after Bonferroni correction
dim(final_ttest_entire[final_ttest_entire$corrected_pValues_Bonf<=0.05,]) # how many probes identified as differentially methylated after BH correction

# Step 11. A heatmap of the top 100 differentially mehtylated probes 
library(gplots)
input_heatmap=as.matrix(final_ttest_entire_corrected[1:100,1:8])
pheno$Group
colorbar <- c("green","green","orange","orange","green","green","orange","orange")
heatmap.2(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F)
heatmap.2(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'single'),dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F)
heatmap.2(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'average'),dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F)
# Hierchical clustering using different methods gives different results

# Step 12. A volcano plot and a Manhattan plot of the results of differential methylation analysis 
# For a volcano plot we create beta values of DS and WT samples, then we calculate the mean for each row of each group 
beta <- final_ttest_entire_corrected[,1:8]
beta_groupDS <- beta[,pheno$Group=="DS"]
mean_beta_groupDS <- apply(beta_groupDS,1,mean)
beta_groupWT <- beta[,pheno$Group=="WT"]
mean_beta_groupWT <- apply(beta_groupWT,1,mean)
delta <- mean_beta_groupWT-mean_beta_groupDS
head(delta)

toVolcPlot <- data.frame(delta, -log10(final_ttest_entire_corrected$pValues_ttest_entire))
head(toVolcPlot)
plot(toVolcPlot[,1], toVolcPlot[,2])
plot(toVolcPlot[,1], toVolcPlot[,2],pch=16,cex=0.5)
# We can add a threshold for p value significance
abline(a=-log10(0.01),b=0,col="red")
# We can highlight with yellow the probes with p value <0.01 and delta >0.01
toHighlight <- toVolcPlot[toVolcPlot[,1]>0.1 & toVolcPlot[,2]>(-log10(0.01)),]
plot(toVolcPlot[,1], toVolcPlot[,2],pch=16,cex=0.5)
points(toHighlight[,1], toHighlight[,2],pch=16,cex=0.7,col="yellow")
# We can highlight with red the probes with delta >0.01
toHighlight <- toVolcPlot[abs(toVolcPlot[,1])>0.1 & toVolcPlot[,2]>(-log10(0.01)),]
plot(toVolcPlot[,1], toVolcPlot[,2],pch=16,cex=0.5)
points(toHighlight[,1], toHighlight[,2],pch=16,cex=0.7,col="red")

# To create a Manhattan plot we need the gap package and Illumina manifest for genome annotation information for each probe
library(gap)
load("/home/aigerim/Downloads/DRD/Lesson 2-20200514/Illumina450Manifest_clean.RData")
final_ttest_entire_corrected <- data.frame(rownames(final_ttest_entire_corrected),final_ttest_entire_corrected)
# To merge two dataframes we need to use a column that has the same name in both 
colnames(final_ttest_entire_corrected)[1] <- "IlmnID"
colnames(final_ttest_entire_corrected)
final_ttest_entire_corrected_annotated <- merge(final_ttest_entire_corrected, Illumina450Manifest_clean,by="IlmnID")
# Create an input for the Manhattan plot + reorder the levels of the chromosome column
input_Manhattan <- data.frame(final_ttest_entire_corrected_annotated$CHR, final_ttest_entire_corrected_annotated$MAPINFO, final_ttest_entire_corrected_annotated$pValues_ttest_entire)
str(input_Manhattan)
levels(input_Manhattan$final_ttest_entire_corrected_annotated.CHR)
input_Manhattan$final_ttest_entire_corrected_annotated.CHR <- factor(input_Manhattan$final_ttest_entire_corrected_annotated.CHR,levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"))
levels(input_Manhattan$final_ttest_entire_corrected_annotated.CHR)
palette <- rainbow(24)
mhtplot(input_Manhattan,control=mht.control(colors=palette))
axis(2,cex=0.5)
abline(a=-log10(0.01),b=0)
