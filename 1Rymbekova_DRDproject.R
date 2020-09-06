rm(list=ls())
setwd("~/Downloads/DRD/Input_data")
suppressMessages(library(minfi)) 

# Step 1. Import a raw data and create RGset containing the RGChannelSet
SampleSheet <- read.table("~/Downloads/DRD/Input_data/Samplesheet_report_2020.csv",sep=",",header=T)
SampleSheet
# Set the directory with the raw data
baseDir <- ("~/Downloads/DRD/Input_data")
targets <- read.metharray.sheet(baseDir)
targets
# Create RGSet with read.metharray.exp function
RGSET <- read.metharray.exp(targets = targets)
save(RGSET,file="RGSET.RData")

# Step 2. Create Red and Green fluorescences dataframes with functions getRed and getGreen
RED <- data.frame(getRed(RGSET))
head(RED)
GREEN <- data.frame(getGreen(RGSET))
head(GREEN)

# Step 3. What are the Red and Green fluorescences for the address 52682510 + check if the address corresponds to a Type I or a Type II probe
RED[rownames(RED)=="52682510",]
GREEN[rownames(GREEN)=="52682510",]
load("/home/aigerim/Downloads/DRD/Lesson 2-20200514/Illumina450Manifest_clean.RData")
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID=="52682510",]
# AddressA_ID with code 52682510 associated to the probe cg01882666, type II 

# Step 4. Create the object Mset.raw to extract methylated and unmethylated signals with getMeth and getUnmeth functions
MSET.RAW <- preprocessRaw(RGSET)
MSET.RAW
save(MSET.RAW,file="MSET_RAW.RData")
METH <- getMeth(MSET.RAW)
str(METH)
UNMETH <- getUnmeth(MSET.RAW)
str(UNMETH)

# Step 5. Perform the quality check of the samples by visualizing the median intesities of Meth and Unmeth channels with getQC function 
# QC plot
QC <- getQC(MSET.RAW)
QC
plotQC(QC)
# All samples have high median intensity signals - good quality
# Check the intensity of negative controls by comparing with the expected intensities for each type of control probes according to the Illumina guide
getProbeInfo(RGSET, type = "Control")
df_TypeControl <- data.frame(getProbeInfo(RGSET, type = "Control"))
str(df_TypeControl)
table(df_TypeControl$Type)
# Plot the intensity values of each type of controls probes in the samples
controlStripPlot(RGSET, controls="NEGATIVE")
# Negative controls values are below 1000 (log2(1000)=10) - all good
# Calculate detection p-values considering the given threshold 0.01
detP <- detectionP(RGSET)
str(detP)
head(detP)
failed <- detP>0.01
head(failed)
table(failed)
summary(failed)
# in total, 2177 probes have a detection p-value higher than the threshold

# Step 6. Calculate raw beta and M values and plot the densities of mean methylation values in DS and WT samples
beta <- getBeta(MSET.RAW)
head(beta)
summary(beta)
M <- getM(MSET.RAW)
head(M)
summary(M)
# Subset both matrices into DS and WT samples
betaDS <- beta[,-c(3,4,5,8)]
head(betaDS)
betaWT <- beta[,-c(1,2,6,7)]
head(betaWT)
MDS <- M[,-c(3,4,5,8)]
head(MDS)
MWT <- M[,-c(1,2,6,7)]
head(MWT)
# Apply mean() function
mean_of_betaDS <- apply(betaDS,1,mean,na.rm=T)
d_mean_of_betaDS <- density(mean_of_betaDS)
mean_of_betaWT <- apply(betaWT,1,mean,na.rm=T)
d_mean_of_betaWT <- density(mean_of_betaWT)
mean_of_MDS <- apply(MDS,1,mean,na.rm=T)
d_mean_of_MDS <- density(mean_of_MDS)
mean_of_MWT <- apply(MWT,1,mean,na.rm=T)
d_mean_of_MWT <- density(mean_of_MWT)
# Plot densities for DS samples
par(mfrow=c(1,2))
plot(d_mean_of_betaDS,main="Density of DS Beta Values",col="orange")
plot(d_mean_of_MDS,main="Density of DS M Values",col="purple")
# Plot densities for WT samples
par(mfrow=c(1,2))
plot(d_mean_of_betaWT,main="Density of WT Beta Values",col="orange")
plot(d_mean_of_MWT,main="Density of WT M Values",col="purple")

# Step 7. Normalize the data using preprocessQuantile, compare raw and normalized data + 6 panels: density plots of beta mean values, density plots of beta standard deviation values, boxplot of beta values for type I and type II probes
beta <- getBeta(MSET.RAW)
M <- getM(MSET.RAW)
load("/home/aigerim/Downloads/DRD/Lesson 2-20200514/Illumina450Manifest_clean.RData")
dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="I",]
dfI <- droplevels(dfI)
dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="II",]
dfII <- droplevels(dfII)

beta_I <- beta[rownames(beta) %in% dfI$IlmnID,]
beta_II <- beta[rownames(beta) %in% dfII$IlmnID,]

mean_of_beta_I <- apply(beta_I,1,mean)
mean_of_beta_II <- apply(beta_II,1,mean)
d_mean_of_beta_I <- density(mean_of_beta_I,na.rm=T)
d_mean_of_beta_II <- density(mean_of_beta_II,na.rm=T)

sd_of_beta_I <- apply(beta_I,1,sd,na.rm=T)
sd_of_beta_II <- apply(beta_II,1,sd,na.rm=T)
d_sd_of_beta_I <- density(sd_of_beta_I,)
d_sd_of_beta_II <- density(sd_of_beta_II)

load("/home/aigerim/Downloads/DRD/Input_data/RGSET.RData")
preprocessQuantile_results <- preprocessQuantile(RGSET)
beta_preprocessQuantile <- getBeta(preprocessQuantile_results)
beta_preprocessQuantile_I <- beta_preprocessQuantile[rownames(beta_preprocessQuantile) %in% dfI$IlmnID,]
beta_preprocessQuantile_II <- beta_preprocessQuantile[rownames(beta_preprocessQuantile) %in% dfII$IlmnID,]
mean_of_beta_preprocessQuantile_I <- apply(beta_preprocessQuantile_I,1,mean)
mean_of_beta_preprocessQuantile_II <- apply(beta_preprocessQuantile_II,1,mean)
d_mean_of_beta_preprocessQuantile_I <- density(mean_of_beta_preprocessQuantile_I,na.rm=T)
d_mean_of_beta_preprocessQuantile_II <- density(mean_of_beta_preprocessQuantile_II,na.rm=T)
sd_of_beta_preprocessQuantile_I <- apply(beta_preprocessQuantile_I,1,sd)
sd_of_beta_preprocessQuantile_II <- apply(beta_preprocessQuantile_II,1,sd)
d_sd_of_beta_preprocessQuantile_I <- density(sd_of_beta_preprocessQuantile_I,na.rm=T)
d_sd_of_beta_preprocessQuantile_II <- density(sd_of_beta_preprocessQuantile_II,na.rm=T)

par(mfrow=c(2,3))
plot(d_mean_of_beta_I,col="blue",main="Raw beta",xlim=c(0,1),ylim=c(0,6))
lines(d_mean_of_beta_II,col="red")
plot(d_sd_of_beta_I,col="blue",main="Raw SD",xlim=c(0,0.6),ylim=c(0,90))
lines(d_sd_of_beta_II,col="red")
boxplot(beta,ylim=c(0,1))
plot(d_mean_of_beta_preprocessQuantile_I,col="blue",main="preprocessQuantile beta",xlim=c(0,1),ylim=c(0,6))
lines(d_mean_of_beta_preprocessQuantile_II,col="red")
plot(d_sd_of_beta_preprocessQuantile_I,col="blue",main="preprocessQuantile SD",xlim=c(0,0.6),ylim=c(0,90))
lines(d_sd_of_beta_preprocessQuantile_II,col="red")
boxplot(beta_preprocessQuantile,ylim=c(0,1))
# With the same scales on x and y axes the difference between the plots of raw and normalized data is more clear
