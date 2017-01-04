# Methylation QC pipeline (quick)


library(RColorBrewer)
library(minfi)
library(ENmix)
library(GEOquery)

# Read IDAT files
#Currently minfi does not support reading compressed IDAT files. 
setwd("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/P160281_MethylationEPIC_NicolaBeer/03.archive/P160281_MethylationEPIC_NicolaBeer.idats")


rgSet <- read.metharray.exp("all_for_R") # Reads an entire methylation array experiment 

rgSet

  # assayData: 1052641 features, 32 samples 
  # array: IlluminaHumanMethylationEPIC
  # annotation: ilm10b2.hg19


head(sampleNames(rgSet))


# we see the samples are named following a standard IDAT naming convention with a 10 digit 
# number which is an array identifier followed by an identifier of the form R01C01.

# The 200526570053_R01C01 means row 1 and column 1 on chip 200526570053. 


# We now get the standard GEO representation to get the phenotype data stored in GEO. 
# Most of the columns in this phenotype data are irrelevant (contains data such as the 
# address of the person who submitted the data); we keep the useful ones.

# Then we clean it.

pheno=read.csv2(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/P160281_MethylationEPIC_NicolaBeer/03.archive/3.Results/Phenotype.csv")

pD <- pheno[,c(1,10)] # we keep sample name in the experiment and the corresponding stages and samples


# add stage column

stage= c("iPSC","iPSC","iPSC","DE","DE","DE","PGT","PGT","PGT","PFG","PFG","PFG","PE","PE","PE","EP","EN6","EN6","EN6","EN7","EN7","124I","141 (14-19)","177","179","182","184","124R","136","137","138","139")
# add sample column
sample=c("neo1.1","SBAd2.1","SBAd3.1","neo1.1","SBAd2.1","SBAd3.1","neo1.1","SBAd2.1","SBAd3.1","neo1.1","SBAd2.1","SBAd3.1","neo1.1","SBAd2.1","SBAd3.1","SBAd3.1","neo1.1","SBAd2.1","SBAd3.1","neo1.1","SBAd3.1","ISL","ISL","ISL","ISL","ISL","ISL","R","R","R","R","R")

  # some samples are missing

pD=cbind(pD,stage,sample)


pD$Sample_Name=gsub( "^.*?_", "", pD$Sample_Name) # take out stage part of experiment id
rownames(pD) <- pD$Sample_Name # add that as rownames

rgSet <- rgSet[,pD$Sample_Name] # reordering before merging



pData(rgSet) <- pD  # merging into pheno data of rgSet

rgSet

# give more descriptive names to our samples:

pD$simple_id=paste(pD$stage,pD$sample,sep="-")
sampleNames(rgSet) <- pD$simple_id



#########trying ENmix tutorial############


mraw <- preprocessRaw(rgSet) #  A MethylSet object contains only the methylated and unmethylated signals. 

#This function matches up the different probes and color channels. Note that the dimension of this object is much smaller than for 
# the RGChannelSet; this is because CpGs measured by type I probes are measured by 2 probes.

# The accessors getMeth and getUnmeth can be used to get the methylated and unmethylated intensities matrices:


head(getMeth(mraw)[,1:3])

head(getUnmeth(mraw)[,1:3])

##########modifying multifreqpoly###########

bincount <- function(x,breaks){
  x <- x[!is.na(x)]
  bc <- table(.bincode(x, breaks, TRUE, TRUE))
  temp=data.frame(id=c(1: (length(breaks)-1)))
  bc=data.frame(id=as.numeric(names(bc)),counts=as.numeric(bc))
  resu=merge(temp,bc,by="id",all.x=TRUE)
  resu$counts[is.na(resu$counts)]=0
  resu$counts[order(resu$id)]
}

multifreqpoly <- function(mat, nbreaks=100, col=1:ncol(mat), xlab="", 
                          ylab="Frequency", legend = list(x = "topright", fill=col,inset=c(-0.2,0),cex=0.6,
                                                          legend = if(is.null(colnames(mat))) paste(1:ncol(mat)) else 
                                                            colnames(mat)),...)
{
  if(!is.matrix(mat)) stop("Warning: input data is not a numeric matrix\n")
  if(is.null(col)) col="black"
  col=rep(col,ceiling(ncol(mat)/length(col)))
  if(nbreaks > nrow(mat)) nbreaks=min(15,round(nrow(mat)/2))
  
  breaks <- seq(min(mat,na.rm=TRUE), max(mat,na.rm=TRUE), 
                diff(range(mat,na.rm=TRUE))/nbreaks)
  mids <- 0.5 * (breaks[-1] + breaks[-length(breaks)])
  counts <- sapply(data.frame(mat),bincount,breaks=breaks)
  plot(range(mids),c(0,max(counts)),type="n",xlab=xlab,ylab=ylab,...)
  for(i in 1:ncol(counts)){lines(mids,counts[,i],col=col[i],...)}
  if(is.list(legend)) do.call(graphics::legend, legend)
}
#####################

#total intensity plot is userful for data quality inspection
#and identification of outlier samples

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
 
multifreqpoly(assayData(mraw)$Meth+assayData(mraw)$Unmeth,xlab="Total intensity")
 
 #frequency plots of beta values

 beta<-getBeta(mraw, "Illumina")
 anno=getAnnotation(rgSet)   # get the EPIC Annotation data

 
 beta1=beta[anno$Type=="I",]
 beta2=beta[anno$Type=="II",]


library(geneplotter)
 
jpeg("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/multifreq_plots_beta_values_methylset.jpg",height=900,width=600)

par(mfrow=c(3,1),mar=c(5, 5, 2, 5), xpd=TRUE)

multifreqpoly(beta,main="Multifreqpoly",xlab="Beta value")

multifreqpoly(beta1,main="Multifreqpoly: Infinium I",xlab="Beta value")

multifreqpoly(beta2,main="Multifreqpoly: Infinium II",xlab="Beta value")
 
dev.off()



# When type I and type II probes are plotted separately the difference in modes between type I and II probes
# can be appreciated. But when all probes are plotted together (Fig 4 top panels), the multidensity
# plot obscures these differences, while they remain readily apparent in the multifreqpoly plot. In
# addition, the multidensity plots appear to suggest that probes range in value from <0 to >1, whereas
# multifreqpoly correctly show the range from 0 to 1
#

###############

#  Data quality measures, including detection P values, number of beads for each methylation read
#  and average intensities for bisulfite conversion probes can be extracted using the function QCinfo
#  from an object of RGChannelSetExtended. According default or user specified quality score
#  thresholds, the QCinfo can also identify and export a list of low quality samples and CpG probes.
#  Outlier samples in total intensity or beta value distribution were often excluded before further
#  analysis. Such samples were tricky to be identified, by default the argument outlier=TRUE will
#  trigger the function to identify these outlier samples automatically. Quality score figures from
#  QCinfo can be used to guide the selection of quality score thresholds. Low quality samples and
#  probes can be filtered out using QCfilter or preprocessENmix.



# throws ERROR: rgSet is not an object of class 'RGChannelSetExtended'
#  qc<-QCinfo(rgSet)
#  #exclude before backgroud correction
# mdat<-preprocessENmix(rgSet, bgParaEst="oob", dyeCorr=TRUE, QCinfo=qc, nCores=6)
# #Or exclude after background correction
# mdat <- QCfilter(mdat,qcinfo=qc,outlier=TRUE)








# detection p-vals Maksimovic paper

detP <- detectionP(rgSet)
head(detP)
dim(detP)

# % of probes in total and for each sample that have a detection p-value above 0.01:
failed <- detP>0.01
colMeans(failed)*100 # % of failed positions per sample
sum(rowMeans(failed)>0.5) # How many positions failed in >50% of samples?


pal <- brewer.pal(5,"Dark2")

jpeg("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/mean_detection_pvals.jpg",height=700,width=700)
par(mfrow=c(1,1),oma=c(4, 2, 3, 3),xpd=TRUE)
barplot(colMeans(detP),col=pal[as.factor(pD$sample)],las=2, cex.names=0.9,ylim = c(0,0.0003))
mtext(side = 2, text = "Mean detection p-values", line = 4.5, cex=1.5)
abline(h=0.01,col="red") 
legend("topright", legend=levels(as.factor(pD$sample)),fill=pal, bg="white",ncol=1,inset=c(-0.05,0))
dev.off()

# remove poor quality samples

keep <- colMeans(detP) <0.01
rgSet = rgSet[,keep]
rgSet  # no samples lost 03/01/2017

# remove from detection p-val table:
detP <- detP[,keep]
dim(detP) # no samples lost 03/01/2017




########################################## OK up to here, check QC below

# following minfi manual

# The RGChannelSet stores also a manifest object that contains the probe design information of the array:

manifest <- getManifest(rgSet)
manifest

head(getProbeInfo(manifest))


#

#  A RatioSet object is a class designed to store Beta values and/or M values instead of the methylated and unmethylated signals. 
# An optional copy number matrix, CN, the sum of the methylated and unmethylated signals, can be also stored. Mapping a MethylSet 
# to a RatioSet may be irreversible, i.e. one cannot be guaranteed to retrieve the methylated and unmethylated signals from a RatioSet. 
# A RatioSet can be created with the function ratioConvert:

RSet <- ratioConvert(mraw, what = "both", keepCN = TRUE)
RSet

  # Preprocessing
  # Method: Raw (no normalization or bg correction)
  # minfi version: 1.18.6
  # Manifest version: 0.3.0


# The functions getBeta, getM and getCN return respectively the Beta value matrix, M value matrix and the Copy Number matrix.

beta <- getBeta(RSet)



anno=getAnnotation(RSet)   # get the EPIC Annotation data


beta1=beta[anno$Type=="I",]
beta2=beta[anno$Type=="II",]



jpeg("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/multifreq_plots_beta_values_ratioset.jpg",height=900,width=600)

par(mfrow=c(3,1),mar=c(5, 5, 2, 5), xpd=TRUE)

multifreqpoly(beta,main="Multifreqpoly",xlab="Beta value")

multifreqpoly(beta1,main="Multifreqpoly: Infinium I",xlab="Beta value")

multifreqpoly(beta2,main="Multifreqpoly: Infinium II",xlab="Beta value")

dev.off()

# we see that the frequencies of the beta values haven't changed 

# same with M values:

mval<-getM(RSet)


mval1=mval[anno$Type=="I",]
mval2=mval[anno$Type=="II",]

# I can't plot directly, because there are NaN and Inf values (that are tried to be used by the function to set the min x value)


mval[is.infinite(mval)] <- NA
mval <- na.omit(mval)

mval1[is.infinite(mval1)] <- NA
mval1 <- na.omit(mval1)

mval2[is.infinite(mval2)] <- NA
mval2 <- na.omit(mval2)

jpeg("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/multifreq_plots_M_values.jpg",height=900,width=600)

par(mfrow=c(3,1),mar=c(5, 5, 2, 5), xpd=TRUE)

multifreqpoly(mval,main="Multifreqpoly",xlab="M value")

multifreqpoly(mval1,main="Multifreqpoly: Infinium I",xlab="M value")

multifreqpoly(mval2,main="Multifreqpoly: Infinium II",xlab="M value")

dev.off()

# remember that all that has been done doesn't include any normalization


####### FIX THE FOLLOWING WITH CORRECT ANNOTATION ###########
# The function mapToGenome applied to a RatioSet object will add genomic coordinates to each probe together with some additional annotation 
# information. The output object is a GenomicRatioSet (class holding M or/and Beta values together with associated genomic coordinates). 
# It is possible to merge the manifest object with the genomic locations by setting the option mergeManifest to TRUE.

GRset <- mapToGenome(RSet,mergeManifest=TRUE)
GRset

# It is possible to merge the manifest object with the genomic locations by
# setting the option mergeManifest to TRUE.


# Note that the GenomicRatioSet extends the class SummarizedExperiment. Here are the main accessors functions to access the data:

beta <- getBeta(GRset)
M <- getM(GRset)
CN <- getCN(GRset)

sampleNames <- sampleNames(GRset)
probeNames <- featureNames(GRset)
pheno <- pData(GRset)


# To return the probe locations as a GenomicRanges objects, one can use the accessor granges:

gr <- granges(GRset)
head(gr, n= 3)

# To access the full annotation, one can use the command getAnnotation:


annotation <- getAnnotation(GRset)
names(annotation)




############### QC ##################################



# minfi provides a simple quality control plot that uses the log median intensity in both the methylated (M) and unmethylated (U) channels. 
# When plotting these two medians against each other, it has been observed that good samples cluster together, while failed samples tend to 
# separate and have lower median intensities. In order to obtain the methylated and unmethylated signals, we need to convert the RGChannelSet 
# to an object containing the methylated and unmethylated signals using the function preprocessRaw. It takes as input a RGChannelSet and converts 
# the red and green intensities to methylated and unmethylated signals according to the special 450K probe design, and returns the converted signals
# in a new object of class MethylSet. It does not perform any normalization.

# The accessors getMeth and getUnmeth can be used to get the methylated and unmethylated intensities matrices:

head(getMeth(mraw)[,1:3])
head(getUnmeth(mraw)[,1:3])

# The functions getQC and plotQC are designed to extract and plot the quality control information from the MethylSet:

qc <- getQC(mraw)
head(qc)

par(mfrow=c(1,1)) 
plotQC(qc)

# Moreover, the function addQC applied to the MethylSet will add the QC information to the phenotype data.

# To further explore the quality of the samples, it is useful to look at the Beta value densities of the samples, with the option to color 
# the densities by sampleID. Already done previously with frequencies.

# densityPlot(mraw, sampGroups = pD$sample)

# or density bean plots:

densityBeanPlot(mraw, sampGroups = pD$sample)





########### # Control probes plot

# The 450k array contains several internal control probes that can be used to assess the quality control of different sample preparation 
# steps (bisulfite conversion, hybridization, etc.). The values of these control probes are stored in the initial RGChannelSet and can be 
# plotted by using the function controlStripPlot and by specifying the control probe type:

controlStripPlot(rgSet, controls="BISULFITE CONVERSION II")


# qcReport(rgSet, pdf= "qcReport.pdf")   # it's not plotting properly


############################################ normalization

# The function preprocessQuantile implements stratified quantile normalization preprocessing for Illumina methylation microarrays. 
# Probes are stratified by region (CpG island, shore, etc.).


# normalize the data using Quantiles; this results in a GenomicRatioSet object 
gset.quantile <- preprocessQuantile(rgSet)

# another function, more apt for data with expected large scale differences:

gset.funnorm <- preprocessFunnorm(rgSet)

# visualise what the data looks like before and after normalisation 
jpeg("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/before_after_quantile_normalization.jpg",height=700,width=1400)
par(mfrow=c(2,1),mar=c(5, 5, 2, 5), xpd=TRUE)
densityPlot(rgSet, sampGroups=pD$sample,main="Raw", legend=FALSE) 
legend("topright", legend = levels(as.factor(pD$sample)), text.col=brewer.pal(8,"Dark2"),inset=c(-0.07,0))

densityPlot(getBeta(gset.quantile), sampGroups=pD$sample, main="Normalized", legend=FALSE)
legend("topright", legend = levels(as.factor(pD$sample)), text.col=brewer.pal(8,"Dark2"),inset=c(-0.07,0))

dev.off()

jpeg("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/before_after_funnorm_normalization.jpg",height=700,width=1400)
par(mfrow=c(2,1),mar=c(5, 5, 2, 5), xpd=TRUE)
densityPlot(rgSet, sampGroups=pD$sample,main="Raw", legend=FALSE) 
legend("topright", legend = levels(as.factor(pD$sample)), text.col=brewer.pal(8,"Dark2"),inset=c(-0.07,0))

densityPlot(getBeta(gset.funnorm), sampGroups=pD$sample, main="Normalized", legend=FALSE)
legend("topright", legend = levels(as.factor(pD$sample)), text.col=brewer.pal(8,"Dark2"),inset=c(-0.07,0))

dev.off()


# MDS plots
library(limma)
colours37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095",
             "#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d",
             "#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977",
             "#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b",
             "#bd5975") # larger selection of colours

pal <- brewer.pal(5,"Dark2")

jpeg("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/mds_common_samples_stages_quantilenorm.jpg",height=700,width=1400)

par(mfrow=c(1,2),mar=c(5, 5, 2, 5), xpd=TRUE)

plotMDS(getM(gset.quantile), top=1000, gene.selection="common", col=colours37[as.factor(pD$stage)])
legend("topright", legend=levels(as.factor(pD$stage)), text.col=colours37, bg="white", cex=0.7)


plotMDS(getM(gset.quantile), top=1000, gene.selection="common", col=pal[as.factor(pD$sample)])
legend("topright", legend=levels(as.factor(pD$sample)), text.col=pal, bg="white", cex=0.7)

dev.off()

jpeg("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/mds_pairwise_samples_stages_quantilenorm.jpg",height=700,width=1400)

par(mfrow=c(1,2),mar=c(5, 5, 2, 5), xpd=TRUE)

plotMDS(getM(gset.quantile), top=1000, gene.selection="pairwise", col=colours37[as.factor(pD$stage)])
legend("topright", legend=levels(as.factor(pD$stage)), text.col=colours37, bg="white", cex=0.7)


plotMDS(getM(gset.quantile), top=1000, gene.selection="pairwise", col=pal[as.factor(pD$sample)])
legend("topright", legend=levels(as.factor(pD$sample)), text.col=pal, bg="white", cex=0.7)
dev.off()


# funnorm

jpeg("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/mds_common_samples_stages_funnorm.jpg",height=700,width=1400)
par(mfrow=c(1,2),mar=c(5, 5, 2, 5), xpd=TRUE)

plotMDS(getM(gset.funnorm), top=1000, gene.selection="common", col=colours37[as.factor(pD$stage)])
legend("topright", legend=levels(as.factor(pD$stage)), text.col=colours37, bg="white", cex=0.7)

plotMDS(getM(gset.funnorm), top=1000, gene.selection="common", col=pal[as.factor(pD$sample)])
legend("topright", legend=levels(as.factor(pD$sample)), text.col=pal, bg="white", cex=0.7)
dev.off()

jpeg("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/mds_pairwise_samples_stages_funnorm.jpg",height=700,width=1400)
par(mfrow=c(1,2),mar=c(5, 5, 2, 5), xpd=TRUE)

plotMDS(getM(gset.funnorm), top=1000, gene.selection="pairwise", col=colours37[as.factor(pD$stage)])
legend("topright", legend=levels(as.factor(pD$stage)), text.col=colours37, bg="white", cex=0.7)

plotMDS(getM(gset.funnorm), top=1000, gene.selection="pairwise", col=pal[as.factor(pD$sample)])
legend("topright", legend=levels(as.factor(pD$sample)), text.col=pal, bg="white", cex=0.7)
dev.off()

# Examine higher dimensions to look at other sources of variation 

jpeg("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/mds_common_samples_stages_secondaryPCs.jpg",height=700,width=2100)
par(mfrow=c(1,3)) 

plotMDS(getM(mSetSq), top=1000, gene.selection="common", col=colours37[as.factor(pD$stage)], dim=c(1,3))
legend("topright", legend=levels(as.factor(pD$stage)), text.col=colours37, cex=0.7, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common", col=colours37[as.factor(pD$stage)], dim=c(2,3))
legend("topright", legend=levels(as.factor(pD$stage)), text.col=colours37, cex=0.7, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common", col=colours37[as.factor(pD$stage)], dim=c(3,4))
legend("topright", legend=levels(as.factor(pD$stage)), text.col=colours37, cex=0.7, bg="white")

dev.off()

# # Violin Plots of methylation levels across the whole genome
# library(vioplot)
# vioplots_samples <- data.frame()
# 
# for(s in colnames(rgSet)){
#   
#   vioplots_samples[,s] = rgSet[,s]
#   
#   
# }
# vioplot(rgSet, names=pD$simple_id, col="gold")
# title("Violin Plots of methylation levels across the whole genome")


############ Filtering

#Poor performing probes are generally filtered out prior to differential methylation analysis. 
# As the signal from these probes is unreliable, by removing them we perform fewer statistical 
# tests and thus incur a reduced multiple testing penalty. We filter out probes that have failed 
# in one or more samples based on detection p-value.

# ensure probes are in the same order in the mSetSq and detP objects 
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]
# remove any probes that have failed in one or more samples 
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)
mean(!keep)*100 # % of failed probes in one or more samples

mSetSqFlt <- mSetSq[keep,] 
mSetSqFlt # lost 4548 probes (0.5%)

# remove probes from sex chromosomes

# if your data includes males and females, remove probes on the sex chromosomes 
keep <- !(featureNames(mSetSqFlt) %in% anno$Name[anno$chr %in% c("chrX","chrY")])
table(keep) 
mSetSqFlt <- mSetSqFlt[keep,] # lost 19077 probes (2.2%)

# removal of probes where common SNPs may affect the CpG. You can either remove all probes 
# affected by SNPs (default), or only those with minor allele frequencies greater than a 
# specified value.

# remove probes with SNPs at CpG site 
test= getSnpInfo(mSetSqFlt, snpAnno = NULL)

mSetSqFlt <- dropLociWithSnps(mSetSqFlt,snpAnno = NULL) # is this using the annotation from the object?

mSetSqFlt # lost 28453 probes (3.4%)

######test
source("/Users/Marta/Documents/WTCHG/R\ scripts/methylation\ from\ matthias/qnorm_function.R")
library(preprocessCore)
`%notin%` <- function(x,y) !(x %in% y)

##get detection p-values:

dp <<- detectionP(RGset, type = "m+u")





###### our pipeline – QN in 6 categories

## Type II probes
TypeII.Name <- getProbeInfo(RGset, type = "II")$Name
TypeII.Green <- getGreen(RGset)[getProbeInfo(RGset, type = "II")$Address,]
TypeII.Red <- getRed(RGset)[getProbeInfo(RGset, type = "II")$Address,]
rownames(TypeII.Red) <- TypeII.Name
colnames(TypeII.Red) <- sampleNames(RGset)
rownames(TypeII.Green) <- TypeII.Name
colnames(TypeII.Green) <- sampleNames(RGset)



## Type I probes, split into green and red channels
TypeI.Green.Name <- getProbeInfo(RGset, type = "I-Green")$Name
TypeI.Green.M <- getGreen(RGset)[getProbeInfo(RGset, type = "I-Green")$AddressB,]
rownames(TypeI.Green.M) <- TypeI.Green.Name
colnames(TypeI.Green.M) <- sampleNames(RGset)
TypeI.Green.U <- getGreen(RGset)[getProbeInfo(RGset, type = "I-Green")$AddressA,]
rownames(TypeI.Green.U) <- TypeI.Green.Name
colnames(TypeI.Green.U) <- sampleNames(RGset)

TypeI.Red.Name <- getProbeInfo(RGset, type = "I-Red")$Name
TypeI.Red.M <- getRed(RGset)[getProbeInfo(RGset, type = "I-Red")$AddressB,]
rownames(TypeI.Red.M) <- TypeI.Red.Name
colnames(TypeI.Red.M) <- sampleNames(RGset)
TypeI.Red.U <- getRed(RGset)[getProbeInfo(RGset, type = "I-Red")$AddressA,]
rownames(TypeI.Red.U) <- TypeI.Red.Name
colnames(TypeI.Red.U) <- sampleNames(RGset)

##remove high missingness samples
d = ifelse(dp<0.1,1,NA)
cr = data.frame(rowSums(is.na(d))/length(d[1,]))
exclude.badcalls = rownames(cbind(cr,rownames(cr))[cbind(cr,rownames(cr))[,1]>0.02,])

exclude.sites = exclude.badcalls
##exclude.sites = unique(rbind(as.matrix(exclude.chrX), as.matrix(exclude.chrY),as.matrix(exclude.cas),as.matrix(exclude.snps),as.matrix(crossmap),as.matrix(exclude.badcalls), as.matrix(exclude.mhc)))

mind = data.frame(colSums(is.na(d))/length(d[,1]))
remove.mind = rownames(cbind(mind,rownames(mind))[cbind(mind,rownames(mind))[,1]>0.02,])
samples = rownames(cbind(mind,rownames(mind))[cbind(mind,rownames(mind))[,1]<0.1,])
#samples=pData(RGset)$Sample_Name

TypeII.Green =subset(TypeII.Green, select=samples)
TypeII.Red = subset(TypeII.Red, select=samples)
TypeI.Green.M = subset(TypeI.Green.M, select=samples)
TypeI.Green.U = subset(TypeI.Green.U, select=samples)
TypeI.Red.M = subset(TypeI.Red.M, select=samples)
TypeI.Red.U = subset(TypeI.Red.U, select=samples)

##set NAs
d = subset(dp, select = samples)
TypeII.Green = TypeII.Green * ifelse(d[rownames(TypeII.Green),]==0,1,NA)
TypeII.Red = TypeII.Red * ifelse(d[rownames(TypeII.Red),]==0,1,NA)
TypeI.Green.M = TypeI.Green.M * ifelse(d[rownames(TypeI.Green.M),]==0,1,NA)
TypeI.Green.U = TypeI.Green.U * ifelse(d[rownames(TypeI.Green.U),]==0,1,NA)
TypeI.Red.M = TypeI.Red.M * ifelse(d[rownames(TypeI.Red.M),]==0,1,NA)
TypeI.Red.U = TypeI.Red.U * ifelse(d[rownames(TypeI.Red.U),]==0,1,NA)

#--------------------------------------------------------------------------------------------------------------------------------
##calculate betas - no QN
TypeII.betas = TypeII.Green/(TypeII.Red+TypeII.Green+100)
TypeI.Green.betas = TypeI.Green.M/(TypeI.Green.M+TypeI.Green.U+100)
TypeI.Red.betas = TypeI.Red.M/(TypeI.Red.M+TypeI.Red.U+100)
beta.noQN=rbind(TypeII.betas,TypeI.Green.betas,TypeI.Red.betas)
colnames(beta.noQN)<-gsub("^X","",colnames(beta.noQN))
beta.noQN=as.matrix(beta.noQN)
##save(beta.noQN, file="beta_noQN.RData")
##rm(beta.noQN, TypeII.betas,TypeI.Green.betas,TypeI.Red.betas)

##QN

#exclude sites
TypeII.Green = TypeII.Green[rownames(TypeII.Green) %notin% as.matrix(exclude.sites),]
TypeII.Red = TypeII.Red[rownames(TypeII.Red) %notin% as.matrix(exclude.sites),]
TypeI.Green.M = TypeI.Green.M[rownames(TypeI.Green.M) %notin% as.matrix(exclude.sites),]
TypeI.Green.U = TypeI.Green.U[rownames(TypeI.Green.U) %notin% as.matrix(exclude.sites),]
TypeI.Red.M = TypeI.Red.M[rownames(TypeI.Red.M) %notin% as.matrix(exclude.sites),]
TypeI.Red.U = TypeI.Red.U[rownames(TypeI.Red.U) %notin% as.matrix(exclude.sites),]


TypeII.Red.norm=normalize.quantiles(TypeII.Red)
TypeII.Green.norm=normalize.quantiles(TypeII.Green)
TypeI.Green.M.norm=normalize.quantiles(TypeI.Green.M)
TypeI.Green.U.norm=normalize.quantiles(TypeI.Green.U)
TypeI.Red.M.norm=normalize.quantiles(TypeI.Red.M)
TypeI.Red.U.norm=normalize.quantiles(TypeI.Red.U)
#rm(TypeII.Red,TypeII.Green,TypeI.Green.M,TypeI.Green.U,TypeI.Red.M,TypeI.Red.U)

#calculate betas
TypeII.betas = TypeII.Green.norm/(TypeII.Red.norm+TypeII.Green.norm+100)
TypeI.Green.betas = TypeI.Green.M.norm/(TypeI.Green.M.norm+TypeI.Green.U.norm+100)
TypeI.Red.betas = TypeI.Red.M.norm/(TypeI.Red.M.norm+TypeI.Red.U.norm+100)
beta=rbind(TypeII.betas,TypeI.Green.betas,TypeI.Red.betas)
colnames(beta)<-gsub("^X","",colnames(beta))
row.names(beta)=c(row.names(TypeII.Green), row.names(TypeI.Green.M), row.names(TypeI.Red.M))#
beta <- as.matrix(beta)





##########end of test


