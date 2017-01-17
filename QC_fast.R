# Methylation QC pipeline (quick)


library(RColorBrewer)
library(minfi)
library(ENmix)
library(GEOquery)
library(wateRmelon) # for BMIQ normalization
library(preprocessCore)  # different normalization functions

currentDate <- Sys.Date() # to save date in name of output files

`%notin%` <- function(x,y) !(x %in% y)

# Read IDAT files
setwd("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/P160281_MethylationEPIC_NicolaBeer/03.archive/P160281_MethylationEPIC_NicolaBeer.idats")
#read cross-reactive probes
cross_reactive_EPIC <- read.csv("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/cross-reactive_probes_Pidsley_2016/EPIC_cross_reactive_probes.csv",sep=";")
rgSet <- read.metharray.exp("all_for_R") # Reads an entire methylation array experiment 

rgSet

  # assayData: 1052641 features, 32 samples 
  # array: IlluminaHumanMethylationEPIC
  # annotation: ilm10b2.hg19

head(sampleNames(rgSet)) # we see the samples are named following a standard IDAT naming convention with a 10 digit number

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


############ probes QC  ####################

# transform to methylset

mraw <- preprocessRaw(rgSet) #  A MethylSet object contains only the methylated and unmethylated signals. 


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

par(mar=c(5.1, 4.1, 1.1, 3.1), xpd=TRUE)
 
multifreqpoly(assayData(mraw)$Meth+assayData(mraw)$Unmeth,xlab="Total intensity",legend=F)

#  A RatioSet object is a class designed to store Beta values and/or M values instead of the methylated and unmethylated signals. 

RSet <- ratioConvert(mraw, what = "both", keepCN = TRUE)
RSet
  # Preprocessing
  # Method: Raw (no normalization or bg correction)
  # minfi version: 1.18.6
  # Manifest version: 0.3.0

# The functions getBeta, getM and getCN return respectively the Beta value matrix, M value matrix and the Copy Number matrix.

# beta <- getBeta(RSet)
# multifreqpoly(beta,main="Multifreqpoly",xlab="Beta value",legend=F)

# # The function mapToGenome applied to a RatioSet object will add genomic coordinates to each probe together with some additional annotation 
# # information. The output object is a GenomicRatioSet 
# 
 GRset <- mapToGenome(RSet,mergeManifest=TRUE)
 GRset

#Poor performing probes are generally filtered out prior to differential methylation analysis. 
# As the signal from these probes is unreliable, by removing them we perform fewer statistical 
# tests and thus incur a reduced multiple testing penalty. We filter out probes that have failed 
# in one or more samples based on detection p-value.

#detP <- detectionP(rgSet) # calculate detection p-vals 

# # % of probes in total and for each sample that have a detection p-value above 0.01:
# failed <- detP>0.01
# colMeans(failed)*100 # % of failed positions per sample
# sum(rowMeans(failed)>0.5) # How many positions failed in >50% of samples?
# 
# # remove poor quality samples
# 
# keep <- colMeans(detP) <0.01
# rgSet = rgSet[,keep]  
# # remove from detection p-val table:
# detP <- detP[,keep]
# dim(detP) # no samples lost 03/01/2017

# # ensure probes are in the same order in the mSetSq and detP objects 
# detP <- detP[match(featureNames(GRset),rownames(detP)),]
# # remove any probes that have failed in one or more samples 
# keep <- rowSums(detP < 0.01) == ncol(GRset) 
# table(keep)
# (table(keep)[1]/table(keep)[2])*100 # % of failed probes in one or more samples
# 
# GRsetFlt <- GRset[keep,] 
# GRsetFlt # lost 4548 probes (0.5%)

# remove probes from sex chromosomes
anno=getAnnotation(GRset)

# if your data includes males and females, remove probes on the sex chromosomes 
discard <- intersect(featureNames(GRset), anno$Name[anno$chr %in% c("chrX","chrY")])
length(discard) # 19681
(length(discard)/length(featureNames(GRset)))*100 # to lose 19681 probes (2.3%)

# removal of probes where common SNPs may affect the CpG. You can either remove all probes 
# affected by SNPs (default), or only those with minor allele frequencies greater than a 
# specified value.

# remove probes with SNPs at CpG site 
# test= getSnpInfo(GRsetFlt, snpAnno = NULL)
# 
# GRsetFlt <- dropLociWithSnps(GRsetFlt,snpAnno = NULL) # is this using the annotation from the object?
# 
# GRsetFlt # lost 28453 probes (3.4%)

# Drop cross-hybridizing probes (aprox 5%).

probeNames <- featureNames(GRset)
discard2 <- intersect(probeNames,as.character(cross_reactive_EPIC[,1])) # keep only probes that are not cross reactive
length(discard2) # 43254
(length(discard2)/length(featureNames(GRset)))*100   #5%

discard=unique(c(discard,discard2)) # 62935
(length(discard)/length(featureNames(GRset)))*100  # 7.10% of the original probes will be lost

############################################ normalization

# QN function 

##get detection p-values:

dp <- detectionP(rgSet, type = "m+u")

## Type II probes
TypeII.Name <- getProbeInfo(rgSet, type = "II")$Name
TypeII.Green <- getGreen(rgSet)[getProbeInfo(rgSet, type = "II")$Address,]
TypeII.Red <- getRed(rgSet)[getProbeInfo(rgSet, type = "II")$Address,]
rownames(TypeII.Red) <- TypeII.Name
colnames(TypeII.Red) <- sampleNames(rgSet)
rownames(TypeII.Green) <- TypeII.Name
colnames(TypeII.Green) <- sampleNames(rgSet)

## Type I probes, split into green and red channels
TypeI.Green.Name <- getProbeInfo(rgSet, type = "I-Green")$Name
TypeI.Green.M <- getGreen(rgSet)[getProbeInfo(rgSet, type = "I-Green")$AddressB,] # address for methylated signal=B
rownames(TypeI.Green.M) <- TypeI.Green.Name
colnames(TypeI.Green.M) <- sampleNames(rgSet)
TypeI.Green.U <- getGreen(rgSet)[getProbeInfo(rgSet, type = "I-Green")$AddressA,] # address for unmethylated signal=A. How does he know?
rownames(TypeI.Green.U) <- TypeI.Green.Name
colnames(TypeI.Green.U) <- sampleNames(rgSet)

TypeI.Red.Name <- getProbeInfo(rgSet, type = "I-Red")$Name
TypeI.Red.M <- getRed(rgSet)[getProbeInfo(rgSet, type = "I-Red")$AddressB,]
rownames(TypeI.Red.M) <- TypeI.Red.Name
colnames(TypeI.Red.M) <- sampleNames(rgSet)
TypeI.Red.U <- getRed(rgSet)[getProbeInfo(rgSet, type = "I-Red")$AddressA,]
rownames(TypeI.Red.U) <- TypeI.Red.Name
colnames(TypeI.Red.U) <- sampleNames(rgSet)

##remove high missingness probes
d = ifelse(dp<0.01,1,NA) # if p-value <0.1, make it 1. Else, make it NA.
cr = data.frame(rowSums(is.na(d))/length(d[1,])) # sums NAs in each row, then divide by nº of samples.
  # In the end, it excludes probes with p-values above 1 in a min of 1 sample. Why the formula, then?
exclude.badcalls = rownames(cbind(cr,rownames(cr))[cbind(cr,rownames(cr))[,1]>0.02,])

# test=names(as.matrix(cr)[which(as.matrix(cr)[,1]>0.02),])  # does the same as exclude.badcalls

exclude.sites = unique(c(exclude.badcalls,discard))  #  4548 failed probes + previous ones = 65321

##exclude.sites = unique(rbind(as.matrix(exclude.chrX), as.matrix(exclude.chrY),as.matrix(exclude.cas),as.matrix(exclude.snps),as.matrix(crossmap),as.matrix(exclude.badcalls), as.matrix(exclude.mhc)))

##remove high missingness samples 

mind = data.frame(colSums(is.na(d))/length(d[,1])) # sums NAs in each column, and divides by nº of probes
remove.mind = rownames(cbind(mind,rownames(mind))[cbind(mind,rownames(mind))[,1]>0.02,])  # remove samples with over 2% failed probes
# test2=names(as.matrix(mind)[which(as.matrix(mind)[,1]>0.02),])  # does the same as remove.mind
# no failed samples

samples = rownames(cbind(mind,rownames(mind))[cbind(mind,rownames(mind))[,1]<0.1,])  # samples with less than 10% failed probes. Why this step?
#samples=pData(rgSet)$Sample_Name

TypeII.Green =subset(TypeII.Green, select=samples)
TypeII.Red = subset(TypeII.Red, select=samples)
TypeI.Green.M = subset(TypeI.Green.M, select=samples)
TypeI.Green.U = subset(TypeI.Green.U, select=samples)
TypeI.Red.M = subset(TypeI.Red.M, select=samples)
TypeI.Red.U = subset(TypeI.Red.U, select=samples)
  # Everything stays the same.

#set NAs -- this step is including NAs that shouldn't be there?????????
# d = subset(dp, select = samples) # no failed samples
# TypeII.Green = TypeII.Green * ifelse(d[rownames(TypeII.Green),]==0,1,NA) # by multiplication, set probes and samples with high p-values to NA
# TypeII.Red = TypeII.Red * ifelse(d[rownames(TypeII.Red),]==0,1,NA)
# TypeI.Green.M = TypeI.Green.M * ifelse(d[rownames(TypeI.Green.M),]==0,1,NA)
# TypeI.Green.U = TypeI.Green.U * ifelse(d[rownames(TypeI.Green.U),]==0,1,NA)
# TypeI.Red.M = TypeI.Red.M * ifelse(d[rownames(TypeI.Red.M),]==0,1,NA)
# TypeI.Red.U = TypeI.Red.U * ifelse(d[rownames(TypeI.Red.U),]==0,1,NA)
##why this step?

#--------------------------------------------------------------------------------------------------------------------------------
##calculate betas - no QN
TypeII.betas = TypeII.Green/(TypeII.Red+TypeII.Green+100)  # beta = methylated (green in type II) / meth (green) + unmeth (red) + 100
TypeI.Green.betas = TypeI.Green.M/(TypeI.Green.M+TypeI.Green.U+100) # here (type I probe)  meth and unmeth have same color
TypeI.Red.betas = TypeI.Red.M/(TypeI.Red.M+TypeI.Red.U+100) # same
beta.noQN=rbind(TypeII.betas,TypeI.Green.betas,TypeI.Red.betas)  # unifying the 3 matrices, I get as a sum the total number of probes
# colnames(beta.noQN)<-gsub("^X","",colnames(beta.noQN)) # Once I rename the sample, this is useless

# plot beta values and see how they change from previous plots


library(geneplotter)

jpeg("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/multifreq_plots_beta_values_Matthias_pipeline_noQN.jpg",height=900,width=600)

par(mfrow=c(3,1),mar=c(5, 5, 2, 5), xpd=TRUE)

multifreqpoly(beta.noQN,main="All probes",xlab="Beta value",legend=F)
# legend("topright", legend = as.factor(pD$sample),fill=as.factor(pD$sample),inset=c(-0.07,0)) # colouring error

multifreqpoly(rbind(TypeI.Red.betas,TypeI.Green.betas),main="Multifreqpoly: Infinium I",xlab="Beta value",legend=F)

multifreqpoly(TypeII.betas,main="Multifreqpoly: Infinium II",xlab="Beta value",legend=F)

dev.off()

  # my previous results and Matthias' are exactly the same

beta.noQN=as.matrix(beta.noQN)
##save(beta.noQN, file="beta_noQN.RData")
##rm(beta.noQN, TypeII.betas,TypeI.Green.betas,TypeI.Red.betas)

#############QN (quantile normalization)

#exclude sites (filtering out failed probes)
TypeII.Green = TypeII.Green[rownames(TypeII.Green) %notin% as.matrix(exclude.sites),]
TypeII.Red = TypeII.Red[rownames(TypeII.Red) %notin% as.matrix(exclude.sites),]
TypeI.Green.M = TypeI.Green.M[rownames(TypeI.Green.M) %notin% as.matrix(exclude.sites),]
TypeI.Green.U = TypeI.Green.U[rownames(TypeI.Green.U) %notin% as.matrix(exclude.sites),]
TypeI.Red.M = TypeI.Red.M[rownames(TypeI.Red.M) %notin% as.matrix(exclude.sites),]
TypeI.Red.U = TypeI.Red.U[rownames(TypeI.Red.U) %notin% as.matrix(exclude.sites),]

# normalizing using quantiles

TypeII.Red.norm=normalize.quantiles(TypeII.Red)
TypeII.Green.norm=normalize.quantiles(TypeII.Green)
TypeI.Green.M.norm=normalize.quantiles(TypeI.Green.M)
TypeI.Green.U.norm=normalize.quantiles(TypeI.Green.U)
TypeI.Red.M.norm=normalize.quantiles(TypeI.Red.M)
TypeI.Red.U.norm=normalize.quantiles(TypeI.Red.U)
#rm(TypeII.Red,TypeII.Green,TypeI.Green.M,TypeI.Green.U,TypeI.Red.M,TypeI.Red.U)

#calculate betas
TypeII.betas = TypeII.Green.norm/(TypeII.Red.norm+TypeII.Green.norm+100)
rownames(TypeII.betas)=rownames(TypeII.Green)
colnames(TypeII.betas)=colnames(TypeII.Green)
TypeI.Green.betas = TypeI.Green.M.norm/(TypeI.Green.M.norm+TypeI.Green.U.norm+100)
rownames(TypeI.Green.betas)=rownames(TypeI.Green.M)
colnames(TypeI.Green.betas)=colnames(TypeI.Green.M)
TypeI.Red.betas = TypeI.Red.M.norm/(TypeI.Red.M.norm+TypeI.Red.U.norm+100)
rownames(TypeI.Red.betas)=rownames(TypeI.Red.M)
colnames(TypeI.Red.betas)=colnames(TypeI.Red.M)

beta=rbind(TypeII.betas,TypeI.Green.betas,TypeI.Red.betas)
# colnames(beta)<-gsub("^X","",colnames(beta))

  # has 65321 less probes than the unnormalized version

jpeg("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/multifreq_plots_beta_values_Matthias_pipeline_filtered_and_QN_without_CR_probes.jpg",height=900,width=600)

par(mfrow=c(3,1),mar=c(5, 5, 2, 5), xpd=TRUE)

multifreqpoly(beta,main="All probes",xlab="Beta value",legend=F)
# legend("topright", legend = as.factor(pD$sample),fill=as.factor(pD$sample),inset=c(-0.07,0)) # colouring error
multifreqpoly(rbind(TypeI.Red.betas,TypeI.Green.betas),main="Multifreqpoly: Infinium I",xlab="Beta value",legend=F)
multifreqpoly(TypeII.betas,main="Multifreqpoly: Infinium II",xlab="Beta value",legend=F)

dev.off()

# most changes in type II probe distribution (in abs value). Also in type I.
# before and after filtering and normalization:

jpeg("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/multifreq_plots_beta_values_Matthias_pipeline_removing_CRprobes_comparison_un-normalized.jpg",height=900,width=600)

par(mfrow=c(2,1),mar=c(5, 5, 2, 5), xpd=TRUE)

multifreqpoly(beta.noQN,main="Before filtering and QN",xlab="Beta value",legend=F)

multifreqpoly(beta,main="After",xlab="Beta value",legend=F)

dev.off()




beta <- as.matrix(beta)
rm(TypeII.Red,TypeII.Green,TypeI.Green.M,TypeI.Green.U,TypeI.Red.M,TypeI.Red.U)

##########end of QNorm function

dim(beta)

m=log2(beta/(1-beta))    #M=log2(Beta/(1-Beta))

par(mfrow=c(1,1))
multifreqpoly(m,main="M values after filtering and normalization",xlab="M value",legend=F)

Nas=which(is.na(m))  # which rows have NA?


#### plot how samples cluster before and after filtering & normalization

before.quantile <- preprocessQuantile(rgSet)

# MDS

colours37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095",
              "#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d",
              "#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977",
              "#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b",
              "#bd5975") # larger selection of colours


jpeg("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/mds_pairwise_stages_before_after_filtering_QN.jpg",height=700,width=1400)

par(mfrow=c(1,2),mar=c(5, 5, 2, 5), xpd=TRUE)

plotMDS(getM(before.quantile), top=1000, gene.selection="pairwise", col=colours37[as.factor(pD$stage)])
legend("topleft", legend=levels(as.factor(pD$stage)), text.col=colours37, bg="white", cex=0.7)

plotMDS(m, top=1000, gene.selection="pairwise", col=colours37[as.factor(pD$stage)])
legend("topright", legend=levels(as.factor(pD$stage)), text.col=colours37, bg="white", cex=0.7)
dev.off()

# PCA

plot_pca=function(x){
  pca1<-prcomp(t(x), retx=TRUE,scale. = F)   # transpose so that samples are in rows, and probes in columns
  
  # plot(pca1, type = "l") #variance vs first 10 components
  summary(pca1)          #importance of each component (important line is "proportion of variance")
  
  percentVar <- (pca1$sdev)^2 / sum(pca1$sdev^2)
  percentVar <- round(100 * percentVar)
  pcs <- as.data.frame(pca1$x)
  pcs <- cbind(pcs,sample=sample,stage=stage)
  pcs$sample=factor(pcs$sample,levels=unique(sample))
  pcs$stage=factor(pcs$stage,levels=unique(stage))
  # pcs$stage <- ordered(pcs$stage, levels = stage)
  
  p <- ggplot(pcs, aes(PC1, PC2, colour=stage, shape=sample)) + 
    geom_point(size = 3) + xlab (paste0( "PC1:" ,percentVar[ 1 ],"% variance")) + 
    ylab (paste0( "PC2: ",percentVar[ 2 ],"% variance" ))
  p <- p + theme(	panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                  panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, colour = "black"), 
                  legend.key = element_blank(),# legend.position = c(0.5,0.5),
                  axis.title.y = element_text(face="bold", angle=90, size=12, vjust=0.2),
                  axis.title.x = element_text(face="bold", size=12, vjust=0),
                  axis.text.x = element_text(face="bold", colour = "black", angle=90, size=12, vjust=0.2, hjust =1 ),
                  axis.text.y = element_text(face="bold", colour = "black"),
                  axis.ticks = element_line(colour = "black"),
                  axis.line = element_line(colour = "black"),
                  legend.position = "bottom",
                  legend.text = element_text(size=7))
  
  return(p)
}

p1=plot_pca(m)
ggsave(paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/PCA_after_filtering_and_QN",currentDate,".jpg",sep=""),p1,width=8,height=8,units="in",dpi=300)

# plot without adult islets (R and ISL)

library(pheatmap)

plot_sdc=function(x){
  sampleDists <- dist(t(x))   
  #This function computes and returns the distance matrix computed by using the specified distance measure to compute the distances between the rows of a data matrix.
  # by default, "euclidean"
  sampleDistMatrix <- as.matrix(sampleDists)
 
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  #png("/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/distance_matrix_voom_peaks.png", type="cairo",
  #    width=7,height=5,units="in",res=200,pointsize = 13)
  
  sdc= pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
  
  # dev.off()
  
  return(sdc)
}

sdc=plot_sdc(m) 

# png("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/sdc_after_filtering_QN.png",height=9,width=9,units="in",res=300)
print(sdc)
# dev.off()

