# QC plots for probes
# solving problems in plotCtrl (ENmix library) 



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


######### fixing function


plotCtrl <- function(rgSet,IDorder=NULL)
{
  if(!is(rgSet, "RGChannelSet")){stop("object needs to be of class 
                                      'RGChannelSet'")}
  if(!is.null(IDorder))
  {
    if(sum(!(IDorder %in% colnames(rgSet)))>0)
    {
      idmissing=IDorder[!(IDorder %in% colnames(rgSet))]
      cat("The IDs were not found in the input data: ",idmissing,"\n")
      stop("Wrong ids in IDorder, please check")
    }else{
      rgSet=rgSet[,IDorder]
    }
  }
  

  
  ctrls<-getProbeInfo(rgSet,type="Control") # control probe info, 635 control probes
  ctrls <- ctrls[ctrls$Address %in% featureNames(rgSet),] # control probe info
  ctrl_r <- getRed(rgSet)[ctrls$Address,] # control probe intensities (red)
  ctrl_g <- getGreen(rgSet)[ctrls$Address,] # control probe intensities (green)
  contlid <- c("STAINING","EXTENSION","HYBRIDIZATION","TARGET REMOVAL",
               "BISULFITE CONVERSION I", "BISULFITE CONVERSION II","SPECIFICITY I",
               "SPECIFICITY II","NON-POLYMORPHIC","NEGATIVE",
               "NORM_A","NORM_C","NORM_G","NORM_T","NORM_ACGT")   # types to be plotted
  col <- as.vector(ctrls$Color) # color of each probe
  col[col == "-99"] <- NA   # fixing each color pallete
  col[col == "Aqua"] <- "aquamarine2"
  col[col == "Crimson"] <- "firebrick2"
  col[col == "Fuchsia"] <- "deeppink1"
  col[col == "Indigo"] <- "darkviolet"
  col[col == "Lime"] <- "yellowgreen"
  col[col == "Olive"] <- "darkolivegreen"
  col[col == "Silver"] <- "azure4"
  col[col == "Teal"] <- "cyan4"
  ctrls$Color <- col
  
  #fine above
  
  # fix missing values in ctrls
  
 
  
  for(ctype in contlid)  # loop through every type of plot
  {
    fn <- ctype;fn <- gsub(" ","_",fn)
    cat("Plotting ",fn,".jpg","\n")  # print what we are doing
    jpeg(paste(fn,".jpg",sep=""),width=1100,height=500,quality=100)
    par(mfrow=c(1,2))
    if(ctype == "NORM_ACGT"){
      cc <- ctrls[ctrls$Type %in% c("NORM_A","NORM_C","NORM_G","NORM_T"),]}
    else{cc <- ctrls[ctrls$Type %in% ctype,]}
    red <- ctrl_r[cc$Address,]
    grn <- ctrl_g[cc$Address,]
    ymax <- max(red,grn)*1.01
    if(ctype == "NEGATIVE")
    {
      par(mar=c(5, 4, 4, 2))
      colnames(red) <- 1:ncol(red)
      colnames(grn) <- 1:ncol(grn)
      boxplot(grn,ylim=c(0,ymax),main=paste(ctype," Green",sep=""),bty="o",
              xlab="Sample",ylab="Intensity",cex.lab=1.2)
      boxplot(red,ylim=c(0,ymax),main=paste(ctype," Red",sep=""),bty="o",
              xlab="Sample",ylab="Intensity",cex.lab=1.2)
    }
    else if(ctype %in% c("NORM_A","NORM_C","NORM_G","NORM_T"))
    {
      par(mar=c(5, 4, 4, 1))
      colnames(red)<-1:ncol(red)
      colnames(grn)<-1:ncol(grn)
      boxplot(grn,ylim=c(0,ymax),col=unique(as.vector(cc$Color)),main=
                paste(ctype," Green",sep=""),bty="o",xlab="Sample",ylab="Intensity",
              cex.lab=1.2)
      boxplot(red,ylim=c(0,ymax),col=unique(as.vector(cc$Color)),main=
                paste(ctype," Red",sep=""),bty="o",xlab="Sample",ylab="Intensity",
              cex.lab=1.2)
    }
    else if(ctype == "NORM_ACGT")
    {
      par(mar=c(5, 4, 4, 8.2))
      label <- c("NORM_A","NORM_C","NORM_G","NORM_T")
      colcode <- c("Red","Green","Blue","Purple")
      idx <- t(replicate(nrow(red),1:ncol(red)));
      loc <- cc$Color;loc[loc == "Red"] <- -0.2
      loc[loc == "Green"] <- -0.1
      loc[loc == "Blue"] <- 0.1
      loc[loc == "Purple"]=0.2;loc=as.numeric(loc)
      loc1 <- matrix(rep(loc,ncol(red)),ncol=ncol(red));idx=idx+loc1
      plot(idx,grn,col=as.vector(cc$Color),ylim=c(0,ymax),main=paste(ctype,
                                                                     " Green",sep=""),bty="o",xlab="Sample",ylab="Intensity",cex.lab=1.2)
      legend(x=(ncol(grn)+0.2)*1.05,y=ymax*0.8,xjust = 0,yjust=0.5, bty="o",
             legend=label, col=colcode,xpd=TRUE,pch=15,cex=0.8)
      plot(idx,red,col=as.vector(cc$Color),ylim=c(0,ymax),main=paste(ctype,
                                                                     " Red",sep=""),bty="o",xlab="Sample",ylab="Intensity",cex.lab=1.2)
      legend(x=(ncol(red)+0.2)*1.05,y=ymax*0.8,xjust = 0,yjust=0.5, bty="o",
             legend=label, col=colcode,xpd=TRUE,pch=15,cex=0.8)
    }
    else
    {
      par(mar=c(5, 4, 4, 8.2))
      idx <- t(replicate(nrow(red),1:ncol(red)))
      if(ctype=="STAINING"){ # fixed so that DNP(20k) and Biotin(5k) - ignored by illumina EPIC table - don't appear
        plot(idx,grn,col=as.vector(cc$Color),ylim=c(0,ymax),main=paste(ctype,
                                                                       " Green",sep=""),bty="o",xlab="Sample",ylab="Intensity",cex.lab=1.2)
        legend(x=ncol(grn)*1.05,y=ymax*0.8,xjust = 0,yjust=0.5, bty="o",legend=
                 as.vector(cc$ExtendedType)[3:6],col=as.vector(cc$Color)[3:6],xpd=TRUE,pch=15,cex=0.8)
        plot(idx,red,col=as.vector(cc$Color),ylim=c(0,ymax),main=paste(ctype,
                                                                       " Red",sep=""),bty="o",xlab="Sample",ylab="Intensity",cex.lab=1.2)
        legend(x=ncol(red)*1.045,y=ymax*0.8,xjust = 0,yjust=0.5, bty="o",legend=
                 as.vector(cc$ExtendedType)[3:6],col=as.vector(cc$Color)[3:6],xpd=TRUE,pch=15,cex=0.8)
      }
      else{
        
        plot(idx,grn,col=as.vector(cc$Color),ylim=c(0,ymax),main=paste(ctype,
                                                                       " Green",sep=""),bty="o",xlab="Sample",ylab="Intensity",cex.lab=1.2)
        legend(x=ncol(grn)*1.05,y=ymax*0.8,xjust = 0,yjust=0.5, bty="o",legend=
                 as.vector(cc$ExtendedType),col=as.vector(cc$Color),xpd=TRUE,pch=15,cex=0.8)
        plot(idx,red,col=as.vector(cc$Color),ylim=c(0,ymax),main=paste(ctype,
                                                                       " Red",sep=""),bty="o",xlab="Sample",ylab="Intensity",cex.lab=1.2)
        legend(x=ncol(red)*1.045,y=ymax*0.8,xjust = 0,yjust=0.5, bty="o",legend=
                 as.vector(cc$ExtendedType),col=as.vector(cc$Color),xpd=TRUE,pch=15,cex=0.8)
      }
    }
    dev.off()
  }
  }
