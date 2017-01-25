# differential methylation analysis (DMA)

# check what I need and save in QC script

library(ChAMP)
currentDate <- Sys.Date() # to save date in name of output files

# load matrix of normalized beta values
beta= read.csv( "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/quantile_normalised_beta_detP_0.01_nocrossreact.csv", header=T,row.names = 1,check.names=F)

# load matrix with sample information
pD=read.csv("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/samples_info.csv", header=T, row.names=1, check.names = F)

stages=c("iPSC","DE","PGT","PFG","PE","EP","EN6","EN7")
#############



# I want to check first the differentiated cells, so I'll exclude the adult islets:

beta_diffcells=beta[,1:21]
beta_diffcells=as.matrix(beta_diffcells)
# create the design matrix, which contains the comparisons I want to make

####### create different matrices for each stage
# design_peak
group_peak_iPSC=c(rep(1,times=3),rep(-1,times=18))
design_peak=pD[1:21,c(3,4)]
design_peak$group=group_peak

# iPSC vs all others. Change for each peak test
group_timecourse <- list()
  group_timecourse[["DE"]]=c(rep(-1,times=3),rep(1,times=3),rep(0,times=15)) #DE vs iPSC
  group_timecourse[["PGT"]]=c(rep(-1,times=3),rep(0,times=3),rep(1,times=3),rep(0,times=12)) #PGT vs iPSC
  group_timecourse[["PFG"]]=c(rep(-1,times=3),rep(0,times=6),rep(1,times=3),rep(0,times=9)) #PFG vs iPSC
  group_timecourse[["PE"]]=c(rep(-1,times=3),rep(0,times=9),rep(1,times=3),rep(0,times=6)) #PE vs iPSC
  group_timecourse[["EP"]]=c(rep(-1,times=3),rep(0,times=12),rep(1,times=1),rep(0,times=5)) #EP vs iPSC. Just one sample, so not really diff meth
  group_timecourse[["EN6"]]=c(rep(-1,times=3),rep(0,times=13),rep(1,times=3),rep(0,times=2)) #EN vs iPSC
  group_timecourse[["EN7"]]=c(rep(-1,times=3),rep(0,times=16),rep(1,times=2)) #BLC vs iPSC

# iPSC overmethylated regions are the undermethylated regions of the above 7

design_timecourse= pD[1:21,c(3,4)]

design_timecourse <- rep(list(design_timecourse),7)
names(design_timecourse) = stages[2:length(stages)]
for(s in stages[2:length(stages)]){
  design_timecourse[[s]]$group=group_timecourse[[s]] # fil list of df with appropriate design
  design_timecourse[[s]]=design_timecourse[[s]][design_timecourse[[s]]$group!=0,] # select only rows with -1 or 1. Check this is necessary
  
}

# can I check contrasts in this way?? what about iPSC?


#######

# Differentially methylated regions

###################### #############
#summary information of probes:
# do the same with piecharts and %
CpG.GUI(CpG=rownames(beta_diffcells),arraytype="EPIC")

CpG=rownames(beta_diffcells)  # get probe names - CpG regions

data("probe.features.epic") #load probe features from EPIC data
cgi.info <- table(probe.features[CpG,"cgi"])
chromsome.info <- table(probe.features[CpG,"CHR"])
feature.info <- table(probe.features[CpG,"feature"])
type.info <- table(probe.features[CpG,"Type"])

# plot as stacked bar chart

p1 <- ggplot(as.data.frame(cgi.info), aes(x="x",y=Freq,fill=Var1)) + geom_bar(stat="identity") +
  ggtitle("Location of probes in genome") +
  ylab ("Amount of probes") +
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.text.y=element_text(size=12,face="bold"),axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        panel.border=element_blank(),plot.title = element_text(size=15,face="bold"),
        axis.title.y=element_text(size=14,face="bold"),axis.title.x=element_blank(),
        legend.text = element_text(size=11,face="bold"),legend.title = element_blank(),
        legend.position = "bottom",legend.direction = "horizontal") +
geom_hline(yintercept=0,size=1)
ggsave(filename=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/probes_in_genome_CpGisland",currentDate,".jpg",sep=""),p1,width=4,height=4,units="in",dpi=300)

p2 <- ggplot(as.data.frame(feature.info), aes(x="x",y=Freq,fill=Var1)) + geom_bar(stat="identity")+
  ggtitle("Location of probes in genome") +
  ylab ("Amount of probes") +
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.text.y=element_text(size=12,face="bold"),axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        panel.border=element_blank(),plot.title = element_text(size=15,face="bold"),
        axis.title.y=element_text(size=14,face="bold"),axis.title.x=element_blank(),
        legend.text = element_text(size=8,face="bold"),legend.title = element_blank(),
        legend.position = "bottom",legend.direction = "horizontal") +
  geom_hline(yintercept=0,size=1)

ggsave(filename=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/probes_in_genome_gene",currentDate,".jpg",sep=""),p2,width=4,height=4,units="in",dpi=300)


p3 <- ggplot(as.data.frame(type.info), aes(x="x",y=Freq,fill=Var1)) + geom_bar(stat="identity") +
  ggtitle("Type of probes") +
  ylab ("Amount of probes") +
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.text.y=element_text(size=12,face="bold"),axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        panel.border=element_blank(),plot.title = element_text(size=15,face="bold"),
        axis.title.y=element_text(size=14,face="bold"),axis.title.x=element_blank(),
        legend.text = element_text(size=11,face="bold"),legend.title = element_blank(),
        legend.position = "bottom",legend.direction = "horizontal") +
  geom_hline(yintercept=0,size=1)
ggsave(filename=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/probe_types",currentDate,".jpg",sep=""),p3,width=4,height=4,units="in",dpi=300)


chromosome.info=as.data.frame(chromsome.info)[2:23,] # subset chromosome data
chromosome.info$Var1=droplevels(chromosome.info$Var1) # drop unused levels
chromosome.info$Var1=factor(chromosome.info$Var1,levels=c(1:22)) # reorder levels for plot

p4 <- ggplot(chromosome.info, aes(x="x",y=Freq,fill=Var1)) + geom_bar(stat="identity") +
  ggtitle("Distribution of probes among chromosomes") +
  ylab ("Amount of probes") +
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.text.y=element_text(size=12,face="bold"),axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        panel.border=element_blank(),plot.title = element_text(size=12,face="bold"),
        axis.title.y=element_text(size=12,face="bold"),axis.title.x=element_blank(),
        legend.key.size = unit(0.3,"cm"),
        legend.text = element_text(size=8,face="bold"),legend.title = element_blank(),
        legend.position = "bottom",legend.direction = "horizontal") +
  guides(fill=guide_legend(nrow=2)) +
  geom_hline(yintercept=0,size=1)
ggsave(filename=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/probe_in_chr",currentDate,".jpg",sep=""),p4,width=4,height=4,units="in",dpi=300)

#################
QC.GUI(CpG=rownames(beta_diffcells),arraytype="EPIC")


# need design to test 1 sample vs all others (peak/stage-specific) or using contrasts against iPSC(timecourse/across-stages)
myDMR_peak <- champ.DMR(beta=beta_diffcells,pheno=design_peak$group,method="Bumphunter",arraytype = "EPIC")

# ChAMP only tests two groups. Ensure the comparisons are correct!
myDMR_timecourse <- list()
for(s in stages[2:length(stages)]){
  x = c("iPSC",s) # two stages to contrast
  print(paste(" Testing contrast", paste(x, collapse = "|") ,sep=" ")) # message of progress
  beta_sub=beta_diffcells[ , grepl(paste(x, collapse = "|") , colnames( beta_diffcells ) ) ] # subset beta on contrast stages
  myDMR_timecourse[[s]] <- champ.DMR(beta=beta_sub, 
                                     maxGap=900, 
                                     cores=4,
                                     pheno=design_timecourse[[s]]$group,
                                     method="Bumphunter",
                                     arraytype = "EPIC") # call function
}

for(s in stages[2:length(stages)]){
 write.csv(myDMR_timecourse[[s]],paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMR/",s,"_timecourse_DMR_CpGs_within_900bp.csv",sep=""), col.names=T,row.names=T, quote=F)
}

DMR.GUI(DMR=myDMR_timecourse[["EN7"]],beta=beta_sub,pheno=design_timecourse$EN7$group,arraytype = "EPIC",runDMP = F) 
# not doing what it's supposed to do


# Differentially methylated probes


# iPSC vs all others. Change for each peak test
group_timecourse <- list()
group_timecourse[["DE"]]=c(rep("iPSC",times=3),rep("DE",times=3),rep(0,times=15)) #DE vs iPSC
group_timecourse[["PGT"]]=c(rep("iPSC",times=3),rep(0,times=3),rep("PGT",times=3),rep(0,times=12)) #PGT vs iPSC
group_timecourse[["PFG"]]=c(rep("iPSC",times=3),rep(0,times=6),rep("PFG",times=3),rep(0,times=9)) #PFG vs iPSC
group_timecourse[["PE"]]=c(rep("iPSC",times=3),rep(0,times=9),rep("PE",times=3),rep(0,times=6)) #PE vs iPSC
group_timecourse[["EP"]]=c(rep("iPSC",times=3),rep(0,times=12),rep("EP",times=1),rep(0,times=5)) #EP vs iPSC. Just one sample, so not really diff meth
group_timecourse[["EN6"]]=c(rep("iPSC",times=3),rep(0,times=13),rep("EN6",times=3),rep(0,times=2)) #EN vs iPSC
group_timecourse[["EN7"]]=c(rep("iPSC",times=3),rep(0,times=16),rep("EN7",times=2)) #BLC vs iPSC

# iPSC overmethylated regions are the undermethylated regions of the above 7

design_timecourse= pD[1:21,c(3,4)]

design_timecourse <- rep(list(design_timecourse),7)
names(design_timecourse) = stages[2:length(stages)]
for(s in stages[2:length(stages)]){
  design_timecourse[[s]]$group=group_timecourse[[s]] # fil list of df with appropriate design
  design_timecourse[[s]]=design_timecourse[[s]][design_timecourse[[s]]$group!=0,] # select only rows non 0 vals.
  
}

myDMP_timecourse <- list()
for(s in stages[2:length(stages)]){
  x = c("iPSC",s) # two stages to contrast
  print(paste(" Testing contrast", paste(x, collapse = "|") ,sep=" ")) # message of progress
  beta_sub=beta_diffcells[ , grepl(paste(x, collapse = "|") , colnames( beta_diffcells ) ) ] # subset beta on contrast stages
  
  if(is.null(dim(beta_diffcells[ , grepl(x[2] , colnames( beta_diffcells ))]))){
    # if second stage (not iPSC) only has one column, do this alternative version of champ.DMP
    # Here there's no average beta value for each CpG in the stage with one sample, just its only beta value
    # As there are not replicates, this is not a true differential methylation analysis. Hence the message:
    message(paste("The stage",x[2],"has only one sample. Doing DMP variation, but do not trust results.",sep=" "))
    message("[===========================]")
    message("[<<<<< ChAMP.DMP VARIATION STARTING >>>>>]")
    message("-----------------------------")
    
    ### setup
    beta=beta_sub
    pheno=factor(design_timecourse[[s]]$group,levels=c("iPSC",s))
    arraytype = "EPIC"
    message(paste("The array type is",arraytype,sep=" "))
    adjPVal = 0.05
    message(paste("The adjusted p-val threshold to report is",adjPVal,sep=" "))
    adjust.method = "BH"
    message(paste("The adjustment method for multiple testing is",adjust.method,sep=" "))
    
    
    #end of setup 
    
    message("\n<< Your pheno information contains following groups. >>")
    sapply(unique(pheno),function(x) message("<",x,">:",sum(pheno==x)," samples."))
    message("[The power of statistics analysis on groups contain very few samples may not be strong.]")
    
    message("You did not assign compare groups. The first two groups: <",unique(pheno)[1],"> and <",unique(pheno)[2],">, will be compared automatically.")
    compare.group <- unique(pheno)[1:2]
    
    p <- pheno[which(pheno %in% compare.group)]
    beta <- beta[,which(pheno %in% compare.group)]
    design <- model.matrix( ~ 0 + p)
    # contrast.matrix requires the initial groups from pheno to be recorded as factors, preferably names (for example, stages)
    contrast.matrix <- makeContrasts(contrasts=paste(colnames(design)[2:1],collapse="-"), levels=colnames(design))
    message("\n<< Contrast Matrix >>")
    print(contrast.matrix)
    
    message("\n<< All beta, pheno and model are prepared successfully. >>")
    
    fit <- lmFit(beta, design)
    fit2 <- contrasts.fit(fit,contrast.matrix)
    tryCatch(fit3 <- eBayes(fit2),
             warning=function(w) 
             {
               stop("limma failed, No sample variance.\n")
             }) # if the contrast matrix is not correct, DMP function will fail here
    
    DMP <- topTable(fit3,coef=1,number=nrow(beta),adjust.method=adjust.method,p.value=adjPVal)
    message("You have found ",sum(DMP$adj.P.Val <= adjPVal), " significant MVPs with a ",adjust.method," adjusted P-value below ", adjPVal,".")
    message("\n<< Calculate DMP successfully. >>")
    
    if(arraytype == "EPIC") data(probe.features.epic) else data(probe.features)
    com.idx <- intersect(rownames(DMP),rownames(probe.features))
    avg <-  cbind(rowMeans(beta[com.idx,which(p==compare.group[1])]),beta[com.idx,which(p==compare.group[2])])
    avg <- cbind(avg,avg[,2]-avg[,1])
    colnames(avg) <- c(paste(compare.group,"AVG",sep="_"),"deltaBeta")
    DMP <- data.frame(DMP[com.idx,],avg,probe.features[com.idx,])
    
    message("[<<<<<< ChAMP.DMP VARIATION ENDED, JUST AS ALL THINGS END IN LIFE >>>>>>]")
    message("[===========================]")
    
  }else{
    
    myDMP_timecourse[[s]] <- champ.DMP(beta = beta_sub, 
                                       pheno=factor(design_timecourse[[s]]$group,levels=c("iPSC",s)), #levels argument necesary, otherwise comparisons will be determined by alphabetical order
                                       arraytype = "EPIC"   ) # call function
  }
  
}

# failing on the calculation of average value of EP, because it has only one column
# write exception for it, but don't trust (as without replicates there's no differential methylation!)

