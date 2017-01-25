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
  print(paste(" Testing contrast", x,sep=" ")) # message of progress
  beta_sub=beta_diffcells[ , grepl(paste(x, collapse = "|") , colnames( beta_diffcells ) ) ] # subset beta on contrast stages
  myDMR_timecourse[[s]] <- champ.DMR(beta=beta_sub, maxGap=900, cores=4,
                                pheno=design_timecourse[[s]]$group,method="Bumphunter",arraytype = "EPIC") # call function
}

for(s in stages[2:length(stages)]){
 write.csv(myDMR_timecourse[[s]],paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMR/",s,"_timecourse_DMR_CpGs_within_900bp.csv",sep=""), col.names=T,row.names=T, quote=F)
}

DMR.GUI(DMR=myDMR_timecourse[["EN7"]],beta=beta_sub,pheno=design_timecourse$EN7$group,arraytype = "EPIC",runDMP = F) 
# not doing what it's supposed to do

