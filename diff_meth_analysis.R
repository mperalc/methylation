# differential methylation analysis (DMA)

# check what I need and save in QC script

library(ChAMP)
currentDate <- Sys.Date() # to save date in name of output files


beta= read.csv( "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/quantile_normalised_beta_detP_0.01_nocrossreact.csv", header=T,row.names = 1,check.names=F)

pD=read.csv("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/samples_info.csv", header=T, row.names=1, check.names = F)

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

group_timecourse=c(rep(-1,times=3),rep(1,times=3),rep(0,times=15)) #DE vs iPSC
design_timecourse= pD[1:21,c(3,4)]
design_timecourse$group=group_timecourse

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
myDMR_timecourse <- champ.DMR(beta=beta_diffcells[,c(1:6)],pheno=design_timecourse$group[1:6],method="Bumphunter",arraytype = "EPIC")

DMR.GUI(DMR=myDMR_timecourse,beta=beta_diffcells,pheno=design_timecourse$group,arraytype = "EPIC") # not doing what it's supposed to do

