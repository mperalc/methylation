### plot methylation and differential methylation stats across genome 


################## setup ###################

library(ggplot2)
library(readr)
library(reshape2)
library(ChAMPdata)
library(plyr)

currentDate <- Sys.Date() # to save date in name of output files


################################## read in files #############################
# beta values for all probes
beta= read_csv( file="/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/quantile_normalised_beta_detP_0.01_nocrossreact.csv")
beta=as.data.frame(beta)
rownames(beta)=beta$X1
beta=beta[2:ncol(beta)]
# 
# beta[beta<=0.001]<-0.001 # replacing extreme values
# beta[beta>=0.999]<-0.999
# 
# M <- log((beta/(1-beta)),2) # calculates M values




stages=c("iPSC","DE","PGT","PFG","PE","EP","EN6","EN7")


mean_betas=colMeans(beta)   # all more of less the same
# Better divided in methylated (>0.8), partially methylated (0.6-0.8), mid(0.4-0.6), partially unmethylated(0.2-0.4), unmethylated (<0.2)
meth_ranks=c("M","PM","P","PU","U")

general_methylation=as.data.frame(matrix(nrow = 5,ncol=ncol(beta)))
colnames(general_methylation)=colnames(beta)
rownames(general_methylation)=meth_ranks

probes=nrow(beta)  # total number of probes around 800,000 (after filtering)

for(c in colnames(beta)){
  
  general_methylation[1,c]=sum(beta[,c]>0.8)
  general_methylation[2,c]=sum((0.6<beta[,c]) & (beta[,c]<0.8))
  general_methylation[3,c]=sum((0.4<beta[,c]) & (beta[,c]<0.6))
  general_methylation[4,c]=sum((0.2<beta[,c]) & (beta[,c]<0.4))
  general_methylation[5,c]=sum(beta[,c]<0.2)
  
}

general_methylation=cbind(rownames(general_methylation),general_methylation)
colnames(general_methylation)[1]="methylation_status"


plot_general_methylation=melt(data = general_methylation)
plot_general_methylation$methylation_status=factor(plot_general_methylation$methylation_status,levels = c("M","PM","P","PU","U"))


p1 <- ggplot(as.data.frame(plot_general_methylation), aes(x=variable,y=value,fill=methylation_status)) + geom_bar(stat="identity") +
  ggtitle("Methylation") +
  ylab ("Amount of probes") +
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.text.y=element_text(size=8,face="bold"),axis.text.x=element_text(size=7,face="bold",angle = 90, hjust = 1,vjust = 0.5),
        panel.border=element_blank(),plot.title = element_text(size=9,face="bold"),
        axis.title.y=element_text(size=8,face="bold"),axis.title.x=element_blank(),
        legend.text = element_text(size=7,face="bold"),legend.title = element_blank(),
        legend.position = "bottom",legend.direction = "horizontal") 



ggsave(filename=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/methylation_levels",currentDate,".jpg",sep=""),p1,width=4,height=4,units="in",dpi=300)


## probe features

data(probe.features.epic)

beta_complete=probe.features[rownames(beta),c(1,6,7)]
beta_complete=cbind(beta,beta_complete)



# select methylated as beta value > 0.6
# select probe names for each sample
# after that, plot each sample separately according to the probe features


# for each CpG island annotation

island_annotation=as.data.frame(matrix(nrow = 4,ncol=ncol(beta)))
colnames(island_annotation)=colnames(beta)
rownames(island_annotation)=unique(probe.features$cgi)


for(c in colnames(beta)){   # in % of total
  
  island_annotation[1,c]=(sum((0.6<beta_complete[,c]) & beta_complete$cgi=="shore") /sum(0.6<beta_complete[,c])) * 100
  island_annotation[2,c]=(sum((0.6<beta_complete[,c]) & beta_complete$cgi=="opensea") /sum(0.6<beta_complete[,c])) * 100
  island_annotation[3,c]=(sum((0.6<beta_complete[,c]) & beta_complete$cgi=="island") /sum(0.6<beta_complete[,c])) * 100
  island_annotation[4,c]=(sum((0.6<beta_complete[,c]) & beta_complete$cgi=="shelf") /sum(0.6<beta_complete[,c])) * 100
  
}

island_annotation=cbind(rownames(island_annotation),island_annotation)
colnames(island_annotation)[1]="island_annotation"



plot_island_annotation=melt(data = island_annotation)
plot_island_annotation$island_annotation=factor(plot_island_annotation$island_annotation,levels = c("opensea", "shelf", "shore", "island" ))


p2 <- ggplot(as.data.frame(plot_island_annotation), aes(x=variable,y=value,fill=island_annotation)) + geom_bar(stat="identity") +
  ggtitle("Methylated (B>0.6) probes by CpG island annotation") +
  ylab ("% of probes") +
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.text.y=element_text(size=8,face="bold"),axis.text.x=element_text(size=7,face="bold",angle = 90, hjust = 1,vjust = 0.5),
        panel.border=element_blank(),plot.title = element_text(size=9,face="bold"),
        axis.title.y=element_text(size=8,face="bold"),axis.title.x=element_blank(),
        legend.text = element_text(size=7,face="bold"),legend.title = element_blank(),
        legend.position = "bottom",legend.direction = "horizontal") 


ggsave(filename=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/methylated_B_over0.6_probes_by_CpG_island_annotation",currentDate,".jpg",sep=""),p2,width=4,height=4,units="in",dpi=300)



# for each genomic annotation

genomic_annotation=as.data.frame(matrix(nrow = 8,ncol=ncol(beta)))
colnames(genomic_annotation)=colnames(beta)
rownames(genomic_annotation)=unique(probe.features$feature)


for(c in colnames(beta)){ # in %
  
  genomic_annotation[1,c]=(sum((0.6<beta_complete[,c]) & beta_complete$feature=="TSS1500")/sum(0.6<beta_complete[,c])) * 100
  genomic_annotation[2,c]=(sum((0.6<beta_complete[,c]) & beta_complete$feature=="IGR")  /sum(0.6<beta_complete[,c])) * 100
  genomic_annotation[3,c]=(sum((0.6<beta_complete[,c]) & beta_complete$feature=="Body") /sum(0.6<beta_complete[,c])) * 100
  genomic_annotation[4,c]=(sum((0.6<beta_complete[,c]) & beta_complete$feature=="3'UTR") /sum(0.6<beta_complete[,c])) * 100
  genomic_annotation[5,c]=(sum((0.6<beta_complete[,c]) & beta_complete$feature=="1stExon") /sum(0.6<beta_complete[,c])) * 100
  genomic_annotation[6,c]=(sum((0.6<beta_complete[,c]) & beta_complete$feature=="TSS200") /sum(0.6<beta_complete[,c])) * 100
  genomic_annotation[7,c]=(sum((0.6<beta_complete[,c]) & beta_complete$feature=="5'UTR") /sum(0.6<beta_complete[,c])) * 100
  genomic_annotation[8,c]=(sum((0.6<beta_complete[,c]) & beta_complete$feature=="ExonBnd") /sum(0.6<beta_complete[,c])) * 100
}


genomic_annotation=cbind(rownames( genomic_annotation), genomic_annotation)
colnames(genomic_annotation)[1]="genomic_annotation"



plot_genomic_annotation=melt(data = genomic_annotation)
plot_genomic_annotation$genomic_annotation=factor(plot_genomic_annotation$genomic_annotation,
                                                  levels = c("IGR", "TSS1500", "TSS200", "3'UTR","Body","1stExon","ExonBnd","5'UTR" ))




p3 <- ggplot(as.data.frame(plot_genomic_annotation), aes(x=variable,y=value,fill=genomic_annotation)) + geom_bar(stat="identity") +
  ggtitle("Methylated (B>0.6) probes by genomic annotation") +
  ylab ("% of probes") +
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.text.y=element_text(size=8,face="bold"),axis.text.x=element_text(size=7,face="bold",angle = 90, hjust = 1,vjust = 0.5),
        panel.border=element_blank(),plot.title = element_text(size=9,face="bold"),
        axis.title.y=element_text(size=8,face="bold"),axis.title.x=element_blank(),
        legend.text = element_text(size=7,face="bold"),legend.title = element_blank(),
        legend.position = "bottom",legend.direction = "horizontal") 


ggsave(filename=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/methylated_B_over0.6_probes_by_genomic_annotation",currentDate,".jpg",sep=""),p3,width=5,height=4,units="in",dpi=300)


# DMR and DMP info for all stages

DMR_timecourse <- list()
DMP_timecourse <- list()

for(s in stages[2:length(stages)]){
  DMR_timecourse[[s]]<- as.data.frame(read_csv(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMR/",s,"_timecourse_DMR_CpGs_within_900bp.csv",sep="")))
  rownames(DMR_timecourse[[s]]) <- DMR_timecourse[[s]]$X1
  DMR_timecourse[[s]] = DMR_timecourse[[s]][2:ncol(DMR_timecourse[[s]])]
}

for(s in stages[2:length(stages)]){
  DMP_timecourse[[s]]<- as.data.frame(read_csv(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMP/",s,"_timecourse_DMPs_2017-01-27.csv",sep="")))
  rownames(DMP_timecourse[[s]]) <- DMP_timecourse[[s]]$X1
  DMP_timecourse[[s]] = DMP_timecourse[[s]][2:ncol(DMP_timecourse[[s]])]
}

# plot variation in DMP and DMR with stages
DMP_timecourse_amount=sapply(DMP_timecourse,nrow)
DMR_timecourse_amount=sapply(DMR_timecourse,nrow)

DMP_timecourse_amount=melt(DMP_timecourse_amount)  # melt for ggplot
DMP_timecourse_amount=cbind(DMP_timecourse_amount,rownames(DMP_timecourse_amount))
colnames(DMP_timecourse_amount)[2]="stage_iPSC"
DMP_timecourse_amount$stage_iPSC=factor(DMP_timecourse_amount$stage_iPSC,levels= stages[2:length(stages)])
DMR_timecourse_amount=melt(DMR_timecourse_amount)
DMR_timecourse_amount=cbind(DMR_timecourse_amount,rownames(DMR_timecourse_amount))
colnames(DMR_timecourse_amount)[2]="stage_iPSC"
DMR_timecourse_amount$stage_iPSC=factor(DMR_timecourse_amount$stage_iPSC,levels= stages[2:length(stages)])


p4 <- ggplot(as.data.frame(DMP_timecourse_amount), aes(x=stage_iPSC,y=value)) + geom_bar(stat="identity") +
  ggtitle("Amount of DMPs: contrasts stage - iPSC") +
  ylab ("Number of DMPs") +
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.text.y=element_text(size=12,face="bold"),axis.text.x=element_text(size=12,face="bold",angle = 90, hjust = 1,vjust = 0.5),
        panel.border=element_blank(),plot.title = element_text(size=15,face="bold"),
        axis.title.y=element_text(size=12,face="bold"),axis.title.x=element_blank(),
        legend.text = element_text(size=11,face="bold"),legend.title = element_blank(),
        legend.position = "bottom",legend.direction = "horizontal") 

ggsave(filename=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMPs_timecourse",currentDate,".jpg",sep=""),p4,width=5,height=4,units="in",dpi=300)


p5 <- ggplot(as.data.frame(DMR_timecourse_amount), aes(x=stage_iPSC,y=value)) + geom_bar(stat="identity") +
  ggtitle("Amount of DMRs: contrasts stage - iPSC") +
  ylab ("Number of DMRs") +
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.text.y=element_text(size=12,face="bold"),axis.text.x=element_text(size=12,face="bold",angle = 90, hjust = 1,vjust = 0.5),
        panel.border=element_blank(),plot.title = element_text(size=15,face="bold"),
        axis.title.y=element_text(size=12,face="bold"),axis.title.x=element_blank(),
        legend.text = element_text(size=11,face="bold"),legend.title = element_blank(),
        legend.position = "bottom",legend.direction = "horizontal") 

ggsave(filename=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMRs_timecourse_maxdist900",currentDate,".jpg",sep=""),p5,width=5,height=4,units="in",dpi=300)


# plot variation in DMP with stages and features
for (s in stages[2:length(stages)]){
  DMP_timecourse[[s]]= cbind(DMP_timecourse[[s]],rep(s,nrow(DMP_timecourse[[s]])))
  colnames(DMP_timecourse[[s]])[20]="stage" 
  DMP_timecourse[[s]]=DMP_timecourse[[s]][c(7,12,13,20)]
}


DMP_timecourse_extended=do.call("rbind", DMP_timecourse)

DMP_by_islands=count(DMP_timecourse_extended, c("stage", "cgi"))
DMP_by_islands$cgi=factor(DMP_by_islands$cgi,levels = c("opensea", "shelf", "shore", "island" ))


DMP_by_genomic=count(DMP_timecourse_extended, c("stage", "feature"))
DMP_by_genomic$feature= factor(DMP_by_genomic$feature,
                           levels = c("IGR", "TSS1500", "TSS200", "3'UTR","Body","1stExon","ExonBnd","5'UTR" ))
                           
DMP_by_chr=count(DMP_timecourse_extended, c("stage", "CHR"))
DMP_by_chr$CHR= factor(DMP_by_chr$CHR, levels = c(1:22))

p6 <- ggplot(as.data.frame(DMP_by_islands), aes(x=stage,y=freq,fill=cgi)) + geom_bar(stat="identity") +
  ggtitle(expression(atop("Amount of DMPs in island features:", atop(italic("contrasts stage - iPSC"), "")))) +
  ylab ("Number of DMPs") +
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.text.y=element_text(size=12,face="bold"),axis.text.x=element_text(size=12,face="bold",angle = 90, hjust = 1,vjust = 0.5),
        panel.border=element_blank(),plot.title = element_text(size=15,face="bold"),
        axis.title.y=element_text(size=12,face="bold"),axis.title.x=element_blank(),
        legend.text = element_text(size=11,face="bold"),legend.title = element_blank(),
        legend.position = "bottom",legend.direction = "horizontal") 

ggsave(filename=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMPs_timecourse_amount_in_island_features",currentDate,".jpg",sep=""),p6,width=5,height=4,units="in",dpi=300)


p7 <- ggplot(as.data.frame(DMP_by_genomic), aes(x=stage,y=freq,fill=feature)) + geom_bar(stat="identity") +
  ggtitle(expression(atop("Amount of DMPs in genomic features:", atop(italic("contrasts stage - iPSC"), "")))) +
  ylab ("Number of DMPs") +
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.text.y=element_text(size=12,face="bold"),axis.text.x=element_text(size=12,face="bold",angle = 90, hjust = 1,vjust = 0.5),
        panel.border=element_blank(),plot.title = element_text(size=15,face="bold"),
        axis.title.y=element_text(size=12,face="bold"),axis.title.x=element_blank(),
        legend.text = element_text(size=11,face="bold"),legend.title = element_blank(),
        legend.position = "bottom",legend.direction = "horizontal") 

ggsave(filename=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMPs_timecourse_amount_in_genomic_features",currentDate,".jpg",sep=""),p7,width=5,height=4,units="in",dpi=300)


p8 <- ggplot(as.data.frame(DMP_by_chr), aes(x=stage,y=freq,fill=CHR)) + geom_bar(stat="identity") +
  ggtitle(expression(atop("Amount of DMPs in each chr:", atop(italic("contrasts stage - iPSC"), "")))) +
  ylab ("Number of DMPs") +
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.text.y=element_text(size=12,face="bold"),axis.text.x=element_text(size=12,face="bold",angle = 90, hjust = 1,vjust = 0.5),
        panel.border=element_blank(),plot.title = element_text(size=15,face="bold"),
        axis.title.y=element_text(size=12,face="bold"),axis.title.x=element_blank(),
        legend.text = element_text(size=11,face="bold"),legend.title = element_blank(),
        legend.position = "bottom",legend.direction = "horizontal",
        legend.background = element_blank(),legend.key=element_blank(),
        legend.box.background=element_blank(),
        legend.key.width = unit(0.3, "cm"),legend.key.height = unit(0.3, "cm") ) 

ggsave(filename=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMPs_timecourse_amount_in_chr",currentDate,".jpg",sep=""),p8,width=7,height=7,units="in",dpi=300)


# plot variation in DMR with stages and chr

#######DO


for (s in stages[2:length(stages)]){
  DMR_timecourse[[s]]= cbind(DMP_timecourse[[s]],rep(s,nrow(DMP_timecourse[[s]])))
  colnames(DMR_timecourse[[s]])[20]="stage" 
  DMR_timecourse[[s]]=DMR_timecourse[[s]][c(7,12,13,20)]
}



