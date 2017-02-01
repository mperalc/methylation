# plot DMPs and DMRs


################## setup ###################

library(ggplot2)
library(readr)
library(reshape2)
library(grid)
library(gridExtra)
library(gtable)
library(ChAMPdata)
library(digest)

currentDate <- Sys.Date() # to save date in name of output files


################################## read in files #############################
# beta values for all probes
beta= read.csv( "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/quantile_normalised_beta_detP_0.01_nocrossreact.csv", header=T,row.names = 1,check.names=F)
# deltaBeta values of all probes (stage (n) - iPSC), and beta values averages for every stage
deltaBetas <- as.data.frame(read_csv(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/deltaBeta_allprobes.csv"))
rownames(deltaBetas) <- deltaBetas$X1
deltaBetas = deltaBetas[2:ncol(deltaBetas)]

## probe features

data(probe.features.epic)


stages=c("iPSC","DE","PGT","PFG","PE","EP","EN6","EN7")

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



#######################
############### region to plot ##########

region <- c(11, 2450360,2540221)  # Chr, position start, position finish
timecourse_stage = "EN7"  # stage from timecourse to plot (meth values compared to iPSC stage)



######### to re-run with each change of region ########


x = c("iPSC",timecourse_stage) # two stages to subset
beta_sub=beta[ , grepl(paste(x, collapse = "|") , colnames(beta) ) ] # subset beta on iPSC and the stage selected

beta_sub = cbind(beta_sub,probe.features[rownames(beta_sub),c("CHR","MAPINFO","gene","feature","cgi")]) # get relevant info from probes


#Data for first plot (probes and DMR segment)
plot_betas=beta_sub[which(beta_sub$CHR == region[1] & beta_sub$MAPINFO > region[2] & beta_sub$MAPINFO < region[3]) , ]
plot_betas=melt(plot_betas,id.vars = c("CHR","MAPINFO","gene","feature","cgi") )
plot_betas$variable=gsub( "-.*$", "", plot_betas$variable) # take out everything after "-" (sample name)

###check meaning of columns in DMR!!!!
#fix condition to select specific row
selected_DMR=DMR_timecourse[[timecourse_stage]][which(DMR_timecourse[[timecourse_stage]]$BumphunterDMR.seqnames==paste("chr",region[1],sep="") & DMR_timecourse[[timecourse_stage]]$BumphunterDMR.start > region[2] & DMR_timecourse[[timecourse_stage]]$BumphunterDMR.end < region[3]),]

if(selected_DMR$BumphunterDMR.value>0){  # selecting colour of segment based on over/under methylation. CHANGE!!
  color_DMR <- "red" #over-methylated in red
}else{ #under-methylated in blue
  color_DMR <- "blue"
  
}

scaleFUN <- function(x) sprintf("%.1f", x)  # to round y axis to 1 decimal place



library(RColorBrewer)
myColors <- brewer.pal(8,"Set2") # select colors from pallette
myColors <- myColors[c(1,8)] # green and dark grey
names(myColors) <- levels(plot_betas$variable) # assign levels ("stage x" and "iPSC")
colScale <- scale_colour_manual(name = "variable",values = myColors)

p1 <- ggplot(plot_betas, aes(x = MAPINFO,y=value,group=variable)) 

  if(nrow(selected_DMR)!=0){
   p1 <- p1 + geom_segment(show.legend = F,aes(x=selected_DMR$BumphunterDMR.start,xend=selected_DMR$BumphunterDMR.end,y=0,yend=0),col=color_DMR)  # color outside aes doesn't include it in legend
  }

p1 <- p1 + geom_point(aes(col=variable)) + # this with colScale sets colors to variables and includes in label
   colScale +
 scale_y_continuous(labels=scaleFUN) +  # rounding y axis
  ggtitle(paste("DMR",timecourse_stage,"vs iPSC",sep=" ")) +
  xlab (paste("Location in chr",region[1],sep=" ")) +
  ylab ("Beta value") +
  expand_limits(y = c(0,1)) +
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), 
        panel.border=element_rect(size=2),axis.text.y=element_text(size=12,face="bold"),
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.title.y=element_text(size=14,face="bold"),plot.title = element_text(size=16,face="bold"),
        legend.text = element_text(size=11,face="bold"),legend.title = element_blank(),
        legend.position = "bottom",legend.direction = "horizontal",
        legend.background = element_blank(),legend.key=element_blank(),
        legend.box.background=element_rect(size=0.5),legend.margin=unit(-0.5, "cm"),
        legend.key.width = unit(0, "cm") ) 
 
# to save legend:
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
# 
 mylegend<-g_legend(p1)
 


 
 ### add gene annotation
 
 library(Homo.sapiens)
 
 library(dplyr)
 library(biomaRt)
 
 ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") 
 # I use genome coordinates from GRCh37 (not GRCh38, which is the latest release). Keep that in mind.
 
 
 
 chr.region=paste(region, collapse = ":")
 filterlist <- list(chr.region,"protein_coding","lincRNA")  
 
 
 
 results <- getBM(attributes = c('ensembl_gene_id','entrezgene', 'external_gene_name', 'gene_biotype', 'chromosome_name', "start_position", "end_position"),
                  filters = c("chromosomal_region","biotype"),values = filterlist, mart = ensembl)
 
 results_refseqmRNA <- getBM(attributes=c('refseq_mrna','external_gene_name',"chromosome_name","transcript_start","transcript_end","transcription_start_site"), filters = c("chromosomal_region","biotype"), values = filterlist, mart = ensembl)
 
 results_refseqmRNA[which(results_refseqmRNA$refseq_mrna==""),1] <- NA   # changing empty vectors of first column to NA, to remove them later
 
 results_refseqmRNA=results_refseqmRNA[which(!is.na(results_refseqmRNA$refseq_mrna)),]  # taking out non refseq mrnas
 
 tss=unique(results_refseqmRNA$transcript_start)  # one for each gene
 
 ##########fix this tss thing
 
 
 # plot white y axis
 # plot white dots
 # plot geom segment as results or as region,  depending on whether it starts or ends outside plotting region
 
 myColors <- brewer.pal(8,"Dark2") # select colors from pallette
 myColors <- myColors 
 names(myColors) <- levels(results$external_gene_name) # assign levels ("stage x" and "iPSC")
 colScale <- scale_colour_manual(name = "external_gene_name",values = myColors)
 
 p2 <-   ggplot(plot_betas, aes(x = MAPINFO,y=value)) + geom_point(col="white") + colScale
          
   
    values <- nrow(results)   # how many genes do I have
 for (i in 1:values) {  # sequentially plot new genes. Doesn't work for more than 8
   p2 <- p2 + geom_segment(data=results[i,], 
                           show.legend = T,
                           x=results[i,6],xend=results[i,7], 
                           y=(mean(plot_betas$value)*sin(i/4))+0.1,yend=(mean(plot_betas$value)*sin(i/4))+0.1,  # oscilate around 0.5
                           aes(col=external_gene_name)) +
     geom_segment(show.legend = F,x=tss,xend=tss+10,
                  y=(mean(plot_betas$value)*sin(i/4))+0.1-0.08,yend=(mean(plot_betas$value)*sin(i/4))+0.1+0.08)  # making a tick
            
              
    }
   
  p2 <- p2 + ylab ("Genes") +
         scale_y_continuous(labels=scaleFUN) +  # rounding y axis
         expand_limits(y = c(0,1)) +
         theme_bw() + 
         theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
         panel.border=element_blank(),axis.text.y=element_text(size=12,face="bold",colour="white"),
         axis.title.x = element_blank(), axis.text.x =element_blank(),
         axis.ticks.y = element_blank(),axis.title.y=element_text(size=14,face="bold"),
         plot.title = element_text(size=16,face="bold"),
         legend.text = element_text(size=11,face="bold"),legend.title = element_blank(),
         legend.position = "bottom",legend.direction = "horizontal",
         legend.background = element_blank(),legend.key=element_blank(),
         legend.box.background=element_blank(),legend.margin=unit(-0.5, "cm"),
         legend.key.width = unit(0.5, "cm") ) 
 
   
   mylegend2<-g_legend(p2)
   # 
  
 
#Data for deltabeta
plot_deltabeta=deltaBetas[ ,grepl(timecourse_stage, colnames(deltaBetas))] # subset beta on iPSC and the stage selected

plot_deltabeta = cbind(plot_deltabeta ,beta_sub[,(ncol(beta_sub)-4):ncol(beta_sub)]) # join location of probes and deltabeta for selected stage
plot_deltabeta=plot_deltabeta[which(plot_deltabeta$CHR == region[1] & plot_deltabeta$MAPINFO > region[2] & plot_deltabeta$MAPINFO < region[3]) ,]
plot_deltabeta=plot_deltabeta[,c(2:ncol(plot_deltabeta))]
colnames(plot_deltabeta)[1]="deltaBeta"  

#dynamic limits:
#if max deltabeta of the region >0.1, pick max=(max + 0.1) and min -0.15
# if min deltabeta <-0.1, pick min=(min -0.1) and max 0.15
limits_deltabeta=c(-1.0,1.0)
# if(max(plot_deltabeta[1])>=0.1){
#   limits_deltabeta[2]=max(plot_deltabeta[1]) + 0.1
#   }
# if(min(plot_deltabeta[1])<=-0.1){
#   limits_deltabeta[1]=min(plot_deltabeta[1]) - 0.1
#   }


p3 <- ggplot(plot_deltabeta, aes(x = MAPINFO,y=deltaBeta)) + geom_line() + theme_minimal() + 
  xlab (paste("Location in chr",region[1],sep=" ")) +
  ylab ("delta Beta") +
  scale_y_continuous(labels=scaleFUN) +
  expand_limits(y = limits_deltabeta) +
  geom_hline(yintercept=0,linetype="dashed",size=0.1) +
  geom_hline(yintercept=0.1,col="red",size=0.1) +
  geom_hline(yintercept=-0.1,col="red",size=0.1) +
  theme_bw() +
  theme(axis.text=element_text(size=9.50,face="bold"),  # need slightly smaller text on y axis to compensate for the space taken by the "-"
        panel.border=element_rect(size=2,colour = "white"),  # otherwise the ticks on the x axis of all plots will be misaligned
        axis.title=element_text(size=14,face="bold"))

   
        
        
 


# merge all plots
if(nrow(results)<=3){
  size_p2=1.2
}else{
  size_p2=2.2
}

p4 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),ncol=1),
                   arrangeGrob(p2 + theme(legend.position="none"),ncol=1),
                   arrangeGrob(p3 + theme(legend.position="none"),ncol=1),
                   arrangeGrob(mylegend,mylegend2,ncol=2),
                   nrow=4,heights=c(7,size_p2,4.5,1)) 

# blank grob? blank<-rectGrob(gp=gpar(col="white"))

ggsave(plot = p4,
       device="png",
       filename = paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/Methylation_EPIC/DMR_and_DMP_timecourse_iPSC_vs", timecourse_stage,
                        paste(results$external_gene_name, collapse = "_"),
                        currentDate,".png",sep="_"),
       width = 18, height = 6, units = "in",
       dpi = 400)


# table with DMPs

selected_DMP=DMP_timecourse[[timecourse_stage]][which(DMP_timecourse[[timecourse_stage]]$CHR==region[1] & DMP_timecourse[[timecourse_stage]]$MAPINFO > region[2] & DMP_timecourse[[timecourse_stage]]$MAPINFO < region[3]),]
capture.output(selected_DMP, 
               file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMP_timecourse_iPSC_vs",timecourse_stage,
                          paste(results$external_gene_name, collapse = "_"),
                          currentDate,".txt",sep="_")) # to be able to check SNPs

# taking out SNPs column, that makes the .csv file unreadable

selected_DMP2=selected_DMP[,c(1:15,17:19)]

write.csv(selected_DMP2,
          paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMP_timecourse_iPSC_vs",timecourse_stage,
                paste(results$external_gene_name, collapse = "_"),
                currentDate,".csv",sep="_"),
          row.names=T, quote=F)

write.csv(selected_DMP,
          paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMP_timecourse_iPSC_vs",timecourse_stage,
                paste(results$external_gene_name, collapse = "_"),
                currentDate,"plusSNPs.csv",sep="_"),
          row.names=T, quote=F)
