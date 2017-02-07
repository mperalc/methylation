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
diff_type="timecourse"   # type of differential analysis: peak or timecourse
sample_type="islets"   # islets or differentiated cells ("diff")
stages=c("iPSC","DE","PGT","PFG","PE","EP","EN6","EN7")
islets=c("EN7","islets","islets-EN7")


################################## read in files #############################
# beta values for all probes
beta= read.csv( "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/quantile_normalised_beta_detP_0.01_nocrossreact.csv", header=T,row.names = 1,check.names=F)
data(probe.features.epic)  ## probe features

# deltabeta, DMR and DMP info for all stages
if(diff_type=="timecourse"){
  if(sample_type=="diff"){
     # deltaBeta values of all probes (stage (n) - iPSC), and beta values averages for every stage
    deltaBetas <- as.data.frame(read_csv(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/deltaBeta_allprobes.csv"))
    rownames(deltaBetas) <- deltaBetas$X1
    deltaBetas = deltaBetas[2:ncol(deltaBetas)]
    
    DMR_timecourse <- list()
    DMP_timecourse <- list()
    
    for(s in stages[2:length(stages)]){
      DMR_timecourse[[s]]<- as.data.frame(read_csv(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMR/",s,"_timecourse_DMR_CpGs_within_900bp.csv",sep="")))
      rownames(DMR_timecourse[[s]]) <- DMR_timecourse[[s]]$X1
      DMR_timecourse[[s]] = DMR_timecourse[[s]][2:ncol(DMR_timecourse[[s]])]
      
      DMP_timecourse[[s]]<- as.data.frame(read_csv(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMP/",s,"_timecourse_DMPs_2017-01-27.csv",sep="")))
      rownames(DMP_timecourse[[s]]) <- DMP_timecourse[[s]]$X1
      DMP_timecourse[[s]] = DMP_timecourse[[s]][2:ncol(DMP_timecourse[[s]])]
    }
  }
  if(sample_type=="islets"){
     # deltaBetaand beta values averages for every stage
    deltaBetas <- as.data.frame(read_csv(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/deltaBeta_allprobes_islets.csv"))
    rownames(deltaBetas) <- deltaBetas$X1
    deltaBetas = deltaBetas[2:ncol(deltaBetas)]
    
    DMR_timecourse <- list()
    DMP_timecourse <- list()
    
    for(s in islets){
      DMR_timecourse[[s]]<- as.data.frame(read_csv(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMR/",s,"_timecourse_DMR_CpGs_within_900bp_2017-02-07.csv",sep="")))
      rownames(DMR_timecourse[[s]]) <- DMR_timecourse[[s]]$X1
      DMR_timecourse[[s]] = DMR_timecourse[[s]][2:ncol(DMR_timecourse[[s]])]
      
      DMP_timecourse[[s]]<- as.data.frame(read_csv(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMP/",s,"_timecourse_DMPs_2017-02-07.csv",sep="")))
      rownames(DMP_timecourse[[s]]) <- DMP_timecourse[[s]]$X1
      DMP_timecourse[[s]] = DMP_timecourse[[s]][2:ncol(DMP_timecourse[[s]])]
    }
  }
}


#######################
############### region to plot ##########

region <- c(11, 2150000,2185000)  # Chr, position start, position finish


message("--------------------------------------------------")
message("---------------STARTING ANALYSIS------------------")
message("--------------------------------------------------")

if(sample_type=="diff"){
  # timecourse_stage="EN7"
  for(timecourse_stage in stages[2:length(stages)]){  # stage from timecourse to plot (meth values compared to iPSC stage)


######### to re-run with each change of region ########


x = c("iPSC",timecourse_stage) # two stages to subset

message(paste("---------------", x[1]," vs ", x[2], "------------------",sep=""))

beta_sub=beta[ , grepl(paste(x, collapse = "|") , colnames(beta) ) ] # subset beta on iPSC and the stage selected

beta_sub = cbind(beta_sub,probe.features[rownames(beta_sub),c("CHR","MAPINFO","gene","feature","cgi")]) # get relevant info from probes


message("---------------Selecting probes------------------")

#Data for first plot (probes and DMR segment)
plot_betas=beta_sub[which(beta_sub$CHR == region[1] & beta_sub$MAPINFO > region[2] & beta_sub$MAPINFO < region[3]) , ]
plot_betas=melt(plot_betas,id.vars = c("CHR","MAPINFO","gene","feature","cgi") )
plot_betas$variable=gsub( "-.*$", "", plot_betas$variable) # take out everything after "-" (sample name)
plot_betas$variable=factor(plot_betas$variable,levels=c("iPSC",timecourse_stage))

###check meaning of columns in DMR!!!!
#fix condition to select specific row
message("---------------selecting DMRs------------------")

selected_DMR=DMR_timecourse[[timecourse_stage]][which(DMR_timecourse[[timecourse_stage]]$BumphunterDMR.seqnames==paste("chr",region[1],sep="") & DMR_timecourse[[timecourse_stage]]$BumphunterDMR.start > region[2] & DMR_timecourse[[timecourse_stage]]$BumphunterDMR.end < region[3]),]
if(nrow(selected_DMR)==0){
  message("---------------no DMRs found------------------")
 
}else{
  if(selected_DMR$BumphunterDMR.value>0){  # selecting colour of segment based on over/under methylation. 
    message("---------------found overmethylated region------------------")
    color_DMR <- "red" #over-methylated in red
  }else{ #under-methylated in blue
    message("---------------found undermethylated region------------------")
    color_DMR <- "blue"
  }
}

message("---------------Plotting probes------------------")


scaleFUN <- function(x) sprintf("%.1f", x)  # to round y axis to 1 decimal place

library(RColorBrewer)
myColors <- brewer.pal(8,"Set2") # select colors from pallette
myColors <- myColors[c(8,1)] # green and dark grey
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
 


 message("---------------adding gene annotation------------------")
 
 ### add gene annotation
 
 library(Homo.sapiens)
 
 library(dplyr)
 library(biomaRt)
 
 ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") 
 # I use genome coordinates from GRCh37 (not GRCh38, which is the latest release). Keep that in mind.
 
 
 
 chr.region=paste(region, collapse = ":")
 filterlist <- list(chr.region)  
 
 
 
 results <- getBM(attributes = c('ensembl_gene_id','entrezgene', 'external_gene_name', 'gene_biotype', 'chromosome_name', "start_position", "end_position","transcription_start_site"),
                  filters = c("chromosomal_region"),values = filterlist, mart = ensembl)
 
 # results_refseqmRNA <- getBM(attributes=c('refseq_mrna','external_gene_name',"chromosome_name","transcript_start","transcript_end","transcription_start_site"), filters = c("chromosomal_region","biotype"), values = filterlist, mart = ensembl)
 # 
 # results_refseqmRNA[which(results_refseqmRNA$refseq_mrna==""),1] <- NA   # changing empty vectors of first column to NA, to remove them later
 # 
 # results_refseqmRNA=results_refseqmRNA[which(!is.na(results_refseqmRNA$refseq_mrna)),]  # taking out non refseq mrnas
 # 
 # tss=unique(results_refseqmRNA$transcript_start)  
 
 ####### if it throws errors of aesthetics, reduce gene size:  ##############
 
 results=results[which(results$external_gene_name=="SNRPN"),]
 
 ##################
 
 results$external_gene_name=as.factor(results$external_gene_name)
 results_tss <- melt(results[c(3,8)])
 results_tss <- unique(results_tss)
 results_tss_list <- split(results_tss, results_tss$external_gene_name) # make a list for each gene

 message("---------------plotting genes------------------")
 
 # plot white y axis
 # plot white dots
 # plot geom segment as results or as region,  depending on whether it starts or ends outside plotting region
 
 
 myColors <- brewer.pal(8,"Dark2") # select colors from pallette
 names(myColors) <- levels(results$external_gene_name) # assign levels ("stage x" and "iPSC")
 colScale <- scale_colour_manual(name = "external_gene_name",values = myColors)
 
 p2 <-   ggplot(plot_betas, aes(x = MAPINFO,y=value)) + geom_point(col="white") + colScale
          
   
    values <- length(unique(results$external_gene_name))   # how many genes do I have
    n=1  # count for the location of the lines
    for (i in unique(results$external_gene_name)) {  # sequentially plot new genes. Doesn't work for more than 8
      p2 <- p2 + geom_segment(data=results[which(results$external_gene_name==i),], 
                              show.legend = T,
                              x=unique(results[which(results$external_gene_name==i),6]),
                              xend=unique(results[which(results$external_gene_name==i),7]), 
                              y=(mean(plot_betas$value)*sin(n/4))+0.1,yend=(mean(plot_betas$value)*sin(n/4))+0.1,  # oscilate around 0.5
                              aes(col=external_gene_name))
      print(n)
      print(i)
      
      for (t in results_tss_list[[i]]$value) {  # plot all tss for each gene along the same line
        p2 <- p2 +geom_segment(show.legend = F,x=t,xend=t,
                               y=(mean(plot_betas$value)*sin(n/4))+0.1-0.08,yend=(mean(plot_betas$value)*sin(n/4))+0.1+0.08)  # making a tick
        
      }
      n=n+1
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
  
   message("---------------plotting delta beta values------------------")
   
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

message("---------------arranging all plots------------------")

p4 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),ncol=1),
                   arrangeGrob(p2 + theme(legend.position="none"),ncol=1),
                   arrangeGrob(p3 + theme(legend.position="none"),ncol=1),
                   arrangeGrob(mylegend,mylegend2,ncol=2),
                   nrow=4,heights=c(7,size_p2,4.5,1)) 

# blank grob? blank<-rectGrob(gp=gpar(col="white"))

message("---------------saving all plots------------------")


ggsave(plot = p4,
       device="png",
       filename = paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/Methylation_EPIC/DMR_and_DMP_timecourse_iPSC_vs", timecourse_stage,
                        paste(unique(results$external_gene_name), collapse = "_"),
                        paste(region, collapse = "_"),
                        currentDate,".png",sep="_"),
       width = 18, height = 6, units = "in",
       dpi = 400)


message("---------------saving results in tables------------------")

# table with DMPs

selected_DMP=DMP_timecourse[[timecourse_stage]][which(DMP_timecourse[[timecourse_stage]]$CHR==region[1] & DMP_timecourse[[timecourse_stage]]$MAPINFO > region[2] & DMP_timecourse[[timecourse_stage]]$MAPINFO < region[3]),]
capture.output(selected_DMP, 
               file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMP_timecourse_iPSC_vs",timecourse_stage,
                          paste(unique(results$external_gene_name), collapse = "_"),
                          paste(region, collapse = "_"),
                          currentDate,".txt",sep="_")) # to be able to check SNPs

# taking out SNPs column, that makes the .csv file unreadable

selected_DMP2=selected_DMP[,c(1:15,17:19)]

write.csv(selected_DMP2,
          paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMP_timecourse_iPSC_vs",timecourse_stage,
                paste(unique(results$external_gene_name), collapse = "_"),
                paste(region, collapse = "_"),
                currentDate,".csv",sep="_"),
          row.names=T, quote=F)

write.csv(selected_DMP,
          paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMP_timecourse_iPSC_vs",timecourse_stage,
                paste(unique(results$external_gene_name), collapse = "_"),
                paste(region, collapse = "_"),
                currentDate,"plusSNPs.csv",sep="_"),
          row.names=T, quote=F)
}
}
if(sample_type=="islets"){
  for(s in islets){
    if(s=="EN7"){
      x = c("iPSC",s) # two stages to contrast
    }
    if(s=="islets"){
      x = c("iPSC","ISL","R") 
    }
    if(s=="islets-EN7"){
      x = c("EN7","ISL","R") 
    }
    
    ######### to re-run with each change of region ########
    
    
    
    
    message(paste("---------------",x, "------------------",sep=""))
    
    beta_sub=beta[ , grepl(paste(x, collapse = "|") , colnames(beta) ) ] # subset beta on iPSC and the stage selected
    
    beta_sub = cbind(beta_sub,probe.features[rownames(beta_sub),c("CHR","MAPINFO","gene","feature","cgi")]) # get relevant info from probes
    
    
    message("---------------Selecting probes------------------")
    
    #Data for first plot (probes and DMR segment)
    plot_betas=beta_sub[which(beta_sub$CHR == region[1] & beta_sub$MAPINFO > region[2] & beta_sub$MAPINFO < region[3]) , ]
    plot_betas=melt(plot_betas,id.vars = c("CHR","MAPINFO","gene","feature","cgi") ) # melt data by sample
    plot_betas$variable=as.character(plot_betas$variable)
    plot_betas[grepl("-ISL|-R" , plot_betas$variable),"variable"] <- rep("islets",length(plot_betas[grepl("-ISL|-R" , plot_betas$variable),"variable"])) # change names to islets rows
    plot_betas$variable=gsub( "-.*$", "", plot_betas$variable) # take out everything after "-" (sample name)
    plot_betas$variable=factor(plot_betas$variable,levels=unique(plot_betas$variable))
    
    ###check meaning of columns in DMR!!!!
    #fix condition to select specific row
    message("---------------selecting DMRs------------------")
    
    selected_DMR=DMR_timecourse[[s]][which(DMR_timecourse[[s]]$BumphunterDMR.seqnames==paste("chr",region[1],sep="") & DMR_timecourse[[s]]$BumphunterDMR.start > region[2] & DMR_timecourse[[s]]$BumphunterDMR.end < region[3]),]
    if(nrow(selected_DMR)==0){
      message("---------------no DMRs found------------------")
      
    }else{
      selected_DMR$color <- ifelse(selected_DMR$BumphunterDMR.value>0,"red", "blue")
    }
    
    message("---------------Plotting probes------------------")
    
    
    scaleFUN <- function(x) sprintf("%.1f", x)  # to round y axis to 1 decimal place
    
    library(RColorBrewer)
    myColors <- brewer.pal(8,"Set2") # select colors from pallette
    myColors <- myColors[c(8,1)] # green and dark grey
    names(myColors) <- levels(plot_betas$variable) # assign levels ("stage x" and "iPSC")
    colScale <- scale_colour_manual(name = "variable",values = myColors)
    
    p1 <- ggplot(plot_betas, aes(x = MAPINFO,y=value,group=variable))  + geom_point(aes(col=variable)) + # this with colScale sets colors to variables and includes in label
      colScale 
    
    if(nrow(selected_DMR)!=0){
      for(n in 1:nrow(selected_DMR)){
        p1 <- p1 + geom_segment(show.legend = F,x=selected_DMR[n,"BumphunterDMR.start"],xend=selected_DMR[n,"BumphunterDMR.end"],y=0,yend=0,col=selected_DMR[n,"color"])  # color outside aes doesn't include it in legend
      }
    }
    
    p1 <- p1 +
      scale_y_continuous(labels=scaleFUN) +  # rounding y axis
      ggtitle(paste("DMR",unique(plot_betas$variable)[2],"vs",unique(plot_betas$variable)[1],sep=" ")) +
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
    
    
    
    message("---------------adding gene annotation------------------")
    
    ### add gene annotation
    
    library(Homo.sapiens)
    
    library(dplyr)
    library(biomaRt)
    
    ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") 
    # I use genome coordinates from GRCh37 (not GRCh38, which is the latest release). Keep that in mind.
    
    
    
    chr.region=paste(region, collapse = ":")
    filterlist <- list(chr.region)  
    
    
    
    results <- getBM(attributes = c('ensembl_gene_id','entrezgene', 'external_gene_name', 'gene_biotype', 'chromosome_name', "start_position", "end_position","transcription_start_site"),
                     filters = c("chromosomal_region"),values = filterlist, mart = ensembl)
    
    # results_refseqmRNA <- getBM(attributes=c('refseq_mrna','external_gene_name',"chromosome_name","transcript_start","transcript_end","transcription_start_site"), filters = c("chromosomal_region","biotype"), values = filterlist, mart = ensembl)
    # 
    # results_refseqmRNA[which(results_refseqmRNA$refseq_mrna==""),1] <- NA   # changing empty vectors of first column to NA, to remove them later
    # 
    # results_refseqmRNA=results_refseqmRNA[which(!is.na(results_refseqmRNA$refseq_mrna)),]  # taking out non refseq mrnas
    # 
    # tss=unique(results_refseqmRNA$transcript_start)  
    
    ####### if it throws errors of aesthetics, reduce gene size:  ##############
    
    # results=results[which(results$external_gene_name=="SNRPN"),]
    
    ##################
    
    results$external_gene_name=as.factor(results$external_gene_name)
    results_tss <- melt(results[c(3,8)])
    results_tss <- unique(results_tss)
    results_tss_list <- split(results_tss, results_tss$external_gene_name) # make a list for each gene
    
    message("---------------plotting genes------------------")
    
    # plot white y axis
    # plot white dots
    # plot geom segment as results or as region,  depending on whether it starts or ends outside plotting region
    
    
    myColors <- brewer.pal(8,"Dark2") # select colors from pallette
    names(myColors) <- levels(results$external_gene_name) # assign levels ("stage x" and "iPSC")
    colScale <- scale_colour_manual(name = "external_gene_name",values = myColors)
    
    p2 <-   ggplot(plot_betas, aes(x = MAPINFO,y=value)) + geom_point(col="white") + colScale
    
    
    values <- length(unique(results$external_gene_name))   # how many genes do I have
    n=1  # count for the location of the lines
    for (i in unique(results$external_gene_name)) {  # sequentially plot new genes. Doesn't work for more than 8
      p2 <- p2 + geom_segment(data=results[which(results$external_gene_name==i),], 
                              show.legend = T,
                              x=unique(results[which(results$external_gene_name==i),6]),
                              xend=unique(results[which(results$external_gene_name==i),7]), 
                              y=(mean(plot_betas$value)*sin(n/4))+0.1,yend=(mean(plot_betas$value)*sin(n/4))+0.1,  # oscilate around 0.5
                              aes(col=external_gene_name))
      print(n)
      print(i)
      
      for (t in results_tss_list[[i]]$value) {  # plot all tss for each gene along the same line
        p2 <- p2 +geom_segment(show.legend = F,x=t,xend=t,
                               y=(mean(plot_betas$value)*sin(n/4))+0.1-0.08,yend=(mean(plot_betas$value)*sin(n/4))+0.1+0.08)  # making a tick
        
      }
      n=n+1
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
    
    message("---------------plotting delta beta values------------------")
    
    #Data for deltabeta
    if(s=="EN7"){
      x = c(s,"iPSC") # two stages to contrast
      plot_deltabeta=deltaBetas[ ,grepl(paste(x, collapse = "-"), colnames(deltaBetas))] # subset beta on iPSC and the stage selected
    }
    if(s=="islets"){
      x = c(s,"iPSC") 
      plot_deltabeta=deltaBetas[ ,grepl(paste(x, collapse = "-"), colnames(deltaBetas))] 
      
    }
    if(s=="islets-EN7"){
      plot_deltabeta=deltaBetas[ ,grepl(s, colnames(deltaBetas))] 
      
    }
    
    
    plot_deltabeta = cbind(plot_deltabeta ,beta_sub[,(ncol(beta_sub)-4):ncol(beta_sub)]) # join location of probes and deltabeta for selected stage
    plot_deltabeta=plot_deltabeta[which(plot_deltabeta$CHR == region[1] & plot_deltabeta$MAPINFO > region[2] & plot_deltabeta$MAPINFO < region[3]) ,]
    colnames(plot_deltabeta)[1]="deltaBeta"  
    
    limits_deltabeta=c(-1.0,1.0)
 
  
    
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
    
    message("---------------arranging all plots------------------")
    
    p4 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),ncol=1),
                       arrangeGrob(p2 + theme(legend.position="none"),ncol=1),
                       arrangeGrob(p3 + theme(legend.position="none"),ncol=1),
                       arrangeGrob(mylegend,mylegend2,ncol=2),
                       nrow=4,heights=c(7,size_p2,4.5,1)) 
    
    # blank grob? blank<-rectGrob(gp=gpar(col="white"))
    
    message("---------------saving all plots------------------")
    
    
    ggsave(plot = p4,
           device="png",
           filename = paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/Methylation_EPIC/DMR_and_DMP_timecourse_", paste(x, collapse = "-"),
                            paste(unique(results$external_gene_name), collapse = "_"),
                            paste(region, collapse = "_"),
                            currentDate,".png",sep="_"),
           width = 18, height = 6, units = "in",
           dpi = 400)
    
    
    message("---------------saving results in tables------------------")
    
    # table with DMPs
    
    selected_DMP=DMP_timecourse[[s]][which(DMP_timecourse[[s]]$CHR==region[1] & DMP_timecourse[[s]]$MAPINFO > region[2] & DMP_timecourse[[s]]$MAPINFO < region[3]),]
    capture.output(selected_DMP, 
                   file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMP_timecourse_",paste(x, collapse = "-"),
                              paste(unique(results$external_gene_name), collapse = "_"),
                              paste(region, collapse = "_"),
                              currentDate,".txt",sep="_")) # to be able to check SNPs
    
    # taking out SNPs column, that makes the .csv file unreadable
    
    selected_DMP2=selected_DMP[,c(1:15)]
    
    write.csv(selected_DMP2,
              paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMP_timecourse_",paste(x, collapse = "-"),
                    paste(unique(results$external_gene_name), collapse = "_"),
                    paste(region, collapse = "_"),
                    currentDate,".csv",sep="_"),
              row.names=T, quote=F)
    
    write.csv(selected_DMP,
              paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMP_timecourse_iPSC_vs",paste(x, collapse = "-"),
                    paste(unique(results$external_gene_name), collapse = "_"),
                    paste(region, collapse = "_"),
                    currentDate,"plusSNPs.csv",sep="_"),
              row.names=T, quote=F)
  }
}
# calculate mean beta value (only works for one gene) for selected stage (not iPSC)

mean(plot_betas[which(plot_betas$MAPINFO<unique(results$end_position) & plot_betas$MAPINFO>unique(results$start_position) & plot_betas$variable!="iPSC"),7])

# data from adult islets:

islets2=colnames(beta)[22:32]

x_islets = islets2 # stages to subset
beta_sub_islets=beta[ , grepl(paste(x_islets, collapse = "|") , colnames(beta) ) ] # subset beta on iPSC and the stage selected

beta_sub_islets = cbind(beta_sub_islets,probe.features[rownames(beta_sub_islets),c("CHR","MAPINFO","gene","feature","cgi")]) # get relevant info from probes

gene="SNRPN"
region <- c(unique(results[which(results$external_gene_name==gene),5]),
            unique(results[which(results$external_gene_name==gene),6]),
            unique(results[which(results$external_gene_name==gene),7])) # just one gene

plot_betas_islets=beta_sub_islets[which(beta_sub_islets$CHR == region[1] & beta_sub_islets$MAPINFO > region[2] & beta_sub_islets$MAPINFO < region[3]) , ]

plot_betas_islets=melt(plot_betas_islets,id.vars = c("CHR","MAPINFO","gene","feature","cgi") )
plot_betas_islets$variable=factor(plot_betas_islets$variable,levels=islets)


# calculate mean beta value (only works for one gene) for adult islets

mean(plot_betas_islets[which(plot_betas_islets$MAPINFO<unique(results$end_position) & plot_betas_islets$MAPINFO>unique(results$start_position) & plot_betas_islets$feature=="5'UTR"),7])


# for all stages (not islets)


beta_sub_stages=beta[ , grepl(paste(stages, collapse = "|") , colnames(beta) ) ] # subset beta on iPSC and the stage selected

beta_sub_stages = cbind(beta_sub_stages,probe.features[rownames(beta_sub_stages),c("CHR","MAPINFO","gene","feature","cgi")]) # get relevant info from probes



plots_betas_stages=beta_sub_stages[which(beta_sub_stages$CHR == region[1] & beta_sub_stages$MAPINFO > region[2] & beta_sub_stages$MAPINFO < region[3]) , ]

plots_betas_stages=melt(plots_betas_stages,id.vars = c("CHR","MAPINFO","gene","feature","cgi") )
plots_betas_stages$variable=gsub( "-.*$", "", plots_betas_stages$variable) # take out everything after "-" (sample name)

plots_betas_stages$variable=factor(plots_betas_stages$variable,levels=stages)


# calculate mean beta value (only works for one gene) for all stages


mean(plots_betas_stages[which(plots_betas_stages$gene==gene  & (plots_betas_stages$feature=="TSS1500" | plots_betas_stages$feature=="TSS200")),7])

# just EN7, just TSS
mean(plots_betas_stages[which(plots_betas_stages$gene==gene  & (plots_betas_stages$feature=="TSS1500" | plots_betas_stages$feature=="TSS200") & plots_betas_stages$variable=="EN7"),7])

# just EN7, all gene
mean(plots_betas_stages[which(plots_betas_stages$gene==gene   & plots_betas_stages$variable=="EN7"),7])

# adult islets, just TSS
mean(plot_betas_islets[which(plot_betas_islets$gene==gene  & (plot_betas_islets$feature=="TSS1500" | plot_betas_islets$feature=="TSS200")),7])

#specific probes

rowMeans(beta["cg02490034",c(20:21)])  # BLC
rowMeans(beta["cg02490034",c(22:32)]) # islets
