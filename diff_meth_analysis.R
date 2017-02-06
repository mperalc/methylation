# differential methylation analysis (DMA)

# check what I need and save in QC script

library(ChAMP)
currentDate <- Sys.Date() # to save date in name of output files

# load matrix of normalized beta values
beta= read.csv( "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/quantile_normalised_beta_detP_0.01_nocrossreact.csv", header=T,row.names = 1,check.names=F)

# load matrix with sample information
pD=read.csv("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/samples_info.csv", header=T, row.names=1, check.names = F)

stages=c("iPSC","DE","PGT","PFG","PE","EP","EN6","EN7")

islets=colnames(beta)[22:32]  # islets names
diff_type="timecourse"   # type of differential analysis: peak or timecourse
sample_type="diff"   # islets or diff
plots=TRUE
#############   function to perform differential methylation analysis

diff_meth_ChAMP=function(beta,pD,diff_type,sample_type,plots=FALSE){

# select type of samples I want to test:
  
  if(sample_type=="islets"){
    
    beta_diffcells=beta[,22:32]
    
  }
  if(sample_type=="diff"){
    beta_diffcells=beta[,1:21]
  }

  beta_diffcells=as.matrix(beta_diffcells)
  # create the design matrix, which contains the comparisons I want to make
  
  ####### create different matrices for each stage
  
  if(diff_type=="peak"){
    # design_peak
    group_peak_iPSC=c(rep(1,times=3),rep(-1,times=18))  # iPSC vs all others. Change for each peak test
    design_peak=pD[1:21,c(3,4)]
    design_peak$group=group_peak
  }
  if(diff_type=="timecourse"){
    if(sample_type=="diff"){
      group_timecourse <- list()
      group_timecourse[["DE"]]=c(rep(-1,times=3),rep(1,times=3),rep(0,times=15)) #DE vs iPSC
      group_timecourse[["PGT"]]=c(rep(-1,times=3),rep(0,times=3),rep(1,times=3),rep(0,times=12)) #PGT vs iPSC
      group_timecourse[["PFG"]]=c(rep(-1,times=3),rep(0,times=6),rep(1,times=3),rep(0,times=9)) #PFG vs iPSC
      group_timecourse[["PE"]]=c(rep(-1,times=3),rep(0,times=9),rep(1,times=3),rep(0,times=6)) #PE vs iPSC
      group_timecourse[["EP"]]=c(rep(-1,times=3),rep(0,times=12),rep(1,times=1),rep(0,times=5)) #EP vs iPSC. Just one sample, so not really diff meth
      group_timecourse[["EN6"]]=c(rep(-1,times=3),rep(0,times=13),rep(1,times=3),rep(0,times=2)) #EN vs iPSC
      group_timecourse[["EN7"]]=c(rep(-1,times=3),rep(0,times=16),rep(1,times=2)) #BLC vs iPSC
      # iPSC overmethylated regions are the undermethylated regions of the above 7
    }
    # write for islets
  }
  
  if(sample_type=="islets"){
    design= pD[22:32,c(3,4)]
  }else{
    design= pD[1:21,c(3,4)]
  }

  if(diff_type=="timecourse"){
    if(sample_type=="diff"){
      design <- rep(list(design),7)
      names(design_timecourse) = stages[2:length(stages)]
      for(s in stages[2:length(stages)]){
        design_timecourse[[s]]$group=group_timecourse[[s]] # fil list of df with appropriate design
        design_timecourse[[s]]=design_timecourse[[s]][design_timecourse[[s]]$group!=0,] # select only rows with -1 or 1. Check this is necessary
      }
      if(sample_type=="islets"){
        # write
      }
      
    }
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
    write.csv(myDMR_timecourse[[s]],paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMR/",s,"_timecourse_DMR_CpGs_within_900bp_",currentDate,".csv",sep=""), col.names=T,row.names=T, quote=F)
  }
  
  DMR.GUI(DMR=myDMR_timecourse[["EN7"]],beta=beta_sub,pheno=design_timecourse$EN7$group,arraytype = "EPIC",runDMP = F) 
  # not doing what it's supposed to do
  
  
  ######################### Differentially methylated probes#########################
  
  
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
      myDMP_timecourse[[s]] <- DMP 
      message("[<<<<<< ChAMP.DMP VARIATION ENDED, JUST AS ALL THINGS END IN LIFE >>>>>>]")
      message("[===========================]")
      
    }else{
      
      myDMP_timecourse[[s]] <- champ.DMP(beta = beta_sub, 
                                         pheno=factor(design_timecourse[[s]]$group,levels=c("iPSC",s)), #levels argument necesary, otherwise comparisons will be determined by alphabetical order
                                         arraytype = "EPIC"   ) # call function
    }
    
  }
  
  myDMP_timecourse_beta_forplots=myDMP_timecourse  # saving for plotting delta beta etc.
  
  
  
  
  ################ now with M values
  
  
  for(s in stages[2:length(stages)]){
    x = c("iPSC",s) # two stages to contrast
    print(paste(" Testing contrast", paste(x, collapse = "|") ,sep=" ")) # message of progress
    beta_sub=beta_diffcells[ , grepl(paste(x, collapse = "|") , colnames( beta_diffcells ) ) ] # subset beta on contrast stages
    beta_sub[beta_sub<=0.001]<-0.001 # replacing extreme values
    beta_sub[beta_sub>=0.999]<-0.999
    
    M <- log((beta_sub/(1-beta_sub)),2) # calculates M values
    
    M=as.matrix(M) # Y needs to be a matrix
    
    if(is.null(dim(beta_diffcells[ , grepl(x[2] , colnames( beta_diffcells ))]))){
      # if second stage (not iPSC) only has one column, do this alternative version of champ.DMP
      # Here there's no average beta value for each CpG in the stage with one sample, just its only beta value
      # As there are not replicates, this is not a true differential methylation analysis. Hence the message:
      message(paste("The stage",x[2],"has only one sample. Doing DMP variation, but do not trust results.",sep=" "))
      message("[===========================]")
      message("[<<<<< ChAMP.DMP VARIATION STARTING >>>>>]")
      message("-----------------------------")
      
      ### setup
      beta=M
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
      myDMP_timecourse[[s]] <- DMP 
      message("[<<<<<< ChAMP.DMP VARIATION ENDED, JUST AS ALL THINGS END IN LIFE >>>>>>]")
      message("[===========================]")
      
    }else{
      
      myDMP_timecourse[[s]] <- champ.DMP(beta = M, 
                                         pheno=factor(design_timecourse[[s]]$group,levels=c("iPSC",s)), #levels argument necesary, otherwise comparisons will be determined by alphabetical order
                                         arraytype = "EPIC"   ) # call function
    }
    
  }
  
  # calculate deltaBeta (MethDiff in Jaffe's paper - 2016) for all contrasts and probes (for future plots)
  avg <- data.frame(matrix(nrow = nrow(beta_diffcells),ncol = 0))
  
  for(s in stages[1:length(stages)]){
    
    #get average beta for all samples in each stage
    if(is.null(dim(beta_diffcells[,grepl(s, colnames( beta_diffcells ))]))){ #for stages with just one sample
      avg[paste(s,"AVG",sep="-") ] <-  beta_diffcells[,grepl(s,colnames(beta_diffcells))] # not really average, because there's just one value per CpG
      
    }else{
      avg[paste(s,"AVG",sep="-") ] <-  rowMeans(beta_diffcells[,grepl(s,colnames(beta_diffcells))])
      
    } 
    if(s != "iPSC"){
      
      #if not iPSC stage, substract means to get deltaBeta for that second stage
      avg <- cbind(avg,avg[,paste(s,"AVG",sep="-") ]-avg[,paste("iPSC","AVG",sep="-") ])
      colnames(avg)[length(avg)] <- paste(s,"deltaBeta",sep="_")
    }
  }
  
  rownames(avg) <- rownames(beta_diffcells) # add names of CpGs
  
  write.csv(avg,"/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/deltaBeta_allprobes.csv", col.names=T,row.names=T, quote=F)
  #save
  
  
  
  # merge AVG and deltaBeta values with the results of diff meth probes using M values
  
  myDMP_timecourse_merged <- lapply(myDMP_timecourse, "[", c(1:6,10:19)) # subset dataframes without AVG or deltaBeta
  myDMP_timecourse_beta_forplots=lapply(myDMP_timecourse_beta_forplots, "[", c(7:9))  # subset with AVG and deltaBeta
  
  
  for(s in stages[2:length(stages)]){
    
    myDMP_timecourse_merged[[s]] <- cbind(myDMP_timecourse_merged[[s]], myDMP_timecourse_beta_forplots[[s]][row.names(myDMP_timecourse_merged[[s]]),]) #merge by probes from M value df
    
  }
  
  rm(myDMP_timecourse_beta_forplots,myDMP_timecourse)
  
  #save with all probes
  for(s in stages[2:length(stages)]){
    write.csv(myDMP_timecourse_merged[[s]],paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMP/",s,"_timecourse_DMPs_",currentDate,".csv",sep=""), col.names=T,row.names=T, quote=F)
  }
  
  
  
  
  
  
  
  
  
  
}
################ get max logFC o all stages for each CpG: that will be its maximum value ##################
###not sure if this is very relevant, as there is a threshold of methylation/not methylation


### FIX THE FOLLOWING FOR MY DATA#####
#take logFC and adj P values of each, and combine in single dataframe

combined_df <- lapply(DE_list, "[", c(5,9))   # subsetting list of dataframes with columns I want

combined_df <- do.call("cbind", combined_df)   # merging into one dataframe

combined_df <- cbind(DE_list$stagesiPSC[,c(1,2,3,4)],combined_df)

#get max conditions per gene
maxVals <- apply(combined_df,1,function(x) which.max(x[c(5,7,9,11,13,15,17,19)]))

#find significant genes per stage
#rowsums checks that there's at least one positive logFC 
# maxvals selects which one to take as max value for each stage
DE_stages<-list()  

sig_iPSC_stage <- combined_df[combined_df$stagesiPSC.adj.P.Val < 0.01 & maxVals == 1,c(1:6)]
sig_iPSC_stage=sig_iPSC_stage[order(sig_iPSC_stage[6]),] #order by adj p values
rownames(sig_iPSC_stage)=NULL # take out row names
DE_stages [["iPSC"]]<-sig_iPSC_stage

sig_DE_stage <- combined_df[combined_df$stagesDE.adj.P.Val < 0.01 &  maxVals == 2, c(1:4,7,8) ]
sig_DE_stage=sig_DE_stage[order(sig_DE_stage[6]),]
rownames(sig_DE_stage)=NULL 
DE_stages [["DE"]]<-sig_DE_stage

sig_PGT_stage <- combined_df[ combined_df$stagesPGT.adj.P.Val < 0.01 & maxVals == 3,c(1:4,9,10)  ]
sig_PGT_stage=sig_PGT_stage[order(sig_PGT_stage[6]),]
rownames(sig_PGT_stage)=NULL 
DE_stages [["PGT"]]<-sig_PGT_stage

sig_PFG_stage <- combined_df[ combined_df$stagesPFG.adj.P.Val < 0.01 & maxVals == 4, c(1:4,11,12) ]
sig_PFG_stage=sig_PFG_stage[order(sig_PFG_stage[6]),]
rownames(sig_PFG_stage)=NULL 
DE_stages [["PFG"]]<-sig_PFG_stage

sig_PE_stage <- combined_df[ combined_df$stagesPE.adj.P.Val < 0.01 & maxVals == 5,c(1:4,13,14)  ]
sig_PE_stage=sig_PE_stage[order(sig_PE_stage[6]),]
rownames(sig_PE_stage)=NULL 
DE_stages [["PE"]]<-sig_PE_stage

sig_EP_stage <- combined_df[ combined_df$stagesEP.adj.P.Val < 0.01 & maxVals == 6, c(1:4,15,16)  ]
sig_EP_stage=sig_EP_stage[order(sig_EP_stage[6]),]
rownames(sig_EP_stage)=NULL 
DE_stages [["EP"]]<-sig_EP_stage

sig_EN6_stage <- combined_df[ combined_df$stagesEN6.adj.P.Val < 0.01 & maxVals == 7, c(1:4,17,18) ]
sig_EN6_stage=sig_EN6_stage[order(sig_EN6_stage[6]),]
rownames(sig_EN6_stage)=NULL 
DE_stages [["EN6"]]<-sig_EN6_stage

sig_EN7_stage <- combined_df[ combined_df$stagesEN7.adj.P.Val < 0.01 & maxVals == 8,c(1:4,19,20)  ]
sig_EN7_stage=sig_EN7_stage[order(sig_EN7_stage[6]),]
rownames(sig_EN7_stage)=NULL 
DE_stages [["EN7"]]<-sig_EN7_stage

