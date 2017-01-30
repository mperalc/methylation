# plot DMPs and DMRs


################## setup ###################

library(ggplot2)
library(readr)
library(reshape2)
library(grid)
library(gridExtra)
library(gtable)

############### region to plot ##########33

region <- c(7, 5380000,5392500)  # Chr, position start, position finish
timecourse_stage = "DE"  # stage from timecourse to plot (meth values compared to iPSC stage)

######### read in files ########

# deltaBeta values of all probes (stage (n) - iPSC), and beta values averages for every stage
deltaBetas <- as.data.frame(read_csv(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/deltaBeta_allprobes.csv"))
rownames(deltaBetas) <- deltaBetas$X1
deltaBetas = deltaBetas[2:ncol(deltaBetas)]

# beta values for all probes

beta= read.csv( "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/quantile_normalised_beta_detP_0.01_nocrossreact.csv", header=T,row.names = 1,check.names=F)
x = c("iPSC",timecourse_stage) # two stages to subset
beta_sub=beta[ , grepl(paste(x, collapse = "|") , colnames(beta) ) ] # subset beta on iPSC and the stage selected

## probe features
library(ChAMPdata)
data(probe.features.epic)
beta_sub = cbind(beta_sub,probe.features[rownames(beta_sub),c("CHR","MAPINFO","gene","feature","cgi")]) # get relevant info from probes


stages=c("iPSC","DE","PGT","PFG","PE","EP","EN6","EN7")

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

#################

#Data for first plot (probes and DMR segment)
plot_betas=beta_sub[which(beta_sub$CHR == region[1] & beta_sub$MAPINFO > region[2] & beta_sub$MAPINFO < region[3]) , ]
plot_betas=melt(plot_betas,id.vars = c("CHR","MAPINFO","gene","feature","cgi") )
plot_betas$variable=gsub( "-.*$", "", plot_betas$variable) # take out everything after "-" (sample name)

###check meaning of columns in DMR!!!!
#fix condition to select specific row
if(DMR_timecourse[[timecourse_stage]]$BumphunterDMR.value>0){  # selecting colour of segment based on over/under methylation. CHANGE!!
  color_DMR <- "red" #over-methylated in red
}else{ #under-methylated in blue
  color_DMR <- "blue"
  
}

scaleFUN <- function(x) sprintf("%.1f", x)  # to round y axis to 1 decimal place

plot_betas$col <- ifelse(plot_betas$variable=="iPSC", "black", "green") # add color column conditionally

p1 <- ggplot(plot_betas, aes(x = MAPINFO,y=value,group=variable)) +
  geom_segment(show.legend = F,aes(x=5388781,xend=5391498,y=0,yend=0)) + #,col=color_DMR) + #fix coordinates +
  geom_point(stat="identity",aes(colour=col)) +
  scale_colour_identity() +
  scale_y_continuous(labels=scaleFUN) +
  
  # scale_color_manual(values =factor(c("#E69F00","#999999",
  #                              "#FF0000"),levels=c("#E69F00","#999999",
  #                                                  "#FF0000")),labels=factor(c("DE","iPSC",color_DMR),levels=c("DE","iPSC",color_DMR))) +
  #                    #           ifelse(color_DMR=="over-methylated","#F43431","#2CD0FC")),
  #                    # labels=c("DE","iPSC",color_DMR))+
  ggtitle("DMR") +
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
# 
 p1 <- p1 + theme(legend.position="none")

####wrong colors!!!

#Data for second plot (deltabeta)
plot_deltabeta=deltaBetas[ ,grepl(timecourse_stage, colnames(deltaBetas))] # subset beta on iPSC and the stage selected

plot_deltabeta = cbind(plot_deltabeta ,beta_sub[7:11]) # join location of probes and deltabeta for selected stage
plot_deltabeta=plot_deltabeta[which(plot_deltabeta$CHR == region[1] & plot_deltabeta$MAPINFO > region[2] & plot_deltabeta$MAPINFO < region[3]) ,]
plot_deltabeta=plot_deltabeta[,c(2:ncol(plot_deltabeta))]
colnames(plot_deltabeta)[1]="deltaBeta"  

#dynamic limits:
#if max deltabeta of the region >0.1, pick max=(max + 0.1) and min -0.15
# if min deltabeta <-0.1, pick min=(min -0.1) and max 0.15
limits_deltabeta=c(-0.2,0.2)
if(max(plot_deltabeta[1])>=0.1){
  limits_deltabeta[2]=max(plot_deltabeta[1]) + 0.1
  }
if(min(plot_deltabeta[1])<=-0.1){
  limits_deltabeta[1]=min(plot_deltabeta[1]) - 0.1
  }




p2 <- ggplot(plot_deltabeta, aes(x = MAPINFO,y=deltaBeta)) + geom_line() + theme_minimal() + 
  xlab (paste("Location in chr",region[1],sep=" ")) +
  ylab ("delta Beta") +

  scale_y_continuous(labels=scaleFUN) +
  expand_limits(y = limits_deltabeta) +
  geom_hline(yintercept=0,linetype="dashed",size=0.1) +
  geom_hline(yintercept=0.1,col="red",size=0.1) +
  geom_hline(yintercept=0.1,col="red",size=0.1) +
  
  theme(axis.text=element_text(size=12,face="bold"),
  axis.title=element_text(size=14,face="bold"))

p3 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),ncol=1),
                   arrangeGrob(p2 + theme(legend.position="none"),ncol=1),
                   arrangeGrob(mylegend,ncol=1), nrow=3,heights=c(7,3,1))   




# using grid.draw
grid.newpage()
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2),size="first"))    # this works



# other version using gtable:

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)


## merge, basing widths on g1
g <- gtable:::rbind_gtable(g1, g2, "first")
g <- gtable_add_rows(g, unit(-0.5,"cm"), pos=nrow(g1)) # no spacing
grid.newpage()
grid.draw(g)




