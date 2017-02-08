### annotate DMR results


library(readr)
library(reshape2)
library(digest)
library(Homo.sapiens)
library(biomaRt)
library(GenomicRanges)

splitByOverlap <-  function(query, subject, column="ENTREZID", ...)
  {
    olaps <- findOverlaps(query, subject, ...)
    f1 <- factor(subjectHits(olaps),
                 levels=seq_len(subjectLength(olaps)))
    splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)
  }

geneRanges <- function(db, column="ENSEMBL")
{
  g <- genes(db, columns=column)
  col <- mcols(g)[[column]]
  genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
  mcols(genes)[[column]] <- as.character(unlist(col))
  genes
}

currentDate <- Sys.Date() # to save date in name of output files
diff_type="timecourse"   # type of differential analysis: peak or timecourse
#sample_type="islets"   # islets or differentiated cells ("diff")
stages=c("iPSC","DE","PGT","PFG","PE","EP","EN6","EN7")
islets=c("EN7","islets","islets-EN7")

types=c("diff","islets") # to plot both



for(sample_type in types){
  if(diff_type=="timecourse"){
    if(sample_type=="diff"){
      
      DMR_timecourse <- list()
      for(s in stages[2:length(stages)]){
        DMR_timecourse[[s]]<- as.data.frame(read_csv(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMR/",s,"_timecourse_DMR_CpGs_within_900bp.csv",sep="")))
        rownames(DMR_timecourse[[s]]) <- DMR_timecourse[[s]]$X1
        DMR_timecourse[[s]] = DMR_timecourse[[s]][2:ncol(DMR_timecourse[[s]])]
        
        # annotate DMR with genes in genomic ranges
        
        dmr = GRanges(DMR_timecourse[[s]][,1], IRanges(DMR_timecourse[[s]][,2],DMR_timecourse[[s]][,3])) # create IRanges object
        gns = geneRanges(Homo.sapiens, column="SYMBOL") # get gene ranges for all genes in human
      
        symInDMR = splitByOverlap(gns, dmr, "SYMBOL") # gets lists with genes in same order as dmr ranges  
        DMR_timecourse[[s]]$gene_name=unstrsplit(symInDMR, sep=",") # add new column to my query data frame with genes
        write.csv(DMR_timecourse[[s]],paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMR/",s,"_timecourse_DMR_CpGs_within_900bp_genenames_",currentDate,".csv",sep=""), col.names=T,row.names=T, quote=F)
        
      }
    }
    if(sample_type=="islets"){
  
      DMR_timecourse <- list()
      results <- list()
      
      for(s in islets){
        DMR_timecourse[[s]]<- as.data.frame(read_csv(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMR/",s,"_timecourse_DMR_CpGs_within_900bp_2017-02-07.csv",sep="")))
        rownames(DMR_timecourse[[s]]) <- DMR_timecourse[[s]]$X1
        DMR_timecourse[[s]] = DMR_timecourse[[s]][2:ncol(DMR_timecourse[[s]])]
        
        # annotate DMR with genes in genomic ranges
        
        dmr = GRanges(DMR_timecourse[[s]][,1], IRanges(DMR_timecourse[[s]][,2],DMR_timecourse[[s]][,3])) # create IRanges object
        gns = geneRanges(Homo.sapiens, column="SYMBOL") # get gene ranges for all genes in human
        symInDMR = splitByOverlap(gns, dmr, "SYMBOL") # gets lists with genes in same order as dmr ranges  
        DMR_timecourse[[s]]$gene_name=unstrsplit(symInDMR, sep=",") # add new column to my query data frame with genes
        write.csv(DMR_timecourse[[s]],paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Methylation/DMR/",s,"_timecourse_DMR_CpGs_within_900bp_genenames_",currentDate,".csv",sep=""), col.names=T,row.names=T, quote=F)
        
      }
    }
    
    # WRITE FOR PEAK
  }
}
  