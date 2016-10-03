# minfi tutorial

library(minfi)
library(minfiData)
library(sva)

# The starting point of minfi is reading the .IDAT files with the built-in function read.450k.exp. 

# we will load the dataset containing 6 samples from the minfiData package using the sample sheet provided within the package:
  
  baseDir <- system.file("extdata", package="minfiData")
  targets <- read.450k.sheet(baseDir)
  
  RGSet <- read.450k.exp(targets = targets)

  # The class of RGSet is a RGChannelSet object. This is the initial object of a minfi analysis that contains the raw intensities 
  # in the green and red channels. Note that this object contains the intensities of the internal control probes as well. Because 
  # we read the data from a data sheet experiment, the phenotype data is also stored in the RGChannelSet and can be accessed via the 
  # accessor command pData:
  
  phenoData <- pData(RGSet)
  phenoData[,1:6]
  
  # The RGChannelSet stores also a manifest object that contains the probe design information of the array:
  
  manifest <- getManifest(RGSet)
  manifest
  
  head(getProbeInfo(manifest))
  
 #  A MethylSet objects contains only the methylated and unmethylated signals. You create this by
  
  MSet <- preprocessRaw(RGSet) 
  MSet
  
 # This function matches up the different probes and color channels. Note that the dimension of this object is much smaller than for 
  # the RGChannelSet; this is because CpGs measured by type I probes are measured by 2 probes.
  
  # The accessors getMeth and getUnmeth can be used to get the methylated and unmethylated intensities matrices:
    
    head(getMeth(MSet)[,1:3])
    
    head(getUnmeth(MSet)[,1:3])
    
  #  A RatioSet object is a class designed to store Beta values and/or M values instead of the methylated and unmethylated signals. 
    # An optional copy number matrix, CN, the sum of the methylated and unmethylated signals, can be also stored. Mapping a MethylSet 
    # to a RatioSet may be irreversible, i.e. one cannot be guranteed to retrieve the methylated and unmethylated signals from a RatioSet. 
    # A RatioSet can be created with the function ratioConvert:
      
    RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
    RSet
    
   # The functions getBeta, getM and getCN return respectively the Beta value matrix, M value matrix and the Copy Number matrix.
    
    beta <- getBeta(RSet)
    
    # The function mapToGenome applied to a RatioSet object will add genomic coordinates to each probe together with some additional annotation 
    # information. The output object is a GenomicRatioSet (class holding M or/and Beta values together with associated genomic coordinates). 
    # It is possible to merge the manifest object with the genomic locations by setting the option mergeManifest to TRUE.
    
    GRset <- mapToGenome(RSet)
    GRset
    
    # Note that the GenomicRatioSet extends the class SummarizedExperiment. Here are the main accessors functions to access the data:
      
      beta <- getBeta(GRset)
    M <- getM(GRset)
    CN <- getCN(GRset)
    
    sampleNames <- sampleNames(GRset)
    probeNames <- featureNames(GRset)
    pheno <- pData(GRset)
    
    # To return the probe locations as a GenomicRanges objects, one can use the accessor granges:
      
      gr <- granges(GRset)
    head(gr, n= 3)
    
    # To access the full annotation, one can use the command getAnnotation:
      # this function might not work with 850k. check!
      
      annotation <- getAnnotation(GRset)
    names(annotation)
    
    
    
    ############### QC ##################################
    
    
    
    # minfi provides a simple quality control plot that uses the log median intensity in both the methylated (M) and unmethylated (U) channels. 
    # When plotting these two medians against each other, it has been observed that good samples cluster together, while failed samples tend to 
    # separate and have lower median intensities. In order to obtain the methylated and unmethylated signals, we need to convert the RGChannelSet 
    # to an object containing the methylated and unmethylated signals using the function preprocessRaw. It takes as input a RGChannelSet and converts 
    # the red and green intensities to methylated and unmethylated signals according to the special 450K probe design, and returns the converted signals
    # in a new object of class MethylSet. It does not perform any normalization.
    
    # The accessors getMeth and getUnmeth can be used to get the methylated and unmethylated intensities matrices:
      
      head(getMeth(MSet)[,1:3])
    head(getUnmeth(MSet)[,1:3])
    
    # The functions getQC and plotQC are designed to extract and plot the quality control information from the MethylSet:
      
      qc <- getQC(MSet)
    head(qc)
    
    plotQC(qc)
    
    # Moreover, the function addQC applied to the MethylSet will add the QC information to the phenotype data.
    
    # To further explore the quality of the samples, it is useful to look at the Beta value densities of the samples, with the option to color 
    # the densities by group:
      
      densityPlot(MSet, sampGroups = phenoData$Sample_Group)
      
    # or density bean plots:
      
      densityBeanPlot(MSet, sampGroups = phenoData$Sample_Group)
      
      # shinyMethyl is particularly useful to visualize all plots at the same time in an interactive fashion.
      
      
      
      
########### # Control probes plot
      
     # The 450k array contains several internal control probes that can be used to assess the quality control of different sample preparation 
      # steps (bisulfite conversion, hybridization, etc.). The values of these control probes are stored in the initial RGChannelSet and can be 
      # plotted by using the function controlStripPlot and by specifying the control probe type:
        
        controlStripPlot(RGSet, controls="BISULFITE CONVERSION II")
        
    # All the plots above can be exported into a pdf file in one step using the function qcReport:
          
          qcReport(RGSet, pdf= "qcReport.pdf")   # it's not plotting properly
          
  ####### sex prediction
          
# By looking at the median total intensity of the X chromosome-mapped probes, denoted med(X), and the median total intensity of the 
# Y-chromosome-mapped probes, denoted med(Y), one can observe two different clusters of points corresponding to which gender the samples belong to.
          
          predictedSex <- getSex(GRset, cutoff = -2)$predictedSex
          head(predictedSex)     
          
          # To choose the cutoff to separate the two gender clusters, one can plot med(Y) against med(Y)
          # 
          # with the function plotSex:
            
            plotSex(getSex(GRset, cutoff = -2))
            
   # Finally, the function addSex applied to the GenomicRatioSet will add the predicted sex to the phenotype data.
           
  # Note: the function does not handle datasets with only females or only males
            
            
  ###### END of QC
    