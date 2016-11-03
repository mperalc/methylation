# Methylation QC pipeline (quick)



library(minfi)
library(ENmix)
library(GEOquery)

# Read IDAT files
#Currently minfi does not support reading compressed IDAT files. 

# using read.450k.exp() which (in this case) reads all the IDAT files in a directory.
# is deprecated, using read.metharray.exp instead



head(sampleNames(rgSet))


# we see the samples are named following a standard IDAT naming convention with a 10 digit 
# number which is an array identifier followed by an identifier of the form R01C01.

# The 200526570053_R01C01 means row 1 and column 1 on chip 200526570053. This is all good, 
# but does not help us understand which samples are cases and which are controls.


# We now get the standard GEO representation to get the phenotype data stored in GEO. 
# Most of the columns in this phenotype data are irrelevant (contains data such as the 
# address of the person who submitted the data); we keep the useful ones.

# Then we clean it.

pheno=read.csv2(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/Methylation/P160281_MethylationEPIC_NicolaBeer/03.archive/3.Results/Phenotype.csv")

pD <- pheno[,c(1,10)]

# Ask Nicola about this and fix
# add stage column

stage= c("iPSC","iPSC","iPSC","DE","DE","DE","PGT","PGT","PGT","PFG","PFG","PFG","PE","PE","PE","EP","EN6","EN6","EN6","EN7","EN7","124I","141 (14-19)","177","179","182","184","124R","136","137","138","139")
# add sample column
sample=c("neo1.1","SBAd2.1","SBAd3.1","neo1.1","SBAd2.1","SBAd3.1","neo1.1","SBAd2.1","SBAd3.1","neo1.1","SBAd2.1","SBAd3.1","neo1.1","SBAd2.1","SBAd3.1","SBAd3.1","neo1.1","SBAd2.1","SBAd3.1","neo1.1","SBAd3.1","ISL","ISL","ISL","ISL","ISL","ISL","R","R","R","R","R")

pD=cbind(pD,stage,sample)


pD$Sample_Name=gsub( "^.*?_", "", pD$Sample_Name)
rownames(pD) <- pD$Sample_Name
pD <- pD[sampleNames(rgSet),] # reordering before merging


pData(rgSet) <- pD  # merging into pheno data of rgSet

rgSet


# The RGChannelSet stores also a manifest object that contains the probe design information of the array:

manifest <- getManifest(rgSet)
manifest

head(getProbeInfo(manifest))

#  A MethylSet objects contains only the methylated and unmethylated signals. You create this by

MSet <- preprocessRaw(rgSet) 
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


####### FIX THE FOLLOWING WITH CORRECT ANNOTATION ###########
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
# the densities by sampleID

densityPlot(MSet, sampGroups = pD$sample)

# or density bean plots:

densityBeanPlot(MSet, sampGroups = pD$sample)





########### # Control probes plot

# The 450k array contains several internal control probes that can be used to assess the quality control of different sample preparation 
# steps (bisulfite conversion, hybridization, etc.). The values of these control probes are stored in the initial RGChannelSet and can be 
# plotted by using the function controlStripPlot and by specifying the control probe type:

controlStripPlot(rgSet, controls="BISULFITE CONVERSION II")


qcReport(rgSet, pdf= "qcReport.pdf")   # it's not plotting properly

