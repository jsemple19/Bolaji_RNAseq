#!/usr/bin/env Rscript
library(DESeq2)
library(Organism.dplyr)
library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)
library(tximport)
library(GenomicFeatures)
library(ggplot2)
library(RColorBrewer)
library(PoiClaClu)
library(pheatmap)
library(tidyr)
library(EnhancedVolcano)
library(affy)
library(gplots)
library(ggpubr)
library(R.utils)
library(rtracklayer)
library(forcats)
library(DEGreport)
library(stringr)


options(bitmapType='cairo')


####
# input arguments
####
args <- commandArgs(trailingOnly=TRUE)
genomeVer="WS285"
if (length(args)<3) {
  ##   stop("Required arguments: <fastqList_csv_file> <output_folder>", call.=FALSE)
  fastqList_file <- "./fastqList.csv"
  outPath <- "."
  scriptDir<-"."
  genomeDir <- paste0("/Users/semple/Documents/MeisterLab/GenomeVer/",genomeVer)
  print(paste("missing arguments will use default: ", fastqList_file, outPath, genomeDir))
} else {
  fastqList_file <- args[1]
  outPath <- args[2]
  genomeDir <- args[3]
}

source(paste0(scriptDir,"/functions.R"))
source(paste0(scriptDir,"/DESeq2_functions.R"))

print(paste("fastqList_file is: ",fastqList_file))
print(paste("outPath is: ",outPath))
print(paste("genomeDir is: ",genomeDir))

####
## other preset variables
####
namePrefix <- "top1top2" #use _noOsc in prefix to automatically filter oscillating genes
fileNamePrefix <- paste0(namePrefix,"/",namePrefix,"_")
plotPDFs <- F

genomeGR<-GRanges(seqnames=seqnames(Celegans)[1:6], IRanges(start=1, end=seqlengths(Celegans)[1:6]))

ce11seqinfo<-seqinfo(Celegans)

makeDirs(outPath,dirNameList=paste0(c("rds","plots","txt","tracks"),"/",namePrefix))

#################-
# sampleTable------
#################-

# Note, sampleTable should have at least the following three columns:
# "fileName" with path to salmon data
# "sampleName" with unique name of sample
# a column with a grouping variable to designate different groups being compared
# (this can have any name and then just pass the name as a string to varOI)
# It is often easiest to create the sampleTable from the fastqList.csv file
# used for mapping.
fileList<-read.table(fastqList_file,stringsAsFactors=F,header=T, sep = ",")
names(fileList)<-c("fileName1","fileName2","sampleID")

sampleNames<-fileList$sampleID
sampleGroup<-factor(gsub("_[1|2|3]$","",fileList$sampleID))
sampleGroup<-relevel(sampleGroup,ref = "N2")

#TIR1<-factor(ifelse(grepl("^T",fileList$sampleID),"TIR1","noTIR1"),levels=c("noTIR1","TIR1"))
auxin<-factor(ifelse(grepl("_p_",fileList$sampleID),"Auxin","Ctrl"),levels=c("Ctrl","Auxin"))

replicate<-factor(stringr::str_extract(fileList$sampleID,"[1|2|3]$"))

fileNames<-paste0(outPath,"/salmon/mRNA/",sampleNames,"/quant.sf")

sampleTable<-data.frame(fileName=fileNames,sampleName=sampleNames,
                        sampleGroup,auxin,replicate)
sampleTable$degron<-"Ctrl"
sampleTable$degron[grep("^TOP1",sampleTable$sampleGroup)]<-"TOP1deg"
sampleTable$degron[grep("^TOP2",sampleTable$sampleGroup)]<-"TOP2deg"
sampleTable$degron<-factor(sampleTable$degron)

sampleTable<-sampleTable[grep("^T",sampleTable$sampleName),]
sampleTable<-droplevels(sampleTable)
sampleTable

varOI<-"sampleGroup"
groupsOI<-sampleTable[,varOI]
metadata<-getMetadataGR(genomeDir,genomeVer,outPath)
metadata<-tagOscillating(metadata)
#highlight specific genes in MA plots:
highlight<-c(metadata$wormbaseID[metadata$publicID %in% c("top-1","top-2")])
names(highlight)<-c(metadata$publicID[metadata$publicID %in% c("top-1","top-2")])
minBatchSize=3

###############################################################-
# Create metadata object --------------------------------------------------
###############################################################-
### create metadata: sumary counts per gene OBJECT needed for Salmon
###############################################################-


tx2gene<-getTx2Gene(genomeDir,genomeVer)


###############################################################-
# Import into DESeq2 ------------------------------------------------------
###############################################################-

# import the count matrices
txi<-tximport(sampleTable$fileName,type="salmon",tx2gene=tx2gene)

#modelTxt<-"~ degron + degron:replicate + degron:auxin"
#reducedModelTxt=" ~ degron + degron:replicate + auxin"
modelTxt<-"~ replicate + auxin + degron + degron:auxin"
reducedModelTxt=" ~ replicate + auxin + degron"

dds <- DESeqDataSetFromTximport(txi=txi,
                                colData=sampleTable,
                                design= formula(modelTxt))

# TODO replace (lane) + SMC with appropriate for this case
rownames(dds)<-gsub("^Gene:","",rownames(dds))

###############################################################-
### DESeq2 differential expression analysis (using negative binomial distribution)
###############################################################-

#dds<-collapseReplicates(dds,groupby=samples$sampleID,renameCols=T)
#dds <- DESeq(dds)
# This function performs a default analysis through the steps:
#   1. estimation of size factors: estimateSizeFactors
#   2. estimation of dispersion: estimateDispersions
#   3. Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
#   returns a DESeqDataSet object

idx<-match(rownames(dds),metadata$wormbaseID)
# add gene and chormosome names as metadata
featureData <- data.frame(gene=rownames(dds),
                          chr=as.vector(seqnames(metadata))[idx]) #,

rowData(dds) <- DataFrame(mcols(dds), featureData)
if(grepl("_noOsc",fileNamePrefix)){
  dds<-removeOscillating(dds,remove="either")
  print("removing oscillating genes")
}
keep <- rowSums(counts(dds) >= 10) >= minBatchSize
dds <- dds[keep,]
dds<-DESeq(dds)

########################
## get contrastsOI ------
########################
colnames(coef(dds))
colMeans(coef(dds),na.rm=T)
contrastsOI<-list()
modMat<-model.matrix(formula(modelTxt),data=sampleTable)
colnames(modMat)

TIRcontrol<-colMeans(modMat[sampleTable$sampleGroup == "TIR1_m", ])
TIRauxin<-colMeans(modMat[sampleTable$sampleGroup == "TIR1_p", ])

Top1control<-colMeans(modMat[sampleTable$sampleGroup == "TOP1_m",])
Top1auxin<-colMeans(modMat[sampleTable$sampleGroup == "TOP1_p",])

Top2control<-colMeans(modMat[sampleTable$sampleGroup == "TOP2_m",])
Top2auxin<-colMeans(modMat[sampleTable$sampleGroup == "TOP2_p",])

contrastsOI[["TIR1AUXvsTIR1NoAUX"]]<-TIRauxin-TIRcontrol
contrastsOI[["Top1xAUX"]]<-Top1auxin-Top1control - (TIRauxin-TIRcontrol) #to get interaction effect
contrastsOI[["Top2xAUX"]]<-Top2auxin-Top2control- (TIRauxin-TIRcontrol) #to get interaction effect.
contrastsOI[["Top1AUXvsTop1"]]<-Top1auxin-Top1control
contrastsOI[["Top2AUXvsTop2"]]<-Top2auxin-Top2control
contrastsOI[["Top1AUXvsTIRAUX"]]<-Top1auxin-TIRauxin
contrastsOI[["Top2AUXvsTIRAUX"]]<-Top2auxin-TIRauxin
contrastNames<-names(contrastsOI)


######################################################-
# Basic sample stats ------------------------------------------------------
######################################################-

## basic sample stats
writeDDSstatsToFile(dds)


#######-
## sample counts summary: boxplots and density plots -----------------------
#######-

compareSampleCounts(dds,varOI,outPath,fileNamePrefix)


#########-
## Model QC: Size factor estimates and dispersion estimates ----------------
#########-
modelQC(dds,outPath,fileNamePrefix)

#########-
## sample to sample and top expressed genes by sample heatmaps ------------------------
#########-
clusteringSamples_heatmap(dds,outPath,fileNamePrefix)


###########-
## covariates -----------------------------------------------------------------
###########-
checkCovariates(dds,outPath,fileNamePrefix)

###########-
## pca ---------------------------------------------------------------------
###########-
plotPCAbyCovariates(dds,outPath,fileNamePrefix)


##########-
## pairwise correlation between genes in different replicates --------------
##########-
print("plotting pairwise correlation between genes")
plotCorrelations(dds,outPath,fileNamePrefix,groupingVariable=varOI)


##############################################################-
# Significant genes -------------------------------------------------------
##############################################################-

print("finding significant genes")

padjVal=0.05
lfcVal=0.5


# Contrast between different groups: groupsOI vs controlGrp
grp=contrastNames[2]
for(grp in contrastNames){
  res<-NULL
  resLFC<-NULL
  if(is.character(contrastsOI[[grp]])){
    res<-results(dds,contrast=list(contrastsOI[[grp]]),alpha=padjVal)
  } else {
    res<-results(dds,contrast=contrastsOI[[grp]],alpha=padjVal)
  }

  resLFC<-getShrunkenLFC(dds, res, metadata, contrastName=grp, shrinkMethod="ashr",
                           outPath=outPath, fileNamePrefix=fileNamePrefix)
  #######-
  ## plot pvalue and filtering threshold QC -----------------------------
  #######-

  plotPvalQC(resLFC,contrastName=grp,outPath=outPath, fileNamePrefix=fileNamePrefix)


  ##########-
  ## heirarchical clustering of most significantly changed genes -------------
  ##########-
  # remove NAs
  resLFC<-na.omit(resLFC)

  plotHClust_significant(dds,resLFC,contrastName=grp,padjVal, outPath,
                                   fileNamePrefix)

  ##########-
  ## plot individual genes -------
  ##########-
  ## Bar plot of 20 most significantly changed genes

  plotCountsMostSig_barplot(dds,resLFC,contrastName=grp,groupingVariable=varOI,
                            outPath=outPath,fileNamePrefix=fileNamePrefix)


  ###########-
  ## MAplot ALL genes -----
  ############-

  makeMAplots(res,resLFC,contrastName=grp, padjVal, highlight ,outPath,
                        fileNamePrefix)

  #############-
  ## Box plot and barplots by chromosome -----
  #############-

  plotLFCbyChromosome(resLFC,contrastName=grp, padjVal ,outPath,fileNamePrefix)
  plotCountsPerChr(resLFC,contrastName=grp, padjVal, lfcVal,outPath,fileNamePrefix)


  #############-
  ## Volcano plots --------
  #############-
  #https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
  makeVolcanoPlot(resLFC,contrastName=grp, padjVal, lfcVal,outPath,fileNamePrefix)


  #############-
  ## LFC tracks bigwig bedgraph and bed --------
  #############-
  makeLFCtracks(resLFC,contrastName=grp, padjVal, lfcVal,outPath,fileNamePrefix)

}



