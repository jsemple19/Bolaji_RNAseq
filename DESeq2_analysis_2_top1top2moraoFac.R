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
  fastqList_file <- "./fastqList_Morao.csv"
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

namePrefix<-"top1top2moraoF_noOsc" #use _noOsc in prefix to automatically filter oscillating genes
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


moraoList<-read.table("fastqList_Morao.csv",stringsAsFactors=F,header=T, sep = ",")
strain<-factor(gsub("deg","d",gsub("top-","top",sapply(strsplit(moraoList$sampleName,"_"),"[[",2))),
               levels=c("CA1200", "AM05-top1d", "SS01A-top2d"))
auxin<-sapply(strsplit(moraoList$sampleName,"_"),"[[",3)
auxin<-gsub("30","0.5",gsub("min","h",gsub("hours?","h",auxin)))
auxin<-factor(auxin,levels=c("0h-auxin","0.5h-auxin","1h-auxin","2h-auxin"))
moraoTable<-data.frame(fileName=paste0(outPath,"/salmon/mRNA/",moraoList$sampleName,"/quant.sf"),
                      sampleName=factor(moraoList$sampleName),
                      strain=strain,
                      sampleGroup=factor(paste0(strain,"_",auxin),
                          levels=paste(rep(levels(strain),each=4),rep(levels(auxin),3),sep="_")),
                      auxin=factor(auxin),replicate="1")
moraoTable$replicate[duplicated(moraoTable$sampleGroup)]<-"2"
moraoTable$replicate<-factor(moraoTable$replicate)
newSampleNames<-paste0(gsub("in$","_r",moraoTable$sampleGroup),moraoTable$replicate)
moraoTable$sampleName<-newSampleNames
sampleTable<-moraoTable
sampleTable$time<-as.numeric(gsub("h-auxin","",sampleTable$auxin))
sampleTable

varOI<-"sampleGroup"
groupsOI<-sampleTable[,varOI]
metadata<-getMetadataGR(genomeDir,genomeVer,outPath)
metadata<-tagOscillating(metadata)
#highlight specific genes in MA plots:
highlight<-c(metadata$wormbaseID[metadata$publicID %in% c("top-1","top-2")])
names(highlight)<-c(metadata$publicID[metadata$publicID %in% c("top-1","top-2")])
minBatchSize=2

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
modelTxt<-"~ replicate + sampleGroup"

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

auxin0<-colMeans(modMat[sampleTable$sampleGroup == "CA1200_0h-auxin", ])
auxin0.5<-colMeans(modMat[sampleTable$sampleGroup == "CA1200_0.5h-auxin", ])
auxin1<-modMat[sampleTable$sampleGroup == "CA1200_1h-auxin", ]
auxin2<-colMeans(modMat[sampleTable$sampleGroup == "CA1200_2h-auxin", ])

top1aux0<-colMeans(modMat[sampleTable$sampleGroup == "AM05-top1d_0h-auxin", ])
top1aux0.5<-colMeans(modMat[sampleTable$sampleGroup == "AM05-top1d_0.5h-auxin", ])
top1aux1<-colMeans(modMat[sampleTable$sampleGroup == "AM05-top1d_1h-auxin", ])
top1aux2<-colMeans(modMat[sampleTable$sampleGroup == "AM05-top1d_2h-auxin", ])

top2aux0<-colMeans(modMat[sampleTable$sampleGroup == "SS01A-top2d_0h-auxin", ])
top2aux0.5<-colMeans(modMat[sampleTable$sampleGroup == "SS01A-top2d_0.5h-auxin", ])
top2aux1<-colMeans(modMat[sampleTable$sampleGroup == "SS01A-top2d_1h-auxin", ])
top2aux2<-colMeans(modMat[sampleTable$sampleGroup == "SS01A-top2d_2h-auxin", ])


contrastsOI[["Top1vsAux0"]]<-top1aux0-auxin0
contrastsOI[["Top1vsAux0.5"]]<-top1aux0.5-auxin0.5
contrastsOI[["Top1vsAux1"]]<-top1aux1-auxin1
contrastsOI[["Top1vsAux2"]]<-top1aux2-auxin2

contrastsOI[["Top1Aux0.5vsTop1Aux0"]]<-top1aux0.5-top1aux0
contrastsOI[["Top1Aux1vsTop1Aux0"]]<-top1aux1-top1aux0
contrastsOI[["Top1Aux2vsTop1Aux0"]]<-top1aux2-top1aux0

contrastsOI[["Top1xAux0.5"]]<-top1aux0.5-top1aux0 - (auxin0.5 - auxin0)
contrastsOI[["Top1xAux1"]]<-top1aux1-top1aux0 - (auxin1 - auxin0)
contrastsOI[["Top1xAux2"]]<-top1aux2-top1aux0 - (auxin2 - auxin0)

contrastsOI[["Top2vsAux0"]]<-top2aux0-auxin0
contrastsOI[["Top2vsAux0.5"]]<-top2aux0.5-auxin0.5
contrastsOI[["Top2vsAux1"]]<-top2aux1-auxin1
contrastsOI[["Top2vsAux2"]]<-top2aux2-auxin2

contrastsOI[["Top2Aux0.5vsTop2Aux0"]]<-top2aux0.5-top2aux0
contrastsOI[["Top2Aux1vsTop2Aux0"]]<-top2aux1-top2aux0
contrastsOI[["Top2Aux2vsTop2Aux0"]]<-top2aux2-top2aux0


contrastsOI[["Top2xAux0.5"]]<-top2aux0.5-top2aux0 - (auxin0.5 - auxin0)
contrastsOI[["Top2xAux1"]]<-top2aux1-top2aux0 - (auxin1 - auxin0)
contrastsOI[["Top2xAux2"]]<-top2aux2-top2aux0 - (auxin2 - auxin0)

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
plotCorrelations(dds,outPath,fileNamePrefix,groupingVariable="strain")


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



