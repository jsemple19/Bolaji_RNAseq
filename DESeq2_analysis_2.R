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

#library("TxDb.Celegans.UCSC.ce11.refGene")
#library("TxDb.Celegans.UCSC.ce11.ensGene")

options(bitmapType='cairo')
#library(REBayes)
# get funciton for converting gene names from WBID to publicID
#CEFTALL_SRC <- Sys.getenv("CEFTALL_SRC", unset = "./src")

#source(paste0(CEFTALL_SRC, "/convertingGeneNamesFunction1.R"))
#source(paste0(CEFTALL_SRC, "/functions.R"))

#CONDA_ACTIVATE <- Sys.getenv("CONDA_ACTIVATE")
#print(paste("env CONDA_ACTIVATE=", CONDA_ACTIVATE))

####
## input arguments
####
args <- commandArgs(trailingOnly=TRUE)
genomeVer="WS285"
if (length(args)<3) {
  ##   stop("Required arguments: <fastqList_csv_file> <output_folder>", call.=FALSE)
  fastqList_file <- "./fastqList.csv"
  # fastqList_file <- "./data/fastqList_t-all_ringo_test.csv"
  # fastqList_file <- "./data/fastqList_t-all_ringo.csv"
  # fastqList_file <- "./data/fastqList_t-all_ringo_kingston.csv"
  outPath <- "."  # "./out"
  #genomeDir <- paste0("./data/",genomeVer)
  genomeDir <- paste0("/Users/semple/Documents/MeisterLab/GenomeVer/",genomeVer)
  print(paste("missing arguments will use default: ", fastqList_file, outPath, genomeDir))
} else {
  fastqList_file <- args[1]
  outPath <- args[2]
  genomeDir <- args[3]
}

#TODO:
#source(paste0(outPath,"/convertingGeneNamesFunction1.R"))
source(paste0(outPath,"/functions.R"))
clrs<-getClrs()

print(paste("fastqList_file is: ",fastqList_file))
print(paste("outPath is: ",outPath))
print(paste("genomeDir is: ",genomeDir))

####
### other preset variables
####


fileNamePrefix <- "allData_"
plotPDFs <- F

#strain_numbers <- c("200", "284", "285")
#strain_names <- c("WT", "mes2_control", "mes2_HLH1")


genomeGR<-GRanges(seqnames=seqnames(Celegans)[1:6], IRanges(start=1, end=seqlengths(Celegans)[1:6]))

wbseqinfo<-seqinfo(Celegans)
seqnames(wbseqinfo)<-c(gsub("chr","",seqnames(Celegans)))
seqnames(wbseqinfo)<-c(gsub("^M$","MtDNA",seqnames(wbseqinfo)))
genome(wbseqinfo)<-genomeVer
ce11seqinfo<-seqinfo(Celegans)

makeDirs(outPath,dirNameList=c("rds","plots","txt","tracks"))


fileList<-read.table(fastqList_file,stringsAsFactors=F,header=T, sep = ",")
names(fileList)<-c("fileName1","fileName2","sampleID")

sampleNames<-fileList$sampleID
sampleGroup<-factor(gsub("_[1|2|3]$","",fileList$sampleID))
sampleGroup<-relevel(sampleGroup,ref = "N2")

TIR1<-factor(ifelse(grepl("^T",fileList$sampleID),"TIR1","noTIR1"),levels=c("noTIR1","TIR1"))
auxin<-factor(ifelse(grepl("_p_",fileList$sampleID),"Auxin","Ctrl"),levels=c("Ctrl","Auxin"))

replicate<-factor(stringr::str_extract(fileList$sampleID,"[1|2|3]$"))

fileNames<-paste0(outPath,"/salmon/mRNA/",sampleNames,"/quant.sf")

sampleTable<-data.frame(fileName=fileNames,sampleName=sampleNames,
                        sampleGroup,TIR1,auxin,replicate,HS="noHS",date="20211216")


coh1List<-read.table("fastqList_coh1.csv",stringsAsFactors=F,header=T, sep = ",")
coh1Table<-data.frame(fileName=paste0(outPath,"/salmon/mRNA/",coh1List$sampleName,"/quant.sf"),
                      sampleName=coh1List$sampleName, sampleGroup=coh1List$strain,
                      TIR1="noTIR",auxin="Ctrl",replicate=coh1List$repeatNum,
                      HS="HS",date=coh1List$date)
coh1Table$sampleGroup<-ifelse(coh1Table$sampleGroup=="828","coh1","TEVonly")


moraoList<-read.table("fastqList_Morao.csv",stringsAsFactors=F,header=T, sep = ",")
strain<-gsub("deg","d",gsub("top-","top",sapply(strsplit(moraoList$sampleName,"_"),"[[",2)))
auxin<-sapply(strsplit(moraoList$sampleName,"_"),"[[",3)
auxin<-gsub("30","0.5",gsub("min","h",gsub("hours?","h",auxin)))
moraoTable<-data.frame(fileName=paste0(outPath,"/salmon/mRNA/",moraoList$sampleName,"/quant.sf"),
                       sampleName=moraoList$sampleName, sampleGroup=paste0(strain,"_",auxin),
                       TIR1="TIR1",auxin=auxin,replicate="1",HS="noHS",date="Morao2022")
moraoTable$replicate[duplicated(moraoTable$sampleGroup)]<-"2"

sampleTable<-rbind(sampleTable,coh1Table,moraoTable)

groupsOI<-sampleTable$sampleGroup
varsOI<-c("TIR1","auxin","sampleGroup","replicate","date","HS")

###############################################################-
# Create metadata object --------------------------------------------------
###############################################################-
### create metadata: sumary counts per gene OBJECT needed for Salmon
###############################################################-

if(!file.exists(paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer, ".annotations.sqlite"))){
  #dir.create(paste0(genomeDir,"/annotations"),recursive=T)
  si<-seqinfo(Celegans)
  genome(si)<-genomeVer
  seqnames(si)<-gsub("M","MtDNA",gsub("chr","",seqnames(si)))
  wstxdb<-makeTxDbFromGFF(file=paste0(genomeDir,
                                      "/c_elegans.PRJNA13758.",
                                      genomeVer,".annotations.gtf"),
                          format="gtf",organism="Caenorhabditis elegans",
                          chrominfo=si)

  saveDb(wstxdb,paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer,
                       ".annotations.sqlite"))
}

# load a txdb of wormbase data and create a tx2gene object
txdb<-loadDb(paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer,
                    ".annotations.sqlite"))
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
tx2gene$TXNAME<-gsub("Transcript:","",tx2gene$TXNAME)

print("tx2gene columns: "); columns(txdb)
print("tx2gene keytypes: "); keytypes(txdb)
TxptByGene<-transcriptsBy(txdb, by = "gene")
print("TxptByGene length: "); length(TxptByGene)


geneGR<-unlist(range(TxptByGene))
mcols(geneGR)$wormbase<-names(geneGR)
genedf<-as.data.frame(geneGR)
genedf$wormbase<-gsub("^Gene:","",genedf$wormbase)

# download gene id data from simplemine: https://wormbase.org/tools/mine/simplemine.cgi
# for entrez ids, load wormbaseID column into https://david.ncifcrf.gov/conversion.jsp
#geneIDs<-read.delim(paste0(genomeDir,"/c_elegans.PRJNA13758.",genomeVer,".geneIDs.txt"), header=FALSE, sep=",")
geneIDs<-read.delim(gzfile(paste0(genomeDir,"/c_elegans.PRJNA13758.",genomeVer,".geneIDs.txt.gz")), header=FALSE, sep=",")
colnames(geneIDs) <- c("species.ID","WormBase.Gene.ID","Public.Name","Sequence.Name","gene.anotation.status","class")
# can be ignored: david
#david<-read.delim("/Users/semple/Documents/MeisterLab/GenomeVer/annotations/david_wbid2entrez_WS278.txt")

### combine wormbase gene IDs with gene location ###
metadata<-inner_join(geneIDs, genedf,by=c("WormBase.Gene.ID"="wormbase")) %>%
  dplyr::select(WormBase.Gene.ID,Public.Name,Sequence.Name,class,seqnames,start, end, strand) %>%
  collect %>% GenomicRanges::GRanges()

names(mcols(metadata))<-c("wormbaseID","publicID","sequenceID","class")

table(metadata$class)
#i<-which(metadata$wormbaseID %in% david$From)
#j<-match(metadata$wormbaseID[i],david$From)
#metadata$entrezID<-NA
#metadata$entrezID[i]<-david$To[j]

#seqinfo(metadata)<-wbseqinfo
seqlevelsStyle(metadata)<-"ucsc"
seqinfo(metadata)<-ce11seqinfo
metadata<-sort(metadata)


prtncd<-metadata[metadata$class=="protein_coding_gene"]
prtncd$name<-prtncd$wormbaseID
export(prtncd,"proteinCodingGenes.bed")

ncd<-metadata[metadata$class!="protein_coding_gene"]
ncd$name<-ncd$class
export(ncd,"nonCodingGenes.bed")

saveRDS(metadata,paste0(outPath,"/ce11GeneGR_",genomeVer,".rds"))

###############################################################-
# Import into DESeq2 ------------------------------------------------------
###############################################################-

# import the count matrices
txi<-tximport(sampleTable$fileName,type="salmon",tx2gene=tx2gene)

modelTxt<-"~ sampleGroup"
#reducedModelTxt=" ~ replicate + time + treatment"
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

#only take expressed genes
#dds <- dds[ rowSums(counts(dds)) > 1, ]
#print("Number of expressed genes:")
#dim(dds)[1]
#16606 genes from 20127

dds<-DESeq(dds)
#dds<-DESeq(dds, betaPrior = TRUE)


### get contrastsOI
colnames(coef(dds))
contrastsOI<-list()
modMat<-model.matrix(formula(modelTxt),data=sampleTable)
colnames(modMat)
control<-colMeans(modMat[sampleTable$sampleGroup== "N2", ])

TIRcontrol<-colMeans(modMat[sampleTable$sampleGroup == "TIR1_m", ])
TIRauxin<-colMeans(modMat[sampleTable$sampleGroup == "TIR1_p", ])

Top1control<-colMeans(modMat[sampleTable$sampleGroup == "TOP1_m",])
Top1auxin<-colMeans(modMat[sampleTable$sampleGroup == "TOP1_p",])

Top2control<-colMeans(modMat[sampleTable$sampleGroup == "TOP2_m",])
Top2auxin<-colMeans(modMat[sampleTable$sampleGroup == "TOP2_p",])

contrastsOI[["TIR1"]]<-TIRauxin-TIRcontrol
contrastsOI[["Top1_PMauxin"]]<-Top1auxin-Top1control
contrastsOI[["Top2_PMauxin"]]<-Top2auxin-Top2control
contrastsOI[["Top1_auxinBG"]]<-Top1auxin-TIRauxin
contrastsOI[["Top2_auxinBG"]]<-Top2auxin-TIRauxin
contrastNames<-names(contrastsOI)

######################################################-
# Basic sample stats ------------------------------------------------------
######################################################-

## basic sample stats
sink(file=paste0(outPath,"/txt/", fileNamePrefix,
                 "all_logfile.txt"),
     append=FALSE, type="output")
statsPerSample<-data.frame(t(apply(counts(dds),2,summary)))
statsPerSample$totalCounts<-colSums(counts(dds))
rownames(statsPerSample)<-colData(dds)$sampleName
colnames(statsPerSample)<-c("min", "Q1", "median", "mean", "Q3", "max","totalCounts")
statsPerSample$zeros <- apply(counts(dds)==0, 2, sum)
statsPerSample$percZeros <- round(100*statsPerSample$zeros/nrow(counts(dds)),1)
print(statsPerSample)
sink()

grpClrs<-c("grey","green","purple","lightblue","blue","pink","red","yellow","orange")
pdf(file=paste0(outPath,"/plots/",fileNamePrefix,
                "sampleQC.pdf"), width=8,height=8,paper="a4")

#######-
## sample counts summary: boxplots and density plots -----------------------
#######-
## box plots
epsilon <- 1 # pseudo-count to avoid problems with log(0)
#df<-as.data.frame(log2(counts(dds) + epsilon))
#colnames(df)<-colData(dds)$sampleName
#ldf<-gather(df,key=sampleName,value=log2Counts)
#ldf$strain<-gsub("_.*$","",ldf$sampleName)
#p<-ggplot(ldf,aes(x=logFC,y=sampleName))+geom_boxplot(aes(fill=strain))
par(mar=c(5,8,4,2))
boxplot(log2(counts(dds) + epsilon), col=grpClrs[as.numeric(groupsOI)],
        pch=".",horizontal=TRUE, cex.axis=1,
        las=1, ylab=NULL, names=colData(dds)$sampleName,
        xlab=paste0("log2(counts ","+",epsilon,")"))
par(mar=c(5,4,4,2))

## density plots
plotDensity(log2(counts(dds) + epsilon), lty=1, lwd=2, xlab="log2 counts per gene",
            col=grpClrs[as.numeric(groupsOI)])
grid()
legend("topright", legend=colData(dds)$sampleName, lwd=2,cex=0.8,
       col=grpClrs[
         as.numeric(groupsOI)])



#########-
## dispersion estimates ----------------------------------------------------
#########-
plotDispEsts(dds)


#########-
## sample to sample heatmap ------------------------------------------------
#########-
vsd <- vst(dds, blind=TRUE)
colnames(vsd)<-colData(dds)$sampleName
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colData(dds)$sampleName
colnames(sampleDistMatrix) <- colData(dds)$sampleName
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#########-
## Heatmap - most highly expressed genes -----------------------------------
#########-  TODO do we need the top 500?
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:500]
df <- data.frame(colData(dds)[,c("sampleGroup","TIR1","auxin","replicate")])
rownames(df)<-colData(dds)$sampleName
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df, main="Top 500 expressed genes - no clustering")

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=TRUE, annotation_col=df, main="Top 500 expressed genes - clustered by columns")



###########-
## pca ---------------------------------------------------------------------
###########-
p1<-plotPCA(vsd, intgroup=c("TIR1"))
print(p1)
p2<-plotPCA(vsd, intgroup=c("replicate"))
print(p2)
p3<-plotPCA(vsd, intgroup=c("auxin"))
print(p3)
p4<-plotPCA(vsd,intgroup=c("sampleGroup"))
print(p4)
p5<-plotPCA(vsd,intgroup=c("date"))
print(p5)
p6<-plotPCA(vsd,intgroup=c("HS"))
print(p6)

###########-
## size factor ---------------------------------------------------------------------
###########-
print("plotting qc from DEGreport package")
counts <- counts(dds, normalized = TRUE)
design <- as.data.frame(colData(dds))
degCheckFactors(counts)
colnames(counts)<-design$sampleName
rownames(design)<-design$sampleName
degCovariates(log2(counts+0.5),design)
degCorCov(design)
dev.off()



##########-
## pairwise correlation between genes in different replicates --------------
##########-
print("plotting pairwise correlation between genes")
#Define a function to draw a scatter plot for a pair of variables (samples) with density colors
plotFun <- function(x,y){
  dns <- densCols(x,y);
  points(x,y, col=dns, pch=".", panel.first=grid());
  #  abline(a=0, b=1, col="brown")
}

corFun <- function(x,y){
  par(usr = c(0, 1, 0, 1))
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 3*cor(x, y),
       col=colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)[9*(cor(x,y))])
}

df<-as.data.frame(log2(counts(dds) + epsilon))
colnames(df)<-colData(dds)$sampleName
for (grp in groupsOI){
  if(plotPDFs==T){
    pdf(file=paste0(outPath,"/plots/",fileNamePrefix, grp, "_correlations.pdf"),
        width=8,height=8,paper="a4")
  } else {
    png(file=paste0(outPath,"/plots/",fileNamePrefix, grp, "_correlations.png"),
        width=8,height=8,units="in",res=150)
  }
  idx<-colData(dds)$sampleGroup %in% c(grp)
  if(sum(idx)>1){
    pairs(df[,idx], panel=plotFun, lower.panel=corFun, labels=colnames(df)[idx], main=grp)
  }
  dev.off()
}
if(plotPDFs==T){
  pdf(file=paste0(outPath,"/plots/",fileNamePrefix, "all_correlations.pdf"),
      width=12,height=12,paper="a4")
} else {
  png(file=paste0(outPath,"/plots/",fileNamePrefix, "all_correlations.png"),
      width=24,height=24,units="in",res=150)
}
tmp<-data.frame(do.call(rbind,strsplit(colnames(df),"_")))
df<-df[,order(tmp$X3,tmp$X2)]
pairs(df, panel=plotFun, lower.panel=corFun, labels=colnames(df), main="All samples")
dev.off()



##############################################################-
# Significant genes -------------------------------------------------------
##############################################################-

print("finding significant genes")
# TODO
pThresh=0.05
LFCthresh=0.5
#res=list()
#resLFC=list()

# Contrast between different groups: groupsOI vs controlGrp
grp=contrastNames[1]
for(grp in contrastNames){
  res<-results(dds,contrast=contrastsOI[[grp]],alpha=pThresh)
  sink(file=paste0(outPath,"/txt/",fileNamePrefix, grp,
                   "_logfile.txt"), append=TRUE,
       type="output")
  cat(paste0("Number of genes that change expression in ",grp," at different padj cutoffs:\n"))
  print(summary(res))
  print(summary(res,alpha=0.05))
  print(summary(res,alpha=0.01))
  sink()

  ### add metadata
  res$wormbaseID<-rownames(res)
  idx<-match(rownames(res),metadata$wormbaseID)
  res$chr<-factor(seqnames(metadata),levels=paste0("chr",c("I","II","III","IV","V","X")))[idx]
  res$start<-as.vector(start(metadata))[idx]
  res$end<-as.vector(end(metadata))[idx]
  res$strand<-as.vector(strand(metadata))[idx]

  # shrink LFC estimates
  #resultsNames(dds) # to get names of coefficients
  resLFC<-lfcShrink(dds, type="ashr", res=res)
  class(resLFC)
  ### add metadata for further processing
  resLFC$wormbaseID<-rownames(resLFC)
  idx<-match(rownames(resLFC),metadata$wormbaseID)
  resLFC$chr<-factor(seqnames(metadata),levels=paste0("chr",c("I","II","III","IV","V","X")))[idx]
  resLFC$start<-as.vector(start(metadata))[idx]
  resLFC$end<-as.vector(end(metadata))[idx]
  resLFC$strand<-as.vector(strand(metadata))[idx]
  resLFC$publicID<-as.vector(metadata$publicID)[idx]
  resLFC$sequenceID<-as.vector(metadata$sequenceID)[idx]
  #resLFC$entrezID<-as.vector(metadata$entrezID)[idx]
  saveRDS(resLFC,file=paste0(outPath,"/rds/", fileNamePrefix, grp,
                             "_DESeq2_fullResults.rds"))

  #export csv with ordered results
  write.csv(resLFC[order(resLFC$padj),],
            file=paste0(outPath,"/txt/", fileNamePrefix,grp,
                        "_DESeq2_resultsTable.csv"),
            quote=F,row.names=F)

  #######-
  ## plot pvalue QC ----------------------------------------
  #######-

  pdf(file=paste0(outPath,"/plots/",fileNamePrefix,grp,
                  "_pval_qc.pdf"), width=8,height=,paper="a4")
  #counts <- counts(dds, normalized = TRUE)
  #design <- as.data.frame(colData(dds))
  # v=varsOI[1]
  # for (v in varsOI){
  #   print(v)
  #   p<-degQC(counts, design[[v]], pvalue = resLFC[["pvalue"]])
  #   print(p)
  # }
  hist(resLFC$pvalue,breaks=9)
  hist(resLFC$padj,breaks=9)

  #######-
  ## plot results filtering threshold ----------------------------------------
  #######-
  plot(metadata(resLFC)$filterNumRej,
       type="b", ylab="number of rejections",
       xlab="quantiles of filter",main="Threshold for independant filtering of results")
  lines(metadata(resLFC)$lo.fit, col="red")
  abline(v=metadata(resLFC)$filterTheta)
  legend("topright",legend=paste0("Read count \nthreshold: ",
                                  round(metadata(resLFC)$filterThreshold,2)))


  dev.off()

  # remove NAs
  res<-na.omit(res)
  resLFC<-na.omit(resLFC)

  pdf(file=paste0(outPath,"/plots/",fileNamePrefix,grp,
                  "_hclust_mostChanged.pdf"), width=8,height=11,paper="a4")

  ##########-
  ## heirarchical clustering of most significantly changed genes -------------
  ##########-
  # select gene names based on FDR (5%) # TODO maybe try different threshhold: log2FoldChange 0.5 ~= 1.4 fold change
  gene.kept <- rownames(resLFC)[resLFC$padj <= pThresh & !is.na(resLFC$padj) & abs(resLFC$log2FoldChange)>0]
  if(length(gene.kept)>10){
    # Retrieve the normalized counts for gene of interest
    countTable.kept <- log2(counts(dds) + epsilon)[gene.kept, ]
    dim(countTable.kept)
    colnames(countTable.kept)<-colData(dds)$sampleName

    # Perform the hierarchical clustering with
    # A distance based on Pearson-correlation coefficient
    # and average linkage clustering as agglomeration criteria
    heatmap.2(as.matrix(countTable.kept),
              scale="row",
              hclust=function(x)hclust(x,method="average"),
              distfun=function(x) as.dist((1-cor(t(x)))/2),
              margin=c(7,0),
              trace="none",
              density="none",
              labRow="",
              #labCol = names(countTable.kept),
              cexCol=1,
              main=paste0(grp," changed genes (p<",pThresh,", |lfc|>0)"))
  }
  dev.off()


  ##########-
  ## plot individual genes ---------------------------------------------------
  ##########-
  ## Bar plot of 20 most signifficant

  pdf(file=paste0(outPath,"/plots/",fileNamePrefix,grp,
                  "_topGenes_normCounts.pdf"), width=8,height=11,paper="a4")
  par(mfrow=c(5,2))
  selectedGenes <- rownames(resLFC)[order(resLFC$padj)][1:20]
  print("selectedGenes are:")
  print(selectedGenes)
  if(length(selectedGenes)>1){
    for (g in selectedGenes) {
      g.title <- paste0(g, ":", geneIDs[which(geneIDs$WormBase.Gene.ID == g),3])
      barplot(counts(dds, normalized=TRUE)[g,],
              col=c(grpClrs)[as.numeric(groupsOI)],
              main=g.title, las=2, cex.names=1,
              names.arg=colData(dds)$sampleName)
      # las:2 perpendicular; substring redundant prefix from sample name:8
      legend("topleft",legend=levels(groupsOI),
             fill=c(grpClrs), cex=0.7)
    }
  }
  dev.off()


  ##########-
  ## make GRanges for LFC ----------------------------------------------------
  ##########-
  #remove nas
  resGR<-GenomicRanges::GRanges(seqnames=resLFC$chr,
                                IRanges::IRanges(start=resLFC$start,
                                                 end=resLFC$end),
                                strand=resLFC$strand)
  seqlengths(resGR)<-seqlengths(Celegans)[1:6]
  mcols(resGR)<-resLFC[,c("wormbaseID","log2FoldChange","padj")]

  names(mcols(resGR))[names(mcols(resGR))=="log2FoldChange"]<-"score"
  resGR<-sort(resGR,ignore.strand=TRUE)

  #https://github.com/hochwagenlab/hwglabr2/wiki/Apply-function-to-GRanges-scores-in-genomic-tiles
  forBG<-resGR
  mcols(forBG)<-mcols(forBG)[,c("wormbaseID","score")]
  colnames(mcols(forBG))<-c("name","score")
  seqinfo(forBG)<-ce11seqinfo
  export(forBG,paste0(outPath,"/tracks/",fileNamePrefix,grp,"_wt_lfc.bedGraph"),
         format="bedGraph")


  forBW<-disjoin(forBG,ignore.strand=T)
  oldf<-as.data.frame(findOverlaps(forBW,forBG,ignore.strand=T))
  oldf$scorePerSubBp<-forBG$score[oldf$subjectHits]/width(forBG)[oldf$subjectHits]
  oldf$scorePerQuery<-width(forBW)[oldf$queryHits]*oldf$scorePerSubBp
  score<-oldf %>% group_by(queryHits) %>% summarise(score=mean(scorePerQuery))
  forBW$score<-score$score
  export(forBW,paste0(outPath,"/tracks/",fileNamePrefix,grp,"_wt_lfc.bw"),
         format="bigwig")

  #######-
  ## bed file for significant genes ------------------------------------------
  #######-

  idx<-which(resGR$padj<0.05)
  forBed<-resGR[idx]
  mcols(forBed)<-mcols(forBed)[,c("wormbaseID","score")]
  colnames(mcols(forBed))<-c("name","score")
  seqinfo(forBed)<-ce11seqinfo
  #NaIdx<-is.na(forBed$score)
  #forBed$score[NaIdx]<-0
  export(forBed,paste0(outPath,"/tracks/",fileNamePrefix,grp,"_wt_lfc_p",gsub("^0.","",pThresh),".bedGraph"),
         format="bedGraph")


  export(forBed,paste0(outPath,"/tracks/",fileNamePrefix,grp,"_wt_lfc_p",gsub("^0.","",pThresh),".bed"),
         format="bed")


  ###########-
  ## MAplot ALL genes
  ############-

  pdf(file=paste0(outPath,"/plots/",fileNamePrefix,grp, "_MAplots_results.pdf"), width=5,height=5,paper="a4")

  hlh1_WBgene <- "WBGene00001948"
  hlh1_BaseMean <- res[hlh1_WBgene,]$baseMean

  plotMA(res, main=paste0(grp," uncorrected LFC, threshold=", pThresh), ylim=c(-8,8), alpha=pThresh, colSig = "red3")
  points(hlh1_BaseMean, res[hlh1_WBgene,]$log2FoldChange, col="darkblue", pch=1)
  text(hlh1_BaseMean, 0.3+res[hlh1_WBgene,]$log2FoldChange, "HLH-1", col="black")
  #abline(h = res[hlh1_WBgene,]$log2FoldChange, col="darkblue")
  #abline(v = hlh1_BaseMean, col="darkblue")

  plotMA(resLFC, main=paste0(grp," apeglm shrunk LFC, threshold=", pThresh), ylim=c(-8,8), alpha=pThresh, colSig = "red3")
  points(hlh1_BaseMean, resLFC[hlh1_WBgene,]$log2FoldChange, col="darkblue", pch=1)
  text(hlh1_BaseMean, 0.3+resLFC[hlh1_WBgene,]$log2FoldChange, "HLH-1", col="black")
  #plotCounts(dds, gene=which.min(res$padj), intgroup="sampleType")
  dev.off()

  #############-
  ## Box plot by chromosome
  #############-
  pdf(file=paste0(outPath,"/plots/",fileNamePrefix, grp,
                  "_boxPlots_expnByChr.pdf"), width=8,height=5,paper="a4")
  chrName<-factor(resLFC$chr)
  geneCounts<-table(chrName)

  boxplot(log2FoldChange~chrName,data=resLFC,varwidth=TRUE,outline=FALSE,notch=TRUE,
          main=paste0("Expression changes in ", grp), ylab="log2 Fold Change",
          col=c(rep("grey",5),"grey"),xlab="chromosome (number of genes)",
          names=paste(names(geneCounts)," \n(",geneCounts,")",sep=""))
  #stripchart(log2FoldChange~chrType,data=res,method="jitter",vertical=TRUE,pch=20,col="#11115511",cex=0.5,add=TRUE)
  abline(h=0,lty=2,col="blue")
  dev.off()


  #############-
  ## Volcano plots -----------------------------------------------------------
  #############-
  #https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
  #all black plot for manual changing of colours.
  #pdf(file=paste0(outPath,"/plots/",fileNamePrefix, grp,
  #                "_volcanoPlot_allGenes.pdf"), width=8,height=6,paper="a4")
  keyvals<-rep('black', nrow(resLFC))
  names(keyvals)<-rep('NS',nrow(resLFC))
  keyvals[which(resLFC$padj<pThresh & abs(resLFC$log2FoldChange)>1)]<-'red'
    names(keyvals)[which(resLFC$padj<pThresh & abs(resLFC$log2FoldChange)>1)]<-paste0('p<',pThresh,' |lfc|>1')
    sigUp<-sum(resLFC$padj<pThresh & resLFC$log2FoldChange>1)
    sigDown<-sum(resLFC$padj<pThresh & resLFC$log2FoldChange< -1)
    p1<-EnhancedVolcano(resLFC,
                        lab=rownames(resLFC),
                        labSize=0.5,
                        labCol="#11111100",
                        x="log2FoldChange",
                        y="padj",
                        selectLab=rownames(resLFC)[12366],
                        xlim=c(-5.5,5.5),
                        ylim=c(0,65),
                        title= paste0(grp),
                        subtitle=NULL,
                        caption = paste0(nrow(resLFC), ' expressed genes. ',sigUp, " up, ",sigDown," down."),
                        captionLabSize = 12,
                        pCutoff=pThresh,
                        FCcutoff=1.0,
                        xlab=bquote(~Log[2]~'fold change'~.(grp)),
                        ylab=bquote(~-Log[10]~adjusted~italic(P)),
                        #.legend=c('NS','P & Log2 FC'),
                        #legendLabels=c('NS', expression(p-value<pThresh~and~log[2]~FC>1)),
                        legendPosition = 'top',
                        legendLabSize = 12,
                        legendIconSize = 3.0,
                        axisLabSize=14,
                        colCustom=keyvals,
                        #col = c("black", "red"),
                        colAlpha=0.5,
                        pointSize = 1.0)
    #dev.off()
    if(plotPDFs==T){
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, grp,
                             "_volcanoPlot_allGenes.pdf"), plot=p1,
             device="pdf", width=12,height=12,units="cm")
    } else {
      ggsave(filename=paste0(outPath,"/plots/",fileNamePrefix, grp,
                             "_volcanoPlot_allGenes.png"), plot=p1,
             device="png", width=12,height=12,units="cm")
    }




    summaryByChr<-function(resLFC,pThresh,LFCthresh) {
      up<-resLFC[resLFC$padj < pThresh & resLFC$log2FoldChange > LFCthresh,]
      down<-resLFC[resLFC$padj < pThresh & resLFC$log2FoldChange < -LFCthresh, ]
      allChr<-as.data.frame(rbind(up=table(up$chr),down=table(down$chr)))
      allChr$autosomes<-rowSums(allChr[,1:5])
      allChr$total<-rowSums(allChr[,1:6])
      rownames(allChr)<-paste0(rownames(allChr),"_p",pThresh,"_lfc",LFCthresh)
      return(allChr)
    }


    sink(file=paste0(outPath,"/txt/", fileNamePrefix, grp,
                     "_logfile.txt"),append=TRUE, type="output")
    cat("Summary by Chr: \n")
    cat("\np=0.05, LFC=0: \n")
    print(summaryByChr(resLFC,pThres=0.05,LFCthresh=0))
    cat("\np=0.05, LFC=0.5: \n")
    print(summaryByChr(resLFC,pThres=0.05,LFCthresh=0.5))
    cat("\np=0.05, LFC=1: \n")
    print(summaryByChr(resLFC,pThres=0.05,LFCthresh=1))

    cat("\np=0.01, LFC=0: \n")
    print(summaryByChr(resLFC,pThres=0.01,LFCthresh=0))
    cat("\np=0.01, LFC=0.5: \n")
    print(summaryByChr(resLFC,pThres=0.01,LFCthresh=0.5))
    cat("\np=0.01, LFC=1: \n")
    print(summaryByChr(resLFC,pThres=0.01,LFCthresh=1))

    sink()
    #summary(resLFC,alpha=0.05)


    ##### create filtered tables of gene names
    results<-readRDS(paste0(outPath,"/rds/",fileNamePrefix, grp,
                            "_DESeq2_fullResults.rds"))
    #results<-na.omit(results)


    nrow(filterResults(results,padj=0.05,lfc=0,"both","all", writeTable=F))
    nrow(filterResults(results,padj=0.05,lfc=0.5,"both","all", writeTable=F))
    nrow(filterResults(results,padj=0.05,lfc=0.75,"both","all", writeTable=F))
    nrow(filterResults(results,padj=0.05,lfc=1,"both","all", writeTable=F))
    nrow(filterResults(results,padj=0.01,lfc=0,"both","all", writeTable=F))
    nrow(filterResults(results,padj=0.01,lfc=0.5,"both","all", writeTable=F))
    nrow(filterResults(results,padj=0.01,lfc=0.75,"both","all", writeTable=F))
    nrow(filterResults(results,padj=0.01,lfc=1,"both","all", writeTable=F))

}




