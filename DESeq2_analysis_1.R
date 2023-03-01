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
clrs<-getColours()

print(paste("fastqList_file is: ",fastqList_file))
print(paste("outPath is: ",outPath))
print(paste("genomeDir is: ",genomeDir))

####
### other preset variables
####


fileNamePrefix <- "allData_"
plotPDFs <- F


genomeGR<-GRanges(seqnames=seqnames(Celegans)[1:6], IRanges(start=1, end=seqlengths(Celegans)[1:6]))

wbseqinfo<-seqinfo(Celegans)
seqnames(wbseqinfo)<-c(gsub("chr","",seqnames(Celegans)))
seqnames(wbseqinfo)<-c(gsub("^M$","MtDNA",seqnames(wbseqinfo)))
genome(wbseqinfo)<-genomeVer
ce11seqinfo<-seqinfo(Celegans)

makeDirs(outPath,dirNameList=c("rds","plots","txt","tracks"))


### create sample Table -----
fileList<-read.table(fastqList_file,stringsAsFactors=F,header=T, sep = ",")
names(fileList)<-c("fileName1","fileName2","sampleID")

sampleNames<-fileList$sampleID
sampleGroup<-factor(gsub("_[1|2|3]$","",fileList$sampleID))
sampleGroup<-relevel(sampleGroup,ref = "N2")
strain<-gsub("_[p|m]$","",sampleGroup)

TIR1<-factor(ifelse(grepl("^T",fileList$sampleID),"TIR","noTIR"),levels=c("noTIR","TIR"))
auxin<-factor(ifelse(grepl("_p_",fileList$sampleID),"Auxin","Ctrl"),levels=c("Ctrl","Auxin"))

replicate<-factor(stringr::str_extract(fileList$sampleID,"[1|2|3]$"))

fileNames<-paste0(outPath,"/salmon/mRNA/",sampleNames,"/quant.sf")

sampleTable<-data.frame(fileName=fileNames,sampleName=sampleNames,strain=strain,
                        sampleGroup,TIR1,auxin,replicate,HS="noHS",date="20211216")


coh1List<-read.table("fastqList_coh1.csv",stringsAsFactors=F,header=T, sep = ",")
coh1Table<-data.frame(fileName=paste0(outPath,"/salmon/mRNA/",coh1List$sampleName,"/quant.sf"),
                      sampleName=coh1List$sampleName,sampleGroup=coh1List$strain,
                      strain=coh1List$strain,
                      TIR1="noTIR",auxin="noAuxin",replicate=coh1List$repeatNum,
                      HS="HS",date=coh1List$date)
coh1Table$sampleGroup<-ifelse(coh1Table$sampleGroup=="828","coh1","TEVonly")


moraoList<-read.table("fastqList_Morao.csv",stringsAsFactors=F,header=T, sep = ",")
strain<-gsub("deg","d",gsub("top-","top",sapply(strsplit(moraoList$sampleName,"_"),"[[",2)))
auxin<-sapply(strsplit(moraoList$sampleName,"_"),"[[",3)
auxin<-gsub("30","0.5",gsub("min","h",gsub("hours?","h",auxin)))
moraoTable<-data.frame(fileName=paste0(outPath,"/salmon/mRNA/",moraoList$sampleName,"/quant.sf"),
                       sampleName=moraoList$sampleName, sampleGroup=paste0(strain,"_",auxin),
                       strain=strain,
                       TIR1="TIR",auxin=auxin,replicate="R1",HS="noHS",date="Morao2022")
moraoTable$replicate[duplicated(moraoTable$sampleGroup)]<-"R2"

sampleTable<-rbind(sampleTable,coh1Table,moraoTable)

groupsOI<-sampleTable$sampleGroup

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

grpClrs<-clrs#paste0(clrs,"3")
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
par(mar=c(5,10,4,2))
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
rowAvr<-rowMeans(assay(dds))
idx<-order(rowAvr,decreasing = T)[1:10000]
vsd <- vst(dds[idx], blind=TRUE)
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
select <- order(rowMeans(assay(vsd,normalized=TRUE)),
                decreasing=TRUE)[1:500]
df <- data.frame(colData(dds)[,3:ncol(colData(dds))])
rownames(df)<-colData(vsd)$sampleName
#pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df, main="Top 500 expressed genes - no clustering")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=TRUE, annotation_col=df, main="Top 500 expressed genes - clustered by columns")



###########-
## size factor ---------------------------------------------------------------------
###########-
print("plotting qc from DEGreport package")
counts <- counts(dds, normalized = TRUE)
design <- as.data.frame(colData(dds))
degCheckFactors(counts) # This function plots the median ratio normalization used by DESeq2 to visually check whether the median is the best size factor to represent depth
colnames(counts)<-design$sampleName
rownames(design)<-design$sampleName
degCovariates(log2(counts+0.5),design) # correlation of PCs with covariates
degCorCov(design) # correlation between all covariates
dev.off()

###########-
## pca ---------------------------------------------------------------------
###########-
pdf(file=paste0(outPath,"/plots/",fileNamePrefix,
                "PCA.pdf"), width=8,height=8,paper="a4")
colourByVar<-colnames(sampleTable)[grep("fileName|sampleName",colnames(sampleTable),invert=T)]
#numPages<-length(colourByVar)/2
close
plotlist<-list()
for(varToUse in colourByVar){
  p1<-plotPCA(vsd, intgroup=c(varToUse),ntop=5000)+ggtitle(varToUse)+
    theme(legend.text = element_text(size=8))
  plotlist[[varToUse]]<-p1
}
p<-ggarrange(plotlist=plotlist,ncol=1,nrow=2)
print(p)
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

#rowAvr<-rowMeans(assay(dds))
#idx<-order(rowAvr,decreasing = T)[1:10000]
#df<-as.data.frame(log2(counts(dds)[idx,] + epsilon))
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



