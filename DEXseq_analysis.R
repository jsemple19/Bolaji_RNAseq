#https://ycl6.gitbook.io/guide-to-rna-seq-analysis/differential-expression-analysis/differential-transcript-usage/dtu-using-dexseq


library(data.table)
library(GenomicFeatures)
library(tximport)
library(dplyr)
library(reshape2)
library(ggplot2)
library(DRIMSeq)
library(stageR)
library(ggbeeswarm)
library(annotables)
library(tidyverse)
library(AnnotationDbi)
library(ggplot2)
library(ggpubr)
library(BSgenome.Celegans.UCSC.ce11)
library(ggbio)
library(patchwork)
library(gridExtra)
library(DEXSeq)

theme_set(
  theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y=ggtext::element_markdown(),
          axis.title.x=ggtext::element_markdown())
)

outPath="."
genomeVer="WS285"
genomeDir <- paste0("/Users/semple/Documents/MeisterLab/GenomeVer/",genomeVer)
fastqList_file <- "./fastqList_coh1.csv"
fastqList<-read.csv(fastqList_file)

grp="coh1"
filterOscillating=T

source(paste0(outPath,"/functions.R"))
source(paste0(outPath,"/DESeq2_functions.R"))
#source(paste0(outPath,"/Sleuth_functions.R"))


clrs<-getColours()
fileNamePrefix=paste0("DEXseq/all_")
if(filterOscillating){
  fileNamePrefix<-gsub("\\/","_noOsc/noOsc_",fileNamePrefix)
}
makeDirs(outPath,dirname(fileNamePrefix))

print(paste("fastqList_file is: ",fastqList_file))
print(paste("outPath is: ",outPath))
print(paste("genomeDir is: ",genomeDir))

txMeta<-getTxMetadataGR(genomeDir,genomeVer)
#txMeta<-txMeta[txMeta$class=="protein_coding_gene"]
txMeta<-tagOscillating(txMeta)


################-
## DEXseq analysis ------
################-
quantFiles<-file.path("salmon/mRNA",fastqList$sampleName,"quant.sf")
names(quantFiles)<-fastqList$sampleName

# load a txdb of wormbase data and create a tx2gene object
txdb<-loadDb(paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer,
                    ".annotations.sqlite"))
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
tx2gene$TXNAME<-gsub("Transcript:","",tx2gene$TXNAME)
tx2gene$GENEID<-gsub("Gene:","",tx2gene$GENEID)

txi = tximport(quantFiles, type = "salmon", tx2gene = tx2gene,
               txIn = TRUE, txOut = TRUE, countsFromAbundance = "dtuScaledTPM")

## cbuild data frame from counts
cts = txi$counts
cts = cts[rowSums(cts) > 0,]
range(colSums(cts)/1e6)


txdf.sub = tx2gene[match(rownames(cts),tx2gene$TXNAME),]
counts = data.frame(gene_id = txdf.sub$GENEID, feature_id = txdf.sub$TXNAME, cts)
dim(counts)
counts[1:5,]


if(filterOscillating){
  idx<-row.names(counts) %in% txMeta$txptSeqID[txMeta$oscillating=="no"]
  counts<-counts[idx,]
}



# prepare sampleTable and design formula
sampleTable<-data.frame(sample_id=paste0("X",fastqList$sampleName),
                        group=factor(fastqList$strain,levels=c("366","828")),
                        path=quantFiles,
                        lane=fastqList$lane,
                        replicate=fastqList$repeatNum)

#Pre-filter with DRIMSeq dmFilter
# create a dmDSdata object d
d = dmDSdata(counts = counts, samples = sampleTable)
table(table(counts(d)$gene_id))
#10247 genes have only one transcript
#    1     2     3     4     5     6     7     8     9    10    11    12    13
#10247  2194   687   294   146    68    41    37    20    10    13     8     4
#14    15    16    17    22
#3     1     1     1     1

# filtering
n = nrow(sampleTable)
n.small = min(table(sampleTable$group))

d = dmFilter(d,
             min_samps_feature_expr = n.small, min_feature_expr = 5,
             min_samps_feature_prop = n.small, min_feature_prop = 0.05,
             min_samps_gene_expr = n, min_gene_expr = 5)

table(table(counts(d)$gene_id))
#   2    3    4    5    6    7    8
#1381  369  111   47   13    4    2
#~2000 genes left after filtering

## DEXSeq procedure -----
countData = round(as.matrix(counts(d)[,-c(1:2)]))

dxd = DEXSeqDataSet(countData = countData, sampleData = sampleTable,
                    design = ~ lane + replicate + group + exon + group:exon,
                    featureID = counts(d)$feature_id, groupID = counts(d)$gene_id)

dxd
#6229 genes


dxd = estimateSizeFactors(dxd)
dxd = estimateDispersions(dxd)
dxd = testForDEU(dxd, reducedModel = ~lane + replicate + group + exon)

dxr = DEXSeqResults(dxd, independentFiltering = FALSE)

head(dxr)

qval = perGeneQValue(dxr)
dxr.g = data.frame(gene = names(qval), qval)
dxr.t = as.data.frame(dxr[, c("featureID","groupID","pvalue","padj")])

dim(dxr.g)
head(dxr.g)

dim(dxr.t)
head(dxr.t)

dim(dxr.g[dxr.g$qval < 0.05,])

dim(dxr.t[dxr.t$padj < 0.05,])
length(unique(dxr.t$groupID[dxr.t$padj < 0.05]))

write.table(dxr.g[dxr.g$padj <0.05,], file=paste0(outPath,"/",fileNamePrefix,"dxrg_DEXSeq_results.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
write.table(dxr.t[dxr.t$padj <0.05,], file=paste0(outPath,"/",fileNamePrefix,"dxrt_DEXSeq_results.txt"), sep = "\t", quote = F, row.names = F, col.names = T)

dxr.t<-read.table(file=paste0(outPath,"/",fileNamePrefix,"dxrt_DEXSeq_results.txt"), sep = "\t", header = T)

dxr.t[dxr.t$padj <0.05,]
dim(dxr.t)
length(unique(dxr.t$groupID))

## stageR procedure -----
pScreen = qval
pScreen

pConfirmation = matrix(dxr.t$pvalue, ncol=1)
dimnames(pConfirmation) = list(dxr.t$featureID,"transcript")


tx2gene = data.frame(dxr.t[,c("featureID", "groupID")],
                     dxr.t[,c("featureID", "groupID")])
for (i in 1:2) tx2gene[,i] = tx2gene[,i]


# qval used in pScreen, hence pScreenAdjusted=TRUE
stageRObj = stageRTx(pScreen = pScreen,
                     pConfirmation = pConfirmation,
                     pScreenAdjusted = TRUE,
                     tx2gene = tx2gene[1:2])

stageRObj = stageWiseAdjustment(stageRObj, method = "dtu", alpha = 0.05)

dex.padj = getAdjustedPValues(stageRObj, order = FALSE, onlySignificantGenes = TRUE)
dex.padj = merge(tx2gene, dex.padj, by.x = c("groupID","featureID"), by.y = c("geneID","txID"))

length(unique(dex.padj[dex.padj$gene < 0.05,]$groupID))
table(dex.padj$transcript < 0.05)


## Exporting results-----
dex.norm = cbind(as.data.frame(stringr::str_split_fixed(rownames(counts(dxd)), ":", 2)), as.data.frame(counts(dxd, normalized = TRUE)))
colnames(dex.norm) = c("groupID", "featureID", as.character(colData(dxd)$sample_id))
row.names(dex.norm) = NULL

# Per-group normalised mean
dex.mean = as.data.frame(sapply( levels(sampleTable$group),
                                 function(lvl) rowMeans(dex.norm[, 3:ncol(dex.norm)][, sampleTable$group == lvl, drop = FALSE]) ))

# log2 fold change in expression
dex.log2fc = log2(dex.mean[2]/dex.mean[1])
colnames(dex.log2fc) = "log2fc"
rownames(dex.log2fc) = dex.norm$featureID

# Merge to create result data
dexData = cbind(dex.norm[,1:2], dex.mean, dex.norm[, 3:ncol(dex.norm)])
dexData = merge(as.data.frame(txMeta), dexData, by.x = c("wormbaseID","txptSeqID"), by.y = c("groupID","featureID"))
dexData = dexData[order(dexData$wormbaseID, dexData$txptSeqID),]


# Merge to create result data
dexDTU = merge(dex.padj[,c("featureID.1","groupID.1","gene","transcript")], dex.log2fc, by.x = "featureID.1", by.y = "row.names")
dexDTU = merge(txMeta, dexDTU, by.x = c("wormbaseID","txptSeqID"), by.y = c("groupID.1","featureID.1"))
dexDTU = dexDTU[order(dexDTU$wormbaseID, dexDTU$txptSeqID),]

write.table(dexData, file=paste0(outPath,"/",fileNamePrefix,"DTU_DEXSeq-stageR_means_and_counts.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
write.table(dexDTU, file=paste0(outPath,"/",fileNamePrefix,"DTU_DEXSeq-stageR_results.txt"), sep = "\t", quote = F, row.names = F, col.names = T)

dexDTU<-read.table(file=paste0(outPath,"/",fileNamePrefix,"DTU_DEXSeq-stageR_results.txt"), sep = "\t", header=T)

dim(dexDTU)
length(unique(dexDTU$wormbaseID))

unique(dexDTU$wormbaseID) %in% dxr.t$groupID

ggvenn(list(DEXseq=unique(dxr.t$groupID), stageR=unique(dexDTU$wormbaseID)))

########## comparisons of results -----
ggvenn(list(DRIMseqStageR=unique(drimDTU$wormbaseID), DEXseqStageR=unique(dexDTU$wormbaseID)))

ggvenn(list(DEXseq=unique(dxr.t$groupID), DRIMseq=unique(res.t$wormbaseID)))
# Plot the estimated proportions for one of the significant genes, where we can see evidence of switching

ggvenn(list(sleuth=unique(upGenes$wormbaseID),DRIMseqStageR=unique(drimDTU$wormbaseID), DEXseqStageR=unique(dexDTU$wormbaseID)))

ggvenn(list(sleuth=unique(upGenes$wormbaseID),DEXseq=unique(dxr.t$groupID), DRIMseq=unique(res.t$wormbaseID)))


unique(dxr.t$groupID)[unique(dxr.t$groupID) %in% unique(res.t$wormbaseID)]
unique(drimDTU$wormbaseID)[unique(drimDTU$wormbaseID) %in% unique(dexDTU$wormbaseID)]
###########-



# get fountains
fountains<-readRDS("/Users/semple/Documents/MeisterLab/otherPeopleProjects/fountains/detected_fountains_equalQ.RDS")
fountains$fountainName<-paste0("fount",1:length(fountains))

# get daugherty enhancers
daugherty<-readRDS("/Users/semple/Documents/MeisterLab/otherPeopleProjects/fountains/publicData/daugherty2017_L3enhancers_ce11.rds")
activedaugherty<-daugherty[daugherty$L3_chromHMMState=="L3_activeEnhancer"]

# get jaenes enhancers
jaenes<-import("/Users/semple/Documents/MeisterLab/otherPeopleProjects/fountains/publicData/Jaenes2018_enhancers_ce11_stages_L3chromHMM.bed")
activejaenes<-jaenes[jaenes$name=="Active enhancer"]

theme_bed = function(){
  theme(title=element_text(size=8),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank())
}

#' Try bed track plot but if no overlaps in this region - do empty plot keeping title
tryPlot<-function(gr1,gr2,title,fill){
  p1<-autoplot(subsetByOverlaps(gr1,gr2,ignore.strand=T), geom="rect",fill=fill,legend=F) +
    xlim(gr2) + ggtitle(title)
  p2<-autoplot(Celegans, which = gr2) + xlim(gr2)+ ggtitle(title)
  if(is(try(print(p1)),"try-error")) p2 else p1
}

# get transcripts in ce11 format
if(!file.exists(paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer,
                       "_ce11.annotations.sqlite"))){
  makeTxDbsqlite_ce11(genomeDir,genomeVer)
}
txdb_ce11<-loadDb(paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer,
                         "_ce11.annotations.sqlite"))

plotList<-list()
for(gene_id in unique(dexDTU$wormbaseID)){
  publicID<-unique(txMeta$publicID[txMeta$wormbaseID==gene_id])
  p<-plotExpression(dexData, gene_id, sampleTable, isProportion = FALSE)
  gr<-GenomicRanges::reduce(txMeta[txMeta$wormbaseID==gene_id])
  gr<-GenomicRanges::resize(gr,width=width(gr)+2000,fix="center")
  p.txdb<-autoplot(txdb_ce11,aes(fill=strand),which=gr,label.size=4) + xlim(gr)
  p.fount <- tryPlot(fountains,gr,title="Detected fountains",fill="darkblue") + theme_bed()
  p.activeD<-tryPlot(activedaugherty,gr,title="Active enhancers (Daugherty et al.)",fill="darkgreen") + theme_bed()
  p.activeJ<-tryPlot(activejaenes,gr,fill="darkgreen",title="Active enhancers (Jaenes et al.)") + theme_bed()

  p1<-ggarrange(p.txdb@ggplot,p.fount@ggplot,p.activeD@ggplot,p.activeJ@ggplot,p,ncol=1,
                heights=c(10,1,1,1,10))
  p1<-annotate_figure(p1,top=paste0(gene_id,": ",publicID))
  plotList[[gene_id]]<- p1
}

pp<-marrangeGrob(grobs=plotList,ncol=1,nrow=1)

ggplot2::ggsave(paste0(outPath,"/",fileNamePrefix,"_DTU_DEXseq-stageR_counts.pdf"),pp,device="pdf",height=29,width=19,units="cm")

length(unique(dexDTU$wormbaseID))
length(unique(dexDTU$txptSeqID[dexDTU$transcript<0.05]))
head(dexDTU)


### plot the significant transcripts that were discarded by stageR ----
idx<-dxr.t$groupID[dxr.t$padj<0.05] %in% dexDTU$wormbaseID
noStageR_wbid<-unique(dxr.t[dxr.t$padj<0.05,"groupID"][!idx])

plotList<-list()
for(gene_id in noStageR_wbid){
  publicID<-unique(txMeta$publicID[txMeta$wormbaseID==gene_id])
  p<-plotExpression(dexData, gene_id, sampleTable, isProportion = FALSE)
  gr<-GenomicRanges::reduce(txMeta[txMeta$wormbaseID==gene_id])
  gr<-GenomicRanges::resize(gr,width=width(gr)+2000,fix="center")
  p.txdb<-autoplot(txdb_ce11,aes(fill=strand),which=gr,label.size=4) + xlim(gr)
  p.fount <- tryPlot(fountains,gr,title="Detected fountains",fill="darkblue") + theme_bed()
  p.activeD<-tryPlot(activedaugherty,gr,title="Active enhancers (Daugherty et al.)",fill="darkgreen") + theme_bed()
  p.activeJ<-tryPlot(activejaenes,gr,fill="darkgreen",title="Active enhancers (Jaenes et al.)") + theme_bed()

  p1<-ggarrange(p.txdb@ggplot,p.fount@ggplot,p.activeD@ggplot,p.activeJ@ggplot,p,ncol=1,
                heights=c(10,1,1,1,10))
  p1<-annotate_figure(p1,top=paste0(gene_id,": ",publicID))
  plotList[[gene_id]]<- p1
}

pp<-marrangeGrob(grobs=plotList,ncol=1,nrow=1)

ggplot2::ggsave(paste0(outPath,"/",fileNamePrefix,"_DTU_DEXseq_counts.pdf"),pp,device="pdf",height=29,width=19,units="cm")
