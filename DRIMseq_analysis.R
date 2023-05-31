#https://ycl6.gitbook.io/guide-to-rna-seq-analysis/differential-expression-analysis/differential-transcript-usage/dtu-using-drimseq
#http://bioconductor.org/packages/release/bioc/vignettes/DRIMSeq/inst/doc/DRIMSeq.pdf
#https://bioconductor.org/packages/release/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html#ref-Nowicka2016DRIMSeq

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
fileNamePrefix=paste0("DRIMseq/")
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
## DRIMseq analysis ------
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


designFormula<- formula("~ lane + replicate + group")

design = model.matrix(designFormula, data = sampleTable)

set.seed(12345)

d <- dmPrecision(d, design = design)    # Estimate the precision (Higher dispersion is associated with lower precision)
d <- dmFit(d, design = design)          # Fit regression coefficients
d <- dmTest(d, coef = "group828")     # Perform null hypothesis testing on the coefficient of interest

saveRDS(d,file=paste0(outPath,"/",fileNamePrefix,"dmDSdata_coh-1.RDS"))


# Single p-value per gene
res.g = DRIMSeq::results(d)
table(is.na(res.g$pvalue))

# Single p-value per transcript
res.t = DRIMSeq::results(d, level = "feature")
table(is.na(res.t$pvalue))

# set NAs to 1
no.na <- function(x) ifelse(is.na(x), 1, x)
res.g$pvalue <- no.na(res.g$pvalue)
res.t$pvalue <- no.na(res.t$pvalue)
res.g$adj_pvalue <- no.na(res.g$adj_pvalue)
res.t$adj_pvalue <- no.na(res.t$adj_pvalue)


table(res.g$adj_pvalue < 0.05)
table(res.t$adj_pvalue < 0.05)


#res.g[res.g$adj_pvalue <0.05,]
#res.t[res.t$adj_pvalue <0.05,]
res.t<-left_join(res.t,as.data.frame(txMeta),by=join_by(feature_id==txptSeqID))

write.table(res.g[res.g$adj_pvalue <0.05,], file=paste0(outPath,"/",fileNamePrefix,"resg_DRIMSeq_results.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
write.table(res.t[res.t$adj_pvalue <0.05,], file=paste0(outPath,"/",fileNamePrefix,"rest_DRIMSeq_results.txt"), sep = "\t", quote = F, row.names = F, col.names = T)

res.t[res.t$adj_pvalue <0.05,]

## Post-hoc filtering on the standard deviation in proportions -----
filt = smallProportionSD(d)

res.t.filt = DRIMSeq::results(d, level = "feature")
res.t.filt$pvalue[filt] = 1
res.t.filt$adj_pvalue[filt] = 1
res.t.filt$pvalue <- no.na(res.t.filt$pvalue)
res.t.filt$adj_pvalue <- no.na(res.t.filt$adj_pvalue)


table(res.t$adj_pvalue < 0.05)
table(res.t.filt$adj_pvalue < 0.05)
res.t.filt[res.t.filt$adj_pvalue <0.05,]


## stageR procedure -----
pScreen = res.g$pvalue

names(pScreen) = res.g$gene_id

pConfirmation = matrix(res.t.filt$pvalue, ncol = 1)

dimnames(pConfirmation) = list(res.t.filt$feature_id)

stageRObj = stageRTx(pScreen = pScreen,
                     pConfirmation = pConfirmation,
                     pScreenAdjusted = FALSE,
                     tx2gene = tx2gene)

stageRObj = stageWiseAdjustment(stageRObj, method = "dtu", alpha = 0.05)

drim.padj = getAdjustedPValues(stageRObj, order = FALSE, onlySignificantGenes = TRUE)
drim.padj = merge(tx2gene, drim.padj, by.x = c("GENEID","TXNAME"), by.y = c("geneID","txID"))

length(unique(drim.padj[drim.padj$gene < 0.05,]$GENEID))
table(drim.padj$transcript < 0.05)

## Exporting results ------
# Create a data.frame containing counts in long-format data with reshape2::melt
drim.prop = reshape2::melt(counts[counts$feature_id %in% proportions(d)$feature_id,], id = c("gene_id", "feature_id"))
drim.prop = drim.prop[order(drim.prop$gene_id, drim.prop$variable, drim.prop$feature_id),]


# Calculate proportions from counts
drim.prop = drim.prop %>%
  group_by(gene_id, variable) %>%
  mutate(total = sum(value)) %>%
  group_by(variable, .add=TRUE) %>%
  mutate(prop = value/total)


# Convert the data.frame to wide-format data with reshape2::dcast
drim.prop = reshape2::dcast(drim.prop[,c(1,2,3,6)], gene_id + feature_id ~ variable)


# Average proportions calculated from fitted proportions
drim.mean = as.data.frame(sapply( levels(sampleTable$group),
                                  function(lvl) rowMeans(proportions(d)[, 3:ncol(proportions(d))][, sampleTable$group == lvl, drop = FALSE]) ))

# log2 fold change in proportions
drim.log2fcp = log2(drim.mean[2]/drim.mean[1])
colnames(drim.log2fcp) = "log2fcp"
rownames(drim.log2fcp) = proportions(d)$feature_id

# Merge to create result data
drimData = cbind(drim.prop[,1:2], drim.mean, drim.prop[, 3:ncol(drim.prop)])
drimData = merge(data.frame(txMeta), drimData, by.x = c("wormbaseID","txptSeqID"), by.y = c("gene_id","feature_id"))
drimData = drimData[order(drimData$wormbaseID, drimData$txptSeqID),]

head(drim.log2fcp)

drimData[1:10,]



# Merge to create result data
drimDTU = merge(drim.padj[,c("GENEID","TXNAME","gene","transcript")], drim.log2fcp, by.x = "TXNAME", by.y = "row.names")
drimDTU = merge(as.data.frame(txMeta), drimDTU, by.x = c("wormbaseID","txptSeqID"), by.y = c("GENEID", "TXNAME"))
drimDTU = drimDTU[order(drimDTU$wormbaseID, drimDTU$txptSeqID),]

write.table(drimData, file=paste0(outPath,"/",fileNamePrefix,"DTU_DRIMSeq-stageR_means_and_proportions.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
write.table(drimDTU, file=paste0(outPath,"/",fileNamePrefix,"DTU_DRIMSeq-stageR_results.txt"), sep = "\t", quote = F, row.names = F, col.names = T)


# Plot the estimated proportions for one of the significant genes, where we can see evidence of switching
#gene_id = unique(drim.padj[order(drim.padj$transcript, drim.padj$gene),]$GENEID)[1]

#png("plotProportions.DRIMSeq-stageR.1.png", width=6, height=6, units = "in", res = 300)
#plotProportions(d, gene_id = gene_id,
#                group_variable = "group", plot_type = "boxplot1")
#dev.off()

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
for(gene_id in unique(drimDTU$wormbaseID)){
  publicID<-unique(txMeta$publicID[txMeta$wormbaseID==gene_id])
  p<-plotExpression(drimData, gene_id, sampleTable, isProportion = TRUE)
  gr<-GenomicRanges::reduce(txMeta[txMeta$wormbaseID==gene_id])
  gr<-GenomicRanges::resize(gr,width=width(gr)+2000,fix="center")
  p.txdb<-autoplot(txdb_ce11,aes(fill=strand),which=gr,label.size=4) + xlim(gr)
  p.fount <- tryPlot(fountains,gr,title="Detected fountains",fill="darkblue") + theme_bed()
  p.activeD<-tryPlot(activedaugherty,gr,title="Active enhancers (Daugherty et al.)",fill="darkgreen") + theme_bed()
  p.activeJ<-tryPlot(activejaenes,gr,fill="darkgreen",title="Active enhancers (Jaenes et al.)") + theme_bed()


  # plotList[[gene_id]]<-p.txdb@ggplot/p.fount@ggplot/p.activeD@ggplot/p.activeJ@ggplot/p +
  #   plot_layout(heights=c(10,1,1,1,10)) +
  #   plot_annotation(title=paste0(gene_id,": ",publicID))

  p1<-ggarrange(p.txdb@ggplot,p.fount@ggplot,p.activeD@ggplot,p.activeJ@ggplot,p,ncol=1,
                heights=c(10,1,1,1,10))
  p1<-annotate_figure(p1,top=paste0(gene_id,": ",publicID))
  plotList[[gene_id]]<- p1
}

pp<-marrangeGrob(grobs=plotList,ncol=1,nrow=1)

ggplot2::ggsave(paste0(outPath,"/",fileNamePrefix,"DRIMseq_prop.pdf"),pp,device="pdf",height=29,width=19,units="cm")

length(unique(drimDTU$wormbaseID))
length(unique(drimDTU$txptSeqID[drimDTU$transcript<0.05]))
head(drimDTU)
