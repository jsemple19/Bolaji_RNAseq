
library("rtracklayer")
library("GenomicRanges")
library("ggplot2")
library(readxl)
library(fgsea)
library(GSEAplot)
library(BSgenome.Celegans.UCSC.ce11)
library(dplyr)


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

grp="coh1"
filterOscillating=T

source(paste0(outPath,"/functions.R"))
source(paste0(outPath,"/DESeq2_functions.R"))


clrs<-getColours()
plotPDFs <- F

genomeGR<-GRanges(seqnames=seqnames(Celegans)[1:6], IRanges(start=1, end=seqlengths(Celegans)[1:6]))

ce11seqinfo<-seqinfo(Celegans)

scriptName="compareGeneLists_HumanToWorm"
fileNamePrefix=paste0(scriptName,"/")
makeDirs(outPath,dirNameList=c(paste0(c("plots/"),scriptName)))

padjVal=0.05
lfcVal=0
# DEseq2 results files to use
RNAseqRes<-"/Users/semple/Documents/MeisterLab/otherPeopleProjects/Bolaji/BolajiRNAseq_20211216/rds/coh1_noOsc/coh1_noOsc_COH1vsTEVonly_DESeq2_fullResults.rds"

metadata<-getMetadataGR(genomeDir,genomeVer,outPath)
metadata<-tagOscillating(metadata)
write.table(metadata$wormbaseID[metadata$class=="protein_coding_gene"],
            paste0(outPath,"/allGeneNames.csv"),row.names=F,col.names=F,quote=F)

salmon<-readRDS(RNAseqRes)



#######################-
## ortholist coverage
#######################-

ortholist<-read_excel(paste0(outPath,"/publicData/ortholist_master.xlsx"))
length(unique(ortholist$`WormBase ID`)) #8228
length(unique(ortholist$`HGNC Symbol`)) #11433

names(ortholist)<-make.names(names(ortholist))
numProgs=2
filtOrtholist<- ortholist[ortholist$No..of.Programs>=numProgs,]
dim(filtOrtholist)
filtOrtholist<-filtOrtholist[filtOrtholist$WormBase.ID %in% metadata$wormbaseID[metadata$oscillating =="no"],]
dim(filtOrtholist)
length(unique(filtOrtholist$WormBase.ID))

# numProgs="best_rev"
# ortholist<-read.delim(paste0(outPath,"/publicData/DIOPTorthologs_worm-human_20230615.tsv"))
# length(unique(ortholist$wormbaseID)) #8228 (ortholist) -> 11821 (diopt)
# length(unique(ortholist$symbol)) #11433 (ortholist) -> 15749 (diopt)
#
# #names(ortholist)<-make.names(names(ortholist))
# length(unique(ortholist$wormbaseID[ortholist$best_score=="Yes"]))  #11821
# length(unique(ortholist$symbol[ortholist$best_score=="Yes"])) #9894
# length(unique(ortholist$symbol[ortholist$best_score_rev=="Yes"])) #15738
#
# #filtOrtholist<-ortholist[ortholist$confidence %in% c("high"),]
# filtOrtholist<-ortholist[ortholist$best_score_rev=="Yes",]
# dim(filtOrtholist)
# #dim(filtOrtholist)
# idx<-filtOrtholist$wormbaseID %in% metadata$wormbaseID[metadata$oscillating=="no"]
# filtOrtholist<-filtOrtholist[idx,]
# dim(filtOrtholist)
#
# length(unique(filtOrtholist$wormbaseID))


salmon<-as.data.frame(salmon)[as.data.frame(salmon)$wormbaseID %in% filtOrtholist$WormBase.ID,]
#salmon<-as.data.frame(salmon)[as.data.frame(salmon)$wormbaseID %in% filtOrtholist$wormbaseID,]

##############-
## Weiss .... Merkenschlager Cornelia de lange paper
##############-
#https://www.nature.com/articles/s41467-021-23141-9
#Neuronal genes deregulated in Cornelia de Lange Syndrome respond to removal and re-expression of cohesin
#Felix D. Weiss, Lesly Calderon, Yi-Fang Wang, Radina Georgieva, Ya Guo, Nevena Cvetesic, Maninder Kaur, Gopuraja Dharmalingam, Ian D. Krantz, Boris Lenhard, Amanda G. Fisher & Matthias Merkenschlager
#Nature Communications volume 12, Article number: 2919 (2021)
options(tibble.width=Inf)
cdl<-read_excel("/Users/semple/Documents/MeisterLab/otherPeopleProjects/Bolaji/BolajiRNAseq_20211216/publicData/Supplementary Data 3. Gene expression from CdLS NeuN+ RNA-seq.xlsx",skip=3 )

head(cdl)

cdl<-dplyr::inner_join(cdl,filtOrtholist,by=join_by("Gene.Symbol"=="HGNC.Symbol"))
#cdl<-inner_join(cdl,filtOrtholist,by=join_by("Gene.Symbol"=="symbol"))

cdl$padj[is.na(cdl$padj)]<-1

write.table(cdl,"~/Downloads/cdlGenes_withOrthologs.txt",sep="\t",
            quote=F,row.names=F)

cdlUp<-cdl[cdl$padj<0.05 & cdl$log2FoldChange>0,]
cdlDown<-cdl[cdl$padj<0.05 & cdl$log2FoldChange<0,]

sigGenes<-list()
sigGenes[["cdlUp"]]<-cdlUp$WormBase.ID
sigGenes[["cdlDown"]]<-cdlDown$WormBase.ID
#sigGenes[["cdlUp"]]<-cdlUp$wormbaseID
#sigGenes[["cdlDown"]]<-cdlDown$wormbaseID


### scc-1 genes
# scc1<-readRDS("/Users/semple/Documents/MeisterLab/papers/Moushumi1/p0.05_lfc0.5_filtCycChrAX/filtCycChrAX_X.wt.wt.0mM_scc16cs_vs_wt_DESeq2_fullResults_p0.05.rds")
# scc1$padj[is.na(scc1$padj)]<-1
#
# #scc1<-inner_join(data.frame(scc1),filtOrtholist,by=join_by("wormbaseID"=="wormbaseID"))
#
# #scc1Up<-data.frame(Wormbase.ID=getSignificantGenes(scc1, padj=padjVal, lfc=lfcVal,
#                                                    namePadjCol="padj",
#                                                    nameLfcCol="log2FoldChange",
#                                                    direction="gt",
#                                                    chr="all", nameChrCol="chr")$wormbaseID)
#
# scc1Down<-data.frame(Wormbase.ID=getSignificantGenes(scc1, padj=padjVal, lfc=lfcVal,
#                                                      namePadjCol="padj",
#                                                      nameLfcCol="log2FoldChange",
#                                                      direction="lt",
#                                                      chr="all", nameChrCol="chr")$wormbaseID)
#
#
#
# sigGenes[["SCC-1cs.up"]]<-unique(scc1Up$Wormbase.ID)
# sigGenes[["SCC-1cs.down"]]<-unique(scc1Down$Wormbase.ID)




lapply(sigGenes,length)

gseaTbl<-list()
leadEdgeTbl<-NULL

# making ranks (ranks go from high to low)
ranks <- salmon$log2FoldChange
names(ranks) <- salmon$wormbaseID
ranks<-sort(ranks,decreasing=T)

pdf(file=paste0(outPath, "/plots/",fileNamePrefix,"gseaBar", grp,
                "_Vs_CorneliaDeLange_",numProgs,"progs.pdf"),
    title=,
    width=11,height=8,paper="a4r")
par(mfrow=c(2,1))

colidx<-ifelse(names(ranks) %in% sigGenes[["cdlUp"]],"red","black")
barplot(ranks,col=colidx,border=NA,main=paste0("cdlUp orthologs ",numProgs," program"))

colidx<-ifelse(names(ranks) %in% sigGenes[["cdlDown"]],"red","black")
barplot(ranks,col=colidx,border=NA,main=paste0("cdlDown orthologs ",numProgs," programs"))
#barplot(ranks[names(ranks) %in% sigHS[["COH-1cs.up"]]],col="red",border=NA,add=T)
dev.off()
par(mfrow=c(1,1))
#countidx<-ifelse(names(ranks) %in% sigHS[["COH-1cs.up"]],1,0)
#posidx<-1:length(ranks)

#lines(1:length(ranks),cumsum(countidx)/length(sigHS[["COH-1cs.up"]]))

fgseaRes <- fgsea(sigGenes, ranks, minSize=5, maxSize = 5000)

#fgseaRes[order(padj), ]
fgseaRes

gseaTbl<-plotGseaTable(sigGenes, ranks, fgseaRes,gseaParam=0.1,render=F,
                              colwidths = c(5, 3, 0.8, 0, 1.2))
pathName<-unlist(fgseaRes[,"pathway"])[2]
pathway<-sigGenes[[pathName]]

stats<-ranks

#' Custom plotEnrichment
#'
#' Based on code for plotEnrichment function from fgsea package
#' @param pathway list of pathways where each item as named by the pathway it
#' contains and has a vector of gene names.
#' @param stats Statistical feature used to sort the genes e.g. LFC, ordered by
#' magnitude
#' @param gseaParam gsea parameter (?)
#' @param ticksize size of tick
#' @return plot
#' @export
plotEnrichment1<-function (pathway, stats, fgseaRes, pathName, gseaParam = 1, ticksSize = 0.2){
  toPlotLFC<-data.frame(x=1:length(stats),y=stats)
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  g <- ggplot(toPlot, aes(x = x, y = y)) + geom_point(color = "green",
                                                      size = 0.1) +
    geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") +
    geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") +
    geom_hline(yintercept = 0, colour = "black") +
    geom_line(color = "green") + theme_bw() +
    geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), size = ticksSize,
                 alpha=0.6,color="darkgreen") +
    theme(panel.border = element_blank(), panel.grid.minor = element_blank()) +
    labs(x = "Rank", y = "Enrichment score")
  g<-g+labs(title=paste0("COH-1 cleavage enrichment of ",pathName, " genes")) +
    geom_vline(xintercept=sum(sort(ranks)>0), colour="grey40")+
    annotate("text", x=length(ranks)*0.75, y=fgseaRes[pathway==pathName,ES]*0.9,
             label= paste(paste(paste(c("padj","NES","size"),
                                      fgseaRes[pathway==pathName,c(round(padj,3),round(NES,1),size)],
                                      sep=":"),collapse=", ")))
  g1<-ggplot(toPlotLFC,aes(x=x, y=y)) + geom_bar(stat="identity") +
    xlab("Rank of Log2 ( COH-1cs/TEVcontrol ) High->Low") + ylab("Log2FC")
  p<-ggpubr::ggarrange(g1,g,nrow=2,heights=c(1,2))
  p

  }



gseaList<-list()
#barplot(sort(ranks, decreasing = T))
for (pathName in unlist(fgseaRes[,"pathway"])) {
  gseaList[[pathName]]<-plotEnrichment1(sigGenes[[pathName]], ranks,fgseaRes,pathName) #+
    # labs(title=paste0("CdLS enrichment of ",pathName, " genes")) +
    # geom_vline(xintercept=sum(sort(ranks)>0), colour="grey40")+
    # annotate("text", x=length(ranks)*0.75, y=fgseaRes[pathway==pathName,ES]*0.9,
    #          label= paste(paste(paste(c("padj","NES","size"),
    #          fgseaRes[pathway==pathName,c(round(padj,3),round(NES,1),size)],
    #                                   sep=":"),collapse=", ")))
}
pdf(file=paste0(outPath, "/plots/",fileNamePrefix,"gsea", grp,
                "_Vs_CorneliaDeLange_",numProgs,"progs.pdf"),
    title=,
    width=5,height=10,paper="a4")
p<-gridExtra::marrangeGrob(grobs=gseaList, nrow=2, ncol=1,top=paste0("orthologs found with at least ",numProgs," programs"))
print(p)
dev.off()

####################-
## correlation
####################-

sigUp<-data.frame(getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                                                  namePadjCol="padj",
                                                  nameLfcCol="log2FoldChange",
                                                  direction="gt",
                                                  chr="all", nameChrCol="chr"))


sal<-inner_join(as.data.frame(salmon),ortholist[ortholist$No..of.Programs>=numProgs,],by=join_by("wormbaseID"=="WormBase.ID"))

ij<-inner_join(sal,cdl,by=join_by("HGNC.Symbol"=="Gene.Symbol"))

plot(ij$log2FoldChange.y,ij$log2FoldChange.x,xlab=)
abline(v=0,h=0)

# Bin size control + color palette
ggplot(ij, aes(x=log2FoldChange.y, y=log2FoldChange.x) ) +
  geom_bin2d(bins = 150) +
  scale_fill_continuous(type = "viridis") +
  theme_bw() + ylab(label="CdLS LFC") + xlab(label="COH-1cs LFC")+
  geom_hline(yintercept=0) + geom_vline(xintercept=0)


with(ij[which(abs(ij$log2FoldChange.x) > 0.1 | abs(ij$log2FoldChange.y) > 0.1),],cor(log2FoldChange.x,log2FoldChange.y,method = "spearman"))
with(ij[which(abs(ij$log2FoldChange.x) > 0.1 | abs(ij$log2FoldChange.y) > 0.1),],cor(log2FoldChange.x,log2FoldChange.y,method = "pearson"))


##################-
## GSEAplot package -----
##################-

# need to input tables of counts
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150118
#salmon<-readRDS(RNAseqRes)
#salmon<-as.data.frame(salmon)[as.data.frame(salmon)$wormbaseID %in% filtOrtholist$WormBase.ID,]

cnts<-readRDS(paste0(outPath,"/rds/coh1_counts.rds"))

#cnts<-cnts[rownames(cnts) %in% filtOrtholist$WormBase.ID,]
cnts<-cnts[rownames(cnts) %in% filtOrtholist$wormbaseID,]

coh1cols<-grepl("^828",colnames(cnts))
cnts<-cbind(cnts[,coh1cols],cnts[,!coh1cols])


options(timeout = max(300, getOption("timeout")))

cdl<-read_excel("/Users/semple/Documents/MeisterLab/otherPeopleProjects/Bolaji/BolajiRNAseq_20211216/publicData/Supplementary Data 3. Gene expression from CdLS NeuN+ RNA-seq.xlsx",skip=3 )

head(cdl)
cdl$padj[is.na(cdl$padj)]<-1
#cdl<-inner_join(cdl,filtOrtholist,by=join_by("Gene.Symbol"=="HGNC.Symbol"))

cdl<-inner_join(cdl,filtOrtholist,by=join_by("Gene.Symbol"=="symbol"))


cdlUp<-cdl[cdl$padj<0.05 & cdl$log2FoldChange>0,]
cdlDown<-cdl[cdl$padj<0.05 & cdl$log2FoldChange<0,]

sigGenes<-list()
#sigGenes[["cdlUp"]]<-cdlUp$WormBase.ID
#sigGenes[["cdlDown"]]<-cdlDown$WormBase.ID
sigGenes[["cdlUp"]]<-cdlUp$wormbaseID
sigGenes[["cdlDown"]]<-cdlDown$wormbaseID


lapply(sigGenes,length)

#filter table for genes with orthologs in worm and make expression input
expr.input<-cnts
dim(expr.input)

# make phenotype input
phenotypes<-ifelse(grepl("^828", colnames(cnts)),"COH-1cs","TEVCtrl")
pheno.input=create_phenoinput(phenotypes)
pheno.input
gene.set.input<-create_geneset_db(sigGenes)

pp = GSEAplots(input.ds.name=expr.input,
               input.cls.name=pheno.input,
               gene.set.input=gene.set.input,
               doc.string="GSEA_plots",
               nperm=1000,
               abs.val=F,
               bar_percent=0.1,
               gap_percent=0.1,
               under_percent=0.1,
               upper_percent=0.1,
               color_line="black",
               color_tick="black")

names(pp)

pheno.input
pp$report1
pp$report2


phe=2

data1<-pp$ES[[phe]]
enrich_ind=which(data1$EStag==1)
d=data.frame(x=enrich_ind, y=matrix(min(data1$RES)-0.12,length(enrich_ind),1),
             vx=matrix(-0.8,length(enrich_ind),1), vy=matrix(0.04,length(enrich_ind),1))
p <- ggplot(data1, aes(index,RES))+geom_line(col="black")
p <- p+geom_segment(data=d, mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy),col="black",alpha=0.3)
p <- p+theme_classic()
p <- p+ggtitle(names(pp$gene.set.reference.matrix)[[phe]])
print(p)

# calculate signal2noise
s2n<-list()
s2n[[paste0("mean_",pheno.input$phen[1])]]<-rowMeans(expr.input[pheno.input$class.v==0,])
s2n[[paste0("mean_",pheno.input$phen[2])]]<-rowMeans(expr.input[pheno.input$class.v==1,])
s2n[[paste0("sd_",pheno.input$phen[1])]]<-apply(expr.input[pheno.input$class.v==0,],1,sd)
s2n[[paste0("sd_",pheno.input$phen[2])]]<-apply(expr.input[pheno.input$class.v==1,],1,sd)
s2n<-data.frame(do.call(cbind,s2n))
head(s2n)

s2n$sd_COH.1cs <- ifelse(0.2 * abs(s2n$mean_COH.1cs) < s2n$sd_COH.1cs, s2n$sd_COH.1cs, 0.2 * abs(s2n$mean_COH.1cs))
s2n$sd_COH.1cs <- ifelse(s2n$sd_COH.1cs == 0, 0.2, s2n$sd_COH.1cs)
s2n$sd_TEVCtrl <- ifelse(0.2 * abs(s2n$mean_TEVCtrl) < s2n$sd_TEVCtrl, s2n$sd_TEVCtrl, 0.2 * abs(s2n$mean_TEVCtrl))
s2n$sd_TEVCtrl <- ifelse(s2n$sd_TEVCtrl == 0, 0.2, s2n$sd_TEVCtrl)

signal2noise<-(s2n$mean_COH.1cs-s2n$mean_TEVCtrl)/(s2n$sd_COH.1cs+s2n$sd_TEVCtrl)
names(signal2noise)<-row.names(expr.input)

idxOrd<-match(data1$gene.symbol,names(signal2noise))

head(expr.input[idxOrd,],10)
tail(expr.input[idxOrd,],10)

head(signal2noise[idxOrd],10)
tail(signal2noise[idxOrd],10)
head(data1,10)
tail(data1,10)
