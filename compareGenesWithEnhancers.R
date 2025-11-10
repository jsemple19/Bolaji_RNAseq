library(rtracklayer)
library(ggplot2)
library(BSgenome.Celegans.UCSC.ce11)
library(tidyr)
library(ggpubr)
library(dplyr)
library(seqplots)
library(RColorBrewer)
library(ggtext)

outPath="."
scriptPath="."
source(paste0(scriptPath,"/functions.R"))

theme_set(
  theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y=ggtext::element_markdown())
)

scriptName="compareGenesAtEnhancers"
fileNamePrefix=paste0(scriptName,"/")
makeDirs(outPath,dirNameList=c(paste0(c("plots/"),scriptName)))

padjVal=0.05
#lfcVal=0.5
# DEseq2 results files to use
RNAseqDir="/Users/semple/Documents/MeisterLab/otherPeopleProjects/Bolaji/BolajiRNAseq_20211216"

fileList<-data.frame(
  sampleName=c("AMA","FLAVO","COH1cs","TOP1","TOP2","TOP1_0h","TOP1_0.5h","TOP1_1h",
               "TOP1_2h","TOP2_0h","TOP2_0.5h","TOP2_1h", "TOP2_2h"),
  filePath=paste0(RNAseqDir,
                  c("/rds/chem_noOsc/chem_noOsc_AMAvsN2_DESeq2_fullResults.rds",
                    "/rds/chem_noOsc/chem_noOsc_FLAVOvsN2_DESeq2_fullResults.rds",
                    "/rds/coh1_noOsc/coh1_noOsc_COH1vsTEVonly_DESeq2_fullResults.rds",
                    "/rds/top1top2_noOsc/top1top2_noOsc_Top1AUXvsTIRAUX_DESeq2_fullResults.rds",
                    "/rds/top1top2_noOsc/top1top2_noOsc_Top2AUXvsTIRAUX_DESeq2_fullResults.rds",
                    "/rds/top1top2moraoT_noOsc/top1top2moraoT_noOsc_Top1vsAux0_DESeq2_fullResults.rds",
                    "/rds/top1top2moraoT_noOsc/top1top2moraoT_noOsc_Top1vsAux0.5_DESeq2_fullResults.rds",
                    "/rds/top1top2moraoT_noOsc/top1top2moraoT_noOsc_Top1vsAux1_DESeq2_fullResults.rds",
                    "/rds/top1top2moraoT_noOsc/top1top2moraoT_noOsc_Top1vsAux2_DESeq2_fullResults.rds",
                    "/rds/top1top2moraoT_noOsc/top1top2moraoT_noOsc_Top2vsAux0_DESeq2_fullResults.rds",
                    "/rds/top1top2moraoT_noOsc/top1top2moraoT_noOsc_Top2vsAux0.5_DESeq2_fullResults.rds",
                    "/rds/top1top2moraoT_noOsc/top1top2moraoT_noOsc_Top2vsAux1_DESeq2_fullResults.rds",
                    "/rds/top1top2moraoT_noOsc/top1top2moraoT_noOsc_Top2vsAux2_DESeq2_fullResults.rds")))


#' Function for adding count data to plots
count_data <- function (y,ymax=5){
  df <- data.frame(y = ymax, label = length(y))
  return(df)
}



#########################-
## Daugherty ATAC enhancers -----
#########################-

#Genome Res. 2017 Dec;27(12):2096-2107.
# doi: 10.1101/gr.226233.117. Epub 2017 Nov 15.
# Chromatin accessibility dynamics reveal novel functional enhancers in C. elegans
# Aaron C Daugherty 1, Robin W Yeo 1, Jason D Buenrostro 1, William J Greenleaf 1 2, Anshul Kundaje 1 3, Anne Brunet 1 4
# https://pubmed.ncbi.nlm.nih.gov/29141961/
# Suplementary table S3 with added nearest gene and distance to gene

#' boxplot of LFC by enhancer type
daughertyEnh<-readRDS(paste0(outPath,"/publicData/daugherty2017_L3enhancers_ce11.rds"))

table(seqnames(daughertyEnh))
head(daughertyEnh)
length(daughertyEnh)

daughertyEnhByGene<-daughertyEnh %>% as_tibble() %>%
  dplyr::select(L3_chromHMMState,nearestGene,distanceToNearest) %>%
  dplyr::group_by(nearestGene) %>%
  dplyr::summarise(enhancerTypes=paste0(unique(L3_chromHMMState),collapse=","),
            enhancerCount=dplyr::n(),
            enhancerUniqCount=length(unique(L3_chromHMMState)),
            closestEnhancer=min(distanceToNearest),
            furthestEnhancer=max(distanceToNearest))
dim(daughertyEnhByGene)


daughertyEnhByGene$type<-daughertyEnhByGene$enhancerTypes
#daughertyEnhByGene[daughertyEnhByGene$enhancerUniqCount!="1", "type"]<-"L3_mixedEnhancerTypes"
daughertyEnhByGene[daughertyEnhByGene$enhancerUniqCount!="1" &
                grepl("L3_activeEnhancer", daughertyEnhByGene$enhancerTypes), "type"]<-"L3_mixedActive&Repressed"
daughertyEnhByGene[daughertyEnhByGene$enhancerUniqCount!="1" &
                !grepl("L3_activeEnhancer", daughertyEnhByGene$enhancerTypes),"type"]<-"L3_mixedRepressed"
table(daughertyEnhByGene$type)
colnames(daughertyEnhByGene)[colnames(daughertyEnhByGene)=="nearestGene"]<-"wormbaseID"
options(pillar.width = Inf)
head(daughertyEnhByGene)


grp="COH1cs"
#subset<-c("AMA", "COH1cs", "TOP1", "TOP2", "TOP1_1h", "TOP1_2h", "TOP2_1h", "TOP2_2h")
subset=c("COH1cs")

listtbl<-list()
for(grp in subset){
  salmon<-readRDS(file=fileList[fileList$sampleName==grp,"filePath"])
  salmon<-salmon[!is.na(salmon$chr),]
  salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)
  salmongr<-sort(salmongr)
  tbl<-left_join(as_tibble(salmongr),daughertyEnhByGene,by="wormbaseID")
  tbl$type[is.na(tbl$type)]<-"NoEnhancer"
  tmpgr<-resize(width=1,fix="start",makeGRangesFromDataFrame(tbl[tbl$type=="NoEnhancer",]))
  dtn<-distanceToNearest(tmpgr,daughertyEnh)
  tbl[tbl$type=="NoEnhancer","closestEnhancer"]<-mcols(dtn)$distance
  tbl$XvA<-ifelse(tbl$seqnames=="chrX","chrX","autosomes")
  #tbl <- dplyr::filter(tbl,tbl$padj < padjVal)
  #tbl$upVdown<-ifelse(tbl$log2FoldChange > lfcVal,"up","down")
  tbl$sampleName<-grp
  listtbl[[grp]]<-tbl
}

fulltbl<-do.call(rbind,listtbl)
fulltbl$sampleName<-factor(fulltbl$sampleName,levels=subset)
fulltbl$type<-factor(fulltbl$type,levels=c("NoEnhancer",
                                           paste0("L3_",c("activeEnhancer",
                                                        "mixedActive&Repressed",
                                                        "repressedEnhancer",
                                                        "H3K27me3Repressed",
                                                        "mixedRepressed"))))
table(fulltbl$type)
fulltbl %>% dplyr::group_by(type) %>% summarise(minDist=min(closestEnhacer),
                                                avrDist=mean(closestEnhacer),
                                                maxDist=max(closestEnhancer))
#all genes
p1<-ggplot(fulltbl[fulltbl$padj<padjVal,],aes(x=type,y=log2FoldChange,fill=type)) +
  #geom_violin(width=1) +
  geom_boxplot(width=0.5,outlier.colour="grey90",alpha=0.7) +
  ggtitle("LFC of genes closest to Daugherty et al. (2017) L3 enhancers") +
  facet_wrap(.~upVdown) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1,size=10)) +
  stat_summary(fun.data = count_data,fun.args=c(ymax=0.5), geom = "text", position = position_dodge(1),size=5, colour="grey30") +
  geom_hline(yintercept=0,col="grey10",linetype="dashed") +
  coord_cartesian(ylim=c(-0.5,0.5)) +
  scale_fill_manual(values=c("darkgrey","red","purple","blue","darkgreen","#146C94")) +
  ylab(label="Log<sub>2</sub>FC")
p1

# sig genes
p2<-ggplot(fulltbl[fulltbl$padj<padjVal,],aes(x=type,y=log2FoldChange,fill=type)) +
  geom_violin(alpha=0.7) +
  ggtitle("LFC of genes padj<0.05 closest to Daugherty et al. (2017) L3 enhancers")+
  facet_grid(.~sampleName) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1,size=10))+
  stat_summary(fun.data = count_data, fun.args=c(ymax=4),geom = "text", position = position_dodge(1),size=5)+
  geom_hline(yintercept=0,col="grey10",linetype="dashed") +
  scale_fill_manual(values=c("darkgrey","red","purple","blue","darkgreen","#146C94"))+
  coord_cartesian(ylim=c(-4,4)) +
  ylab(label="Log<sub>2</sub>FC")
p2

p<-ggarrange(p1,p2,nrow=2)
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"LFCgenesClosestToEnh.pdf"),plot=p, device="pdf",width=19,height=29,units="cm")


fulltbl$upVdown<-NA
fulltbl$upVdown[fulltbl$log2FoldChange>0]<-"up"
fulltbl$upVdown[fulltbl$log2FoldChange<0]<-"down"
p3<-ggplot(fulltbl[fulltbl$padj<padjVal,],aes(x=upVdown,fill=type)) +
  geom_bar(stat="count",alpha=0.7) +
  #ggtitle("LFC of genes padj<0.05 closest to Daugherty et al. (2017) L3 enhancers")+
  theme(axis.title.x=element_blank(),legend.position="none")+
  #theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1,size=10))+
  #stat_summary(fun.data = count_data, fun.args=c(ymax=4),geom = "text", position =       position_dodge(1),size=5)+
  geom_text(stat='count', aes(label=..count..),position = position_stack(vjust = 0.5))+
  #geom_hline(yintercept=0,col="grey10",linetype="dashed") +
  scale_fill_manual(values=c("darkgrey","red","purple","blue","darkgreen","#146C94"))#+
  #xlab(label=element_blank())
  #coord_cartesian(ylim=c(-4,4)) +
  #ylab(label="Log<sub>2</sub>FC")
p3

fulltbl$shorttype<-fulltbl$type
levels(fulltbl$shorttype)<-gsub("^L3_","",levels(fulltbl$shorttype))
p3a<-ggplot(fulltbl[fulltbl$padj<padjVal,],aes(x=upVdown,fill=shorttype)) +
  geom_bar(stat="count",alpha=1) +
  #ggtitle("LFC of genes padj<0.05 closest to Daugherty et al. (2017) L3 enhancers")+
  theme(axis.title.x=element_blank(),legend.position="none",
        strip.text = element_text(size = 9))+
  facet_grid(.~shorttype)+
  #theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1,size=10))+
  #stat_summary(fun.data = count_data, fun.args=c(ymax=4),geom = "text", position =       position_dodge(1),size=5)+
  geom_text(stat='count', aes(label=..count..),#position = position_stack(vjust = 1),
            nudge_y=15)+
  #geom_hline(yintercept=0,col="grey10",linetype="dashed") +
  scale_fill_manual(values=c("darkgrey","red","purple","blue","darkgreen","#146C94"))#+
#xlab(label=element_blank())
#coord_cartesian(ylim=c(-4,4)) +
#ylab(label="Log<sub>2</sub>FC")
p3a
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"CountgenesClosestToEnh.pdf"),plot=p3a, device="pdf",width=21,height=10,units="cm")

lfcVal=0
resLFC<-fulltbl

#' Do volcano plot
#'
#' Volcano plot with -log10(pvalue) vs log2(foldChange)
#' @param resLFC DESeqResults object after shrinking
#' @param subsetName
#' @param padjVal Adjusted p value to be used for choosing significant genes
#' @return Plot object
#' @export
makeVolcanxoSubPlot<-function(resLFC,subsetName, padjVal){
  # get point colours
  keyvals<-rep('lightgrey', nrow(resLFC))
  names(keyvals)<-rep('NS',nrow(resLFC))

  keyvals[which(resLFC$type == "NoEnhancer")]<-"darkgrey"
  names(keyvals)[which(resLFC$type == "NoEnhancer")]<-"NoEnhancer"
  keyvals[which(resLFC$type == "L3_activeEnhancer")]<-"red"
  names(keyvals)[which(resLFC$type == "L3_activeEnhancer")]<-"L3_activeEnhancer"
  keyvals[which(resLFC$type == "L3_mixedActive&Repressed")]<-"purple"
  names(keyvals)[which(resLFC$type == "L3_mixedActive&Repressed")]<-"L3_mixedActive&Repressed"
  keyvals[which(resLFC$type == "L3_repressedEnhancer")]<-"blue"
  names(keyvals)[which(resLFC$type == "L3_repressedEnhancer")]<-"L3_repressedEnhancer"
  keyvals[which(resLFC$type == "L3_H3K27me3Repressed")]<-"darkgreen"
  names(keyvals)[which(resLFC$type == "L3_H3K27me3Repressed")]<-"L3_H3K27me3Repressed"
  keyvals[which(resLFC$type == "L3_mixedRepressed")]<-"#146C94"
  names(keyvals)[which(resLFC$type == "L3_mixedRepressed")]<-"L3_mixedRepressed"

  #alphaVals<-keyvals
  #alphaVals[names(alphaVals)=="NoEnhancer"]<-0.3
  #alphaVals[names(alphaVals)!="NoEnhancer"]<-0.5

  sigUp<-sum(resLFC$padj[resLFC$type==subsetName]<padjVal &
               resLFC$log2FoldChange[resLFC$type==subsetName]>0)
  sigDown<-sum(resLFC$padj[resLFC$type==subsetName]<padjVal &
                 resLFC$log2FoldChange[resLFC$type==subsetName]<0)

  p4<-EnhancedVolcano(resLFC[resLFC$type==subsetName,],
                      lab=rownames(resLFC[resLFC$type==subsetName,]),
                      labSize=0.5,
                      labCol="#11111100",
                      x="log2FoldChange",
                      y="padj",
                      selectLab=rownames(resLFC[resLFC$type==subsetName,])[1],
                      xlim=c(-2,2),
                      ylim=c(0,40),
                      title= paste0(subsetName),
                      titleLabSize = 12,
                      subtitle=paste0(sigUp," genes up, ",sigDown," genes down."),
                      subtitleLabSize=9,
                      caption=NULL,
                      pCutoff=padjVal,
                      FCcutoff=0,
                      ylab=bquote(~-Log[10]~adjusted~italic(P)),
                      legendPosition = NULL,
                      legendLabSize = 0,
                      legendIconSize = 0,
                      axisLabSize=10,
                      colCustom=keyvals[resLFC$type==subsetName],
                      colAlpha=0.5,#as.numeric(alphaVals)[resLFC$type==subsetName],
                      pointSize = 1)
  return(p4)
}

plotList<-list()
for(subsetName in levels(resLFC$type)){
  plotList[[subsetName]]<-makeVolcanoSubPlot(resLFC,subsetName,padjVal)
}

p<-ggarrange(plotlist=plotList,ncol=3,nrow=2)
p
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"Volcano_genesClosestToEnh.pdf"),plot=p, device="pdf",width=29,height=19,units="cm")


cts<-readRDS(paste0(outPath,"/rds/coh1_counts.rds"))
geomMean<-function(x){
  exp(mean(log(x[x>0])))
}

gm<-data.frame("TEVonly"=apply(cts[,grep("366",colnames(cts))],1,geomMean),
               "COH1cs"=apply(cts[,grep("828",colnames(cts))],1,geomMean))
gm$wormbaseID=rownames(gm)

ggplot(gm,aes(x=log(TEVonly),y=log(COH1cs))) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  geom_abline(slope=1,intercept=0,col="red")
plot(log(gm$TEVonly),log(gm$COH1),pch=16,col="lightgrey")
abline(a=0,b=1)

pdf(file=paste0(outPath, "/plots/",fileNamePrefix,"Scatter_genesClosestToEnh.pdf"),
    paper="a4",width=8,height=11)

par(mfrow=c(3,2))
ss<-gm[rownames(gm) %in% fulltbl$wormbaseID[fulltbl$type=="NoEnhancer" &
                                              fulltbl$padj<padjVal],]
plot(log(ss$TEVonly),log(ss$COH1),pch=16,col="darkgrey",main="NoEnhancer",
     xlim=c(0,12),ylim=c(0,12))
abline(a=0,b=1)

ss<-gm[rownames(gm) %in% fulltbl$wormbaseID[fulltbl$type=="L3_activeEnhancer" &
                                              fulltbl$padj<padjVal],]
plot(log(ss$TEVonly),log(ss$COH1),pch=16,col="red",main="L3_activeEnhancer",
     xlim=c(0,12),ylim=c(0,12))
abline(a=0,b=1)

ss<-gm[rownames(gm) %in% fulltbl$wormbaseID[fulltbl$type=="L3_mixedActive&Repressed"&
                                              fulltbl$padj<padjVal],]
plot(log(ss$TEVonly),log(ss$COH1),pch=16,col="purple",main="L3_mixedActive&Repressed",
     xlim=c(0,12),ylim=c(0,12))
abline(a=0,b=1)

ss<-gm[rownames(gm) %in% fulltbl$wormbaseID[fulltbl$type=="L3_repressedEnhancer" &
                                              fulltbl$padj<padjVal],]
plot(log(ss$TEVonly),log(ss$COH1),pch=16,col="blue",main="L3_repressedEnhancer",
     xlim=c(0,12),ylim=c(0,12))
abline(a=0,b=1)

ss<-gm[rownames(gm) %in% fulltbl$wormbaseID[fulltbl$type=="L3_H3K27me3Repressed" &
                                              fulltbl$padj<padjVal],]
plot(log(ss$TEVonly),log(ss$COH1),pch=16,col="darkgreen",main="L3_H3K27me3Repressed",
     xlim=c(0,12),ylim=c(0,12))
abline(a=0,b=1)

ss<-gm[rownames(gm) %in% fulltbl$wormbaseID[fulltbl$type=="L3_mixedRepressed" &
         fulltbl$padj<padjVal],]
plot(log(ss$TEVonly),log(ss$COH1),pch=16,col="#146C94", main="L3_mixedRepressed",
     xlim=c(0,12),ylim=c(0,12))
abline(a=0,b=1)
dev.off()

par(mfrow=c(1,1))


pdf(file=paste0(outPath, "/plots/",fileNamePrefix,"DensitySig_genesClosestToEnh.pdf"),
    paper="a4",width=8,height=11)

par(mfrow=c(3,2))
ss<-gm[rownames(gm) %in% fulltbl$wormbaseID[fulltbl$type=="NoEnhancer" &
                                              fulltbl$padj<padjVal],]
p1<-ggplot(ss,aes(x=log(TEVonly),y=log(COH1cs))) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  geom_abline(slope=1,intercept=0,col="red") + coord_cartesian(xlim=c(2,10),ylim=c(2,10))+
  ggtitle("NoEnhancer")


ss<-gm[rownames(gm) %in% fulltbl$wormbaseID[fulltbl$type=="L3_activeEnhancer" &
                                              fulltbl$padj<padjVal],]
p2<-ggplot(ss,aes(x=log(TEVonly),y=log(COH1cs))) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  geom_abline(slope=1,intercept=0,col="red") + coord_cartesian(xlim=c(2,10),ylim=c(2,10))+
  ggtitle("L3_activeEnhancer")


ss<-gm[rownames(gm) %in% fulltbl$wormbaseID[fulltbl$type=="L3_mixedActive&Repressed"&
                                              fulltbl$padj<padjVal],]
p3<-ggplot(ss,aes(x=log(TEVonly),y=log(COH1cs))) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  geom_abline(slope=1,intercept=0,col="red") + coord_cartesian(xlim=c(2,10),ylim=c(2,10))+
  ggtitle("L3_mixedActive&Repressed")


ss<-gm[rownames(gm) %in% fulltbl$wormbaseID[fulltbl$type=="L3_repressedEnhancer" &
                                              fulltbl$padj<padjVal],]
p4<-ggplot(ss,aes(x=log(TEVonly),y=log(COH1cs))) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  geom_abline(slope=1,intercept=0,col="red") + coord_cartesian(xlim=c(2,10),ylim=c(2,10))+
  ggtitle("L3_repressedEnhancer")


ss<-gm[rownames(gm) %in% fulltbl$wormbaseID[fulltbl$type=="L3_H3K27me3Repressed" &
                                              fulltbl$padj<padjVal],]
p5<-ggplot(ss,aes(x=log(TEVonly),y=log(COH1cs))) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  geom_abline(slope=1,intercept=0,col="red") + coord_cartesian(xlim=c(2,10),ylim=c(2,10))+
  ggtitle("L3_H3K27me3Repressed")

ss<-gm[rownames(gm) %in% fulltbl$wormbaseID[fulltbl$type=="L3_mixedRepressed" &
                                              fulltbl$padj<padjVal],]
p6<-ggplot(ss,aes(x=log(TEVonly),y=log(COH1cs))) +
  stat_density_2d(aes(fill = ..level..), geom = "raster")+
  geom_abline(slope=1,intercept=0,col="red") + coord_cartesian(xlim=c(2,10),ylim=c(2,10))+
  ggtitle("L3_mixedRepressed")
p6
p<-ggarrange(p1,p2,p3,p4,p5,p6,nrow=3,ncol=2)
print(p)
dev.off()

par(mfrow=c(1,1))

