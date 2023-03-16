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
lfcVal=0.5
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

daughertyEnh<-daughertyEnh %>% as_tibble() %>%
  dplyr::select(L3_chromHMMState,nearestGene,distanceToNearest) %>%
  dplyr::group_by(nearestGene) %>%
  dplyr::summarise(enhancerTypes=paste0(unique(L3_chromHMMState),collapse=","),
            enhancerCount=dplyr::n(),
            enhancerUniqCount=length(unique(L3_chromHMMState)),
            closestEnhancer=min(distanceToNearest),
            furthestEnhancer=max(distanceToNearest))
dim(daughertyEnh)


daughertyEnh$type<-daughertyEnh$enhancerTypes
#daughertyEnh[daughertyEnh$enhancerUniqCount!="1", "type"]<-"L3_mixedEnhancerTypes"
daughertyEnh[daughertyEnh$enhancerUniqCount!="1" &
                grepl("L3_activeEnhancer", daughertyEnh$enhancerTypes), "type"]<-"L3_mixedActive&Repressed"
daughertyEnh[daughertyEnh$enhancerUniqCount!="1" &
                !grepl("L3_activeEnhancer", daughertyEnh$enhancerTypes),"type"]<-"L3_mixedRepressed"
table(daughertyEnh$type)
colnames(daughertyEnh)[colnames(daughertyEnh)=="nearestGene"]<-"wormbaseID"
options(pillar.width = Inf)
head(daughertyEnh)


grp="COH1cs"
#subset<-c("AMA", "COH1cs", "TOP1", "TOP2", "TOP1_1h", "TOP1_2h", "TOP2_1h", "TOP2_2h")
subset=c("COH1cs")

listtbl<-list()
for(grp in subset){
  salmon<-readRDS(file=fileList[fileList$sampleName==grp,"filePath"])
  salmon<-salmon[!is.na(salmon$chr),]
  salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)
  salmongr<-sort(salmongr)
  tbl<-left_join(as_tibble(salmongr),daughertyEnh,by="wormbaseID")
  tbl$type[is.na(tbl$type)]<-"NoEnhancer"
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
#all genes
p1<-ggplot(fulltbl,aes(x=type,y=log2FoldChange,fill=type)) +
  #geom_violin(width=1) +
  geom_boxplot(width=0.5,outlier.colour="grey90",alpha=0.7) +
  ggtitle("LFC of genes closest to Daugherty et al. (2017) L3 enhancers") +
  facet_wrap(.~sampleName) +
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

