library(GenomicRanges)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(eulerr)

outPath="."
scriptPath="."
source(paste0(scriptPath,"/functions.R"))

theme_set(
  theme_bw()+
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
)

scriptName="compareDatasets"
fileNamePrefix=paste0(scriptName,"/")
makeDirs(outPath,dirNameList=c(paste0(c("plots/"),scriptName)))

padjVal=0.05
lfcVal=0.5
# DEseq2 results files to use
RNAseqDir="/Users/semple/Documents/MeisterLab/otherPeopleProjects/Bolaji/BolajiRNAseq_20211216"

fileList<-data.frame(sampleName=c("AMA","FLAVO","COH1cs","TOP1","TOP2","TOP1_0h","TOP1_0.5h","TOP1_1h",
                        "TOP1_2h","TOP2_0h","TOP2_0.5h","TOP2_1h", "TOP2_2h"),
                      filePath=paste0(RNAseqDir,c("/rds/chem_noOsc/chem_noOsc_AMAvsN2_DESeq2_fullResults.rds",
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




###################-
## plot Morao data as trajectories-----
###################-

moraoTop1<-fileList[grepl("top1top2morao",fileList$filePath) & grepl("TOP1",fileList$sampleName),]
moraoTop2<-fileList[grepl("top1top2morao",fileList$filePath) & grepl("TOP2",fileList$sampleName),]


# get ids that are significant in at lest one time point
rl<-getListOfResults(moraoTop1,padjVal)
lapply(rl,dim)
res<-do.call(rbind,rl)
sigIDs<-unique(res$wormbaseID)
# get full tables and filter by IDs
rl<-getListOfResults(moraoTop1)
lapply(rl,dim)
res<-do.call(rbind,rl)
res<-res[res$wormbaseID %in% sigIDs,]
table(res$sampleName)
res$sampleName<-factor(res$sampleName,levels=paste0("TOP1_",c("0h","0.5h","1h","2h")))

p1<-ggplot(res,aes(x=sampleName,y=log2FoldChange,group=wormbaseID))+
  geom_line(alpha=0.1) +ggtitle("TOP1 trajectories")


# get ids that are significant in at lest one time point
rl<-getListOfResults(moraoTop2,padjVal)
lapply(rl,dim)
res<-do.call(rbind,rl)
sigIDs<-unique(res$wormbaseID)
# get full tables and filter by IDs
rl<-getListOfResults(moraoTop2)
lapply(rl,dim)
res<-do.call(rbind,rl)
res<-res[res$wormbaseID %in% sigIDs,]
table(res$sampleName)
res$sampleName<-factor(res$sampleName,levels=paste0("TOP2_",c("0h","0.5h","1h","2h")))

p2<-ggplot(res,aes(x=sampleName,y=log2FoldChange,group=wormbaseID))+
  geom_line(alpha=0.1) +ggtitle("TOP2 trajectories")

p<-ggpubr::ggarrange(p1,p2,nrow=2)

ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"MoraoTOP1TOP2trajectories.pdf"),plot=p, device="pdf",width=19,height=29,units="cm")

########################-
## pairwise correlation of Bolaji and Morao data ------
########################-
moraoTop1<-fileList[grepl("top1top2morao",fileList$filePath) & grepl("TOP1",fileList$sampleName),]
moraoTop2<-fileList[grepl("top1top2morao",fileList$filePath) & grepl("TOP2",fileList$sampleName),]
bolajiTop1<-fileList[grepl("Top1AUX",fileList$filePath),]
bolajiTop2<-fileList[grepl("Top2AUX",fileList$filePath),]

# Plot pairwise correlations
corrPlotFunc<-function(data,x,y,contrastName,minScale,maxScale){
  facet_labels<-pairwiseTable %>% group_by(contrastName) %>%
    summarise(count=n(),label=paste0("n=",count))
    ggplot(pairwiseTable,aes(x=log2FoldChange.x,y=log2FoldChange.y))+
    geom_point(size=1,alpha=0.4) + facet_wrap(.~contrastName)+
    coord_cartesian(xlim=c(minScale,maxScale), ylim=c(minScale,maxScale)) +
    geom_smooth(method=lm,se=F,fullrange=T, size=0.7,show.legend = T) +
    geom_hline(yintercept=0,lty=3,col="grey70",) +
    geom_vline(xintercept=0,lty=3,col="grey70") +
    ggpubr::stat_cor(aes(label = ..r.label..), method="spearman",
                     cor.coef.name = c("Rs"), output.type = "text",
                     show.legend=F,size=3,
                     label.x=minScale+0.1, label.y=c(maxScale-0.1,maxScale-0.5)) +
    ggpubr::stat_cor(aes(label = ..r.label..), method="pearson",
                     cor.coef.name = c("Rp"), output.type = "text",
                     show.legend=F,size=3,
                     label.x=minScale+0.1, label.y=c(maxScale-0.6,maxScale-1.1),colour="red") +
    xlab(label=element_blank()) + ylab(label=element_blank()) +
    geom_text(data=facet_labels,aes(label=label),x=minScale+0.5,y=minScale+0.5,size=3,hjust=0.5)
}


# top1
rl<-getListOfResults(rbind(bolajiTop1,moraoTop1))
lapply(rl,dim)
res<-do.call(rbind,rl)

pairwiseTable<-NULL
for(grp in names(rl)[2:5]){
  tmp<-res[res$sampleName=="TOP1",]
  tmp<-left_join(tmp,res[res$sampleName==grp,],by="wormbaseID")
  tmp$contrastName<-paste0("TOP1(Bolaji) vs ",grp,"(Morao)")
  if(is.null(pairwiseTable)){
    pairwiseTable<-tmp
  } else {
    pairwiseTable<-rbind(pairwiseTable,tmp)
  }
}

pairwiseTable$contrastName<-factor(pairwiseTable$contrastName,levels=unique(pairwiseTable$contrastName))

minScale<- -3
maxScale<- 3
p1<-corrPlotFunc(pairwiseTable,log2FoldChange.x,log2FoldChange.y,contrastName,minScale,maxScale)

#top2
rl<-getListOfResults(rbind(bolajiTop2,moraoTop2))
lapply(rl,dim)
res<-do.call(rbind,rl)

pairwiseTable<-NULL
for(grp in names(rl)[2:5]){
  tmp<-res[res$sampleName=="TOP2",]
  tmp<-left_join(tmp,res[res$sampleName==grp,],by="wormbaseID")
  tmp$contrastName<-paste0("TOP2(Bolaji) vs ",grp,"(Morao)")
  if(is.null(pairwiseTable)){
    pairwiseTable<-tmp
  } else {
    pairwiseTable<-rbind(pairwiseTable,tmp)
  }
}

pairwiseTable$contrastName<-factor(pairwiseTable$contrastName,levels=unique(pairwiseTable$contrastName))

minScale<- -3
maxScale<- 3
p2<-corrPlotFunc(pairwiseTable,log2FoldChange.x,log2FoldChange.y,contrastName,minScale,maxScale)
p2

p<-ggpubr::ggarrange(p1,p2,nrow=2)

ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"corr_Bolaji-Morao.png"),plot=p, device="png",width=19,height=29,units="cm")

########################-
## pairwise correlation of coh-1 with Bolaji and Morao data ------
########################-
moraoTop1<-fileList[grepl("top1top2morao",fileList$filePath) & grepl("TOP1",fileList$sampleName),]
moraoTop2<-fileList[grepl("top1top2morao",fileList$filePath) & grepl("TOP2",fileList$sampleName),]
bolajiTop1<-fileList[grepl("Top1AUX",fileList$filePath),]
bolajiTop2<-fileList[grepl("Top2AUX",fileList$filePath),]
coh1<-fileList[grepl("COH1cs",fileList$sampleName),]


rl<-getListOfResults(rbind(coh1,bolajiTop1,moraoTop1[2:4,]))
lapply(rl,dim)
res<-do.call(rbind,rl)

pairwiseTable<-NULL
for(grp in names(rl)[2:5]){
  tmp<-res[res$sampleName=="COH1cs",]
  tmp<-inner_join(tmp,res[res$sampleName==grp,],by="wormbaseID")
  tmp$contrastName<-paste0("COH1cs vs ",grp)
  if(is.null(pairwiseTable)){
    pairwiseTable<-tmp
  } else {
    pairwiseTable<-rbind(pairwiseTable,tmp)
  }
}

pairwiseTable$contrastName<-factor(pairwiseTable$contrastName,levels=unique(pairwiseTable$contrastName))

minScale<- -3
maxScale<- 3
p1<-corrPlotFunc(pairwiseTable,log2FoldChange.x,log2FoldChange.y,contrastName,minScale,maxScale)
p1

#top2
rl<-getListOfResults(rbind(coh1,bolajiTop2,moraoTop2[2:4,]))
lapply(rl,dim)
res<-do.call(rbind,rl)

pairwiseTable<-NULL
for(grp in names(rl)[2:5]){
  tmp<-res[res$sampleName=="COH1cs",]
  tmp<-inner_join(tmp,res[res$sampleName==grp,],by="wormbaseID")
  tmp$contrastName<-paste0("COH1cs vs ",grp)
  if(is.null(pairwiseTable)){
    pairwiseTable<-tmp
  } else {
    pairwiseTable<-rbind(pairwiseTable,tmp)
  }
}

pairwiseTable$contrastName<-factor(pairwiseTable$contrastName,levels=unique(pairwiseTable$contrastName))

p2<-corrPlotFunc(pairwiseTable,log2FoldChange.x,log2FoldChange.y,contrastName,minScale,maxScale)
#p2<-p2+ggtitle("Top2 (all genes)")
p2

p<-ggpubr::ggarrange(p1,p2,nrow=2)
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"corr_COH1vsBolaji-Morao.png"),plot=p, device="png",width=19,height=29,units="cm")

########################-
## pairwise correlation of coh-1 significant gene with Bolaji and Morao data ------
########################-
moraoTop1<-fileList[grepl("top1top2morao",fileList$filePath) & grepl("TOP1",fileList$sampleName),]
moraoTop2<-fileList[grepl("top1top2morao",fileList$filePath) & grepl("TOP2",fileList$sampleName),]
bolajiTop1<-fileList[grepl("Top1AUX",fileList$filePath),]
bolajiTop2<-fileList[grepl("Top2AUX",fileList$filePath),]
coh1<-fileList[grepl("COH1cs",fileList$sampleName),]


rl<-getListOfResults(rbind(coh1,bolajiTop1,moraoTop1[2:4,]))
lapply(rl,dim)
res<-do.call(rbind,rl)

toKeep<-res$wormbaseID[res$sampleName=="COH1cs" & res$padj<0.05]
res<-res[res$wormbaseID %in% toKeep,]


pairwiseTable<-NULL
for(grp in names(rl)[2:5]){
  tmp<-res[res$sampleName=="COH1cs",]
  tmp<-inner_join(tmp,res[res$sampleName==grp,],by="wormbaseID")
  tmp$contrastName<-paste0("COH1cs vs ",grp)
  if(is.null(pairwiseTable)){
    pairwiseTable<-tmp
  } else {
    pairwiseTable<-rbind(pairwiseTable,tmp)
  }
}

pairwiseTable$contrastName<-factor(pairwiseTable$contrastName,levels=unique(pairwiseTable$contrastName))

minScale<- -3
maxScale<- 3
p1<-corrPlotFunc(pairwiseTable,log2FoldChange.x,log2FoldChange.y,contrastName,minScale,maxScale)
p1

#top2
rl<-getListOfResults(rbind(coh1,bolajiTop2,moraoTop2[2:4,]))
lapply(rl,dim)
res<-do.call(rbind,rl)

toKeep<-res$wormbaseID[res$sampleName=="COH1cs" & res$padj<0.05]
res<-res[res$wormbaseID %in% toKeep,]

pairwiseTable<-NULL
for(grp in names(rl)[2:5]){
  tmp<-res[res$sampleName=="COH1cs",]
  tmp<-inner_join(tmp,res[res$sampleName==grp,],by="wormbaseID")
  tmp$contrastName<-paste0("COH1cs vs ",grp)
  if(is.null(pairwiseTable)){
    pairwiseTable<-tmp
  } else {
    pairwiseTable<-rbind(pairwiseTable,tmp)
  }
}

pairwiseTable$contrastName<-factor(pairwiseTable$contrastName,levels=unique(pairwiseTable$contrastName))

p2<-corrPlotFunc(pairwiseTable,log2FoldChange.x,log2FoldChange.y,contrastName,minScale,maxScale)
#p2<-p2+ggtitle("Top2 (all genes)")
p2

p<-ggpubr::ggarrange(p1,p2,nrow=2)
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"corr_COH1sigvsBolaji-Morao.png"),plot=p, device="png",width=19,height=29,units="cm")




########################-
## pairwise correlation of ama-1 with other data ------
########################-
moraoTop1<-fileList[grepl("top1top2morao",fileList$filePath) & grepl("TOP1",fileList$sampleName),]
moraoTop2<-fileList[grepl("top1top2morao",fileList$filePath) & grepl("TOP2",fileList$sampleName),]
bolajiTop1<-fileList[grepl("Top1AUX",fileList$filePath),]
bolajiTop2<-fileList[grepl("Top2AUX",fileList$filePath),]
coh1<-fileList[grepl("COH1cs",fileList$sampleName),]
ama<-fileList[grepl("AMA",fileList$sampleName),]


rl<-getListOfResults(rbind(ama,moraoTop1[2:4,],moraoTop2[2:4,],bolajiTop1,bolajiTop2,coh1))
lapply(rl,dim)
res<-do.call(rbind,rl)

pairwiseTable<-NULL
for(grp in names(rl)[2:length(rl)]){
  tmp<-res[res$sampleName=="AMA",]
  tmp<-inner_join(tmp,res[res$sampleName==grp,],by="wormbaseID")
  tmp$contrastName<-paste0("AMA vs ",grp)
  if(is.null(pairwiseTable)){
    pairwiseTable<-tmp
  } else {
    pairwiseTable<-rbind(pairwiseTable,tmp)
  }
}

pairwiseTable$contrastName<-factor(pairwiseTable$contrastName,levels=unique(pairwiseTable$contrastName))

minScale<- -5
maxScale<- 5
p1<-corrPlotFunc(pairwiseTable,log2FoldChange.x,log2FoldChange.y,contrastName,minScale,maxScale)
p1

ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"corr_AMAvsCOH1-Bolaji-Morao.png"),plot=p1, device="png",width=19,height=19,units="cm")



###################-
## Compare coh-1cs analysis
###################-

mdas<-readRDS("/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/2021_RNAseq_MDas/rds/p0.05_lfc0.5_filtChrAX/filtChrAX_X.wt.wt.0mM_coh1cs_vs_wt_DESeq2_fullResults_p0.05.rds")

salmon<-readRDS(file=fileList[fileList$sampleName=="COH1cs","filePath"])

sigGenesUp<-list()
sigGenesUp[["mdas"]]<-mdas$wormbaseID[!is.na(mdas$padj) & mdas$padj<0.05 & mdas$log2FoldChange>0]
sigGenesUp[["fount"]]<-salmon$wormbaseID[!is.na(salmon$padj) & salmon$padj<0.05 & salmon$log2FoldChange>0]


fit<-euler(sigGenesUp)
p1<-plot(fit, quantities=list(type="counts"),
         main=list(label=paste0("Up regulated genes: ", "padj<",padjVal,"\n",
                                paste(lapply(row.names(fit$ellipses), function(x){
                                  paste(x, sum(fit$original.values[grep(x,names(fit$original.values))]))
                                }), collapse="  ")), fontsize=8))
p1


sigGenesDown<-list()
sigGenesDown[["mdas"]]<-mdas$wormbaseID[!is.na(mdas$padj) & mdas$padj<0.05 & mdas$log2FoldChange<0]
sigGenesDown[["fount"]]<-salmon$wormbaseID[!is.na(salmon$padj) & salmon$padj<0.05 & salmon$log2FoldChange<0]

fit<-euler(sigGenesDown)
p2<-plot(fit, quantities=list(type="counts"),
         main=list(label=paste0("Down regulated genes: ", "padj<",padjVal,"\n",
                                paste(lapply(row.names(fit$ellipses), function(x){
                                  paste(x, sum(fit$original.values[grep(x,names(fit$original.values))]))
                                }), collapse="  ")), fontsize=8))
p2

# p<-ggarrange(p1,p2,nrow=2)
# ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"venn_coh1_mdasVfount.pdf"),plot=p, device="pdf",width=19,height=29,units="cm")

c1<-inner_join(data.frame(mdas),data.frame(salmon),by="wormbaseID")
dim(c1)
minScale=-3
maxScale=3
p3<-ggplot(c1,aes(x=log2FoldChange.x,y=log2FoldChange.y)) +
  geom_point(size=1,alpha=0.4) +
  coord_cartesian(xlim=c(minScale,maxScale), ylim=c(minScale,maxScale)) +
  geom_smooth(method=lm,se=F,fullrange=T, size=0.7,show.legend = T) +
  geom_hline(yintercept=0,lty=3,col="grey70",) +
  geom_vline(xintercept=0,lty=3,col="grey70") +
  ggpubr::stat_cor(aes(label = ..r.label..), method="spearman",
                   cor.coef.name = c("Rs"), output.type = "text",
                   show.legend=F,size=3,
                   label.x=minScale+0.1, label.y=c(maxScale-0.1,maxScale-0.5)) +
  ggpubr::stat_cor(aes(label = ..r.label..), method="pearson",
                   cor.coef.name = c("Rp"), output.type = "text",
                   show.legend=F,size=3,
                   label.x=minScale+0.1, label.y=c(maxScale-0.6,maxScale-1.1),colour="red") +
  xlab(label="mdas LFC") + ylab("fount LFC")

p<-ggarrange(p1,p2,p3,nrow=2,ncol=2)
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"venn_coh1_mdasVfount.pdf"),plot=p, device="pdf",width=19,height=19,units="cm")
