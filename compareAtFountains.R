library(GenomicRanges)
library(ggplot2)
library(ggpubr)
library(rstatix)

outPath="."
scriptPath="."
source(paste0(scriptPath,"/functions.R"))

theme_set(
  theme_bw()+
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
)

scriptName="compareAtFountains"
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

fountains<-readRDS("/Users/semple/Documents/MeisterLab/otherPeopleProjects/fountains/detected_fountains_equalQ.RDS")
fountains$fountainName<-paste0("fount",1:length(fountains))

nonFount<-gaps(fountains)
nonFount<-resize(nonFount,width=2000,fix="center")

stl<-resultsByGRoverlap(fileList,fountains)
fountRes<-do.call(rbind,lapply(stl,as.data.frame))
fountRes$regionType<-"fountain"

stl<-resultsByGRoverlap(fileList,nonFount)
nonFountRes<-do.call(rbind,lapply(stl,as.data.frame))
nonFountRes$regionType<-"between"

df<-rbind(fountRes,nonFountRes)

#' Pairwise faceted violin/boxplot with stats
#'
#' Takes a df with columns sampleName, regionType and log2FoldChange and
#' plots violin/boxplot for the two regionTyps faceted by sampleName.
#' t_test is performed on each pair and significance indicated. Total number
#' of regions used also indicated in grey on bottom.
#' @param df with sampleName, regionType and log2FoldChange columns
#' @param gr regions used to plot (needed to calculate width)
#' @return plot
#' @export
pairwiseBoxPlotFunc<-function(df,gr,ymin=-1,ymax=1){
  df$sampleName<-factor(df$sampleName, levels=unique(df$sampleName))
  df$regionType<-factor(df$regionType)

  stat_df<-df %>% rstatix::group_by(sampleName,drop=T) %>%
    rstatix::t_test(log2FoldChange~regionType) %>%
    rstatix::adjust_pvalue(method="BH") %>%
    rstatix::add_significance("p.adj") %>%
    rstatix::add_xy_position(x="regionType") %>%
    mutate(y.position=0.9*ymax)

  label_df<-df%>% group_by(sampleName,regionType) %>% summarise(count=n())

  print(stat_df,width=Inf)

  p<-ggplot(df,aes(x=regionType,y=log2FoldChange)) +
    geom_violin(width=0.9,mapping=aes(fill=regionType))+
    geom_boxplot(width=0.2,mapping=aes(fill=regionType),outlier.shape=NA)+
    scale_fill_manual(values=c("pink","lightblue"))+
    coord_cartesian(ylim=c(ymin,ymax))+
    facet_wrap(.~sampleName,nrow=3) +
    geom_hline(yintercept=0,colour="red",linetype=2)+
    stat_pvalue_manual(stat_df, label = "p.adj.signif", remove.bracket=F,hide.ns = T,
                       color="purple", bracket.size=0.5,bracket.short=0.1,tip.length=0.001)+
    geom_text(data=label_df,aes(label=count),y=ymin,size=3,colour="darkgrey") +
    ggtitle(paste0("LFC of all genes in fountains (",width(gr[1])/1000,"kb)"))
  return(p)
}

p<-pairwiseBoxPlotFunc(df,fountains)
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"LFCgenesInFountains",width(fountains[1])/1000,"kb_boxplots.pdf"),plot=p, device="pdf",width=29,height=19,units="cm")

fountains_6kb<-resize(fountains,width=6000,fix="center")
nonFount_6kb<-resize(nonFount,width=6000,fix="center")

stl<-resultsByGRoverlap(fileList,fountains_6kb)
fountRes<-do.call(rbind,lapply(stl,as.data.frame))
fountRes$regionType<-"fountain"

stl<-resultsByGRoverlap(fileList,nonFount_6kb)
nonFountRes<-do.call(rbind,lapply(stl,as.data.frame))
nonFountRes$regionType<-"between"

df<-rbind(fountRes,nonFountRes)

p<-pairwiseBoxPlotFunc(df,fountains_6kb)
p
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"LFCgenesInFountains",width(fountains_6kb[1])/1000,"kb_boxplots.pdf"),plot=p, device="pdf",width=29,height=19,units="cm")

## sig genes 2kb
stl<-resultsByGRoverlap(fileList,fountains,padjVal=0.05)
stl<-stl[-which(lapply(stl,length)<10)]
fountRes<-do.call(rbind,lapply(stl,as.data.frame))
fountRes$regionType<-"fountain"

stl<-resultsByGRoverlap(fileList,nonFount,padjVal=0.05)
stl<-stl[-which(lapply(stl,length)<10)]
nonFountRes<-do.call(rbind,lapply(stl,as.data.frame))
nonFountRes$regionType<-"between"

df<-rbind(fountRes,nonFountRes)

p<-pairwiseBoxPlotFunc(df,fountains,ymin=-3,ymax=3)
p
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"LFCsigGenesInFountains",width(fountains[1])/1000,"kb_boxplots.pdf"),plot=p, device="pdf",width=29,height=19,units="cm")


## sig genes 2kb
stl<-resultsByGRoverlap(fileList,fountains_6kb,padjVal=0.05)
stl<-stl[-which(lapply(stl,length)<10)]
fountRes<-do.call(rbind,lapply(stl,as.data.frame))
fountRes$regionType<-"fountain"

stl<-resultsByGRoverlap(fileList,nonFount_6kb,padjVal=0.05)
stl<-stl[-which(lapply(stl,length)<10)]
nonFountRes<-do.call(rbind,lapply(stl,as.data.frame))
nonFountRes$regionType<-"between"

df<-rbind(fountRes,nonFountRes)

p<-pairwiseBoxPlotFunc(df,fountains_6kb,ymin=-3,ymax=3)
p
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"LFCsigGenesInFountains",width(fountains_6kb[1])/1000,"kb_boxplots.pdf"),plot=p, device="pdf",width=29,height=19,units="cm")
