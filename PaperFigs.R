library(GenomicRanges)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(plyranges)
library(rtracklayer)
library(dplyr)
library(ggtext)

outPath="."
scriptPath="."
source(paste0(scriptPath,"/functions.R"))

theme_set(
  theme_bw()+
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
)

scriptName="PaperFigs"
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

fileList<-fileList[fileList$sampleName=="COH1cs",]


fountains<-readRDS("/Users/semple/Documents/MeisterLab/otherPeopleProjects/fountains/detected_fountains_equalQ.RDS")
fountains$fountainName<-paste0("fount",1:length(fountains))

nonFount<-gaps(fountains)
nonFount<-resize(nonFount,width=2000,fix="center")

#' Pairwise faceted violin/boxplot with stats
#'
#' Takes a df with columns sampleName, regionType and log2FoldChange and
#' plots violin/boxplot for the two regionTyps faceted by sampleName.
#' t_test is performed on each pair and significance indicated. Total number
#' of regions used also indicated in grey on bottom.
#' @param df with sampleName, regionType and log2FoldChange columns
#' @param gr regions used to plot (needed to calculate width)
#' @param yvar Name of variable to plot on y axis (default is Log2FoldChange)
#' @param ymin bottom of y scale used by coord_cartesian
#' @param ymax top of y scale used by coord_cartesian
#' @param facet_by Name of variable to facet by (default is sampleName). Can
#' be set to NULL to eliminate facetting
#' @param geneSet Name of gene set used in plot title e.g. "all" or "significant"
#' @return plot
#' @export
pairwiseBoxPlotFunc<-function(df,gr,yvar="log2FoldChange",ymin=-1,ymax=1,facet_by="sampleName",geneSet="all"){
  df$sampleName<-factor(df$sampleName, levels=unique(df$sampleName))
  df$regionType<-factor(df$regionType)

  stat_df<-df %>% rstatix::group_by(sampleName,drop=T) %>%
    rstatix::mutate(column=get(yvar)) %>%
    rstatix::t_test(column~regionType) %>%
    rstatix::adjust_pvalue(method="BH") %>%
    rstatix::add_significance("p.adj") %>%
    rstatix::add_xy_position(x="regionType") %>%
    rstatix::mutate(y.position=0.9*ymax)
  if(!is.null(facet_by)){
    label_df<-df%>% dplyr::group_by(.data[[facet_by]],regionType) %>% dplyr::summarise(count=dplyr::n())
  } else {
    label_df<-df%>% dplyr::group_by(regionType) %>% dplyr::summarise(count=dplyr::n())
  }
  print(stat_df,width=Inf)

  p<-ggplot(df,aes(x=regionType,y=get(yvar))) +
    geom_violin(width=0.9,mapping=aes(fill=regionType))+
    geom_boxplot(width=0.2,mapping=aes(fill=regionType),outlier.shape=NA)+
    scale_fill_manual(values=c("pink","lightblue"))+
    coord_cartesian(ylim=c(ymin,ymax)) + ylab(yvar) +
    stat_pvalue_manual(stat_df, label = "p.adj", remove.bracket=F,hide.ns = T,
                       color="purple", bracket.size=0.5,bracket.short=0.1,tip.length=0.001)+
    geom_text(data=label_df,aes(label=count),y=ymin,size=3,colour="darkgrey") +
    ggtitle(paste0(yvar," of ",geneSet," genes in fountains (",width(gr[1])/1000,"kb)"))
  if(ymin<0){
    p<-p+geom_hline(yintercept=0,colour="red",linetype=2)
  }
  if(!is.null(facet_by)){
    p<-p+facet_wrap(.~get(facet_by),nrow=3)
  }
  return(p)
}

# all genes 2kb
stl<-resultsByGRoverlap(fileList,fountains)
fountRes<-do.call(rbind,lapply(stl,as.data.frame))
fountRes$regionType<-"fountain"

stl<-resultsByGRoverlap(fileList,nonFount)
nonFountRes<-do.call(rbind,lapply(stl,as.data.frame))
nonFountRes$regionType<-"between"

df<-rbind(fountRes,nonFountRes)

p<-pairwiseBoxPlotFunc(df,fountains,ymin=-1.5,ymax=1.5)
p
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"LFCgenesInFountains",width(fountains[1])/1000,"kb_boxplots.pdf"),plot=p, device="pdf",width=12,height=12,units="cm")

# all genes 6kb
fountains_6kb<-resize(fountains,width=6000,fix="center")
nonFount_6kb<-resize(nonFount,width=6000,fix="center")

stl<-resultsByGRoverlap(fileList,fountains_6kb)
fountRes<-do.call(rbind,lapply(stl,as.data.frame))
fountRes$regionType<-"fountain"

stl<-resultsByGRoverlap(fileList,nonFount_6kb)
nonFountRes<-do.call(rbind,lapply(stl,as.data.frame))
nonFountRes$regionType<-"between"

df<-rbind(fountRes,nonFountRes)

p<-pairwiseBoxPlotFunc(df,fountains_6kb,ymin=-1.5,ymax=1.5)
p
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"LFCgenesInFountains",width(fountains_6kb[1])/1000,"kb_boxplots.pdf"),plot=p, device="pdf",width=12,height=12,units="cm")

## sig genes 2kb
stl<-resultsByGRoverlap(fileList,fountains,padjVal=0.05)
stl<-stl[!sapply(stl,length)<10]
fountRes<-do.call(rbind,lapply(stl,as.data.frame))
fountRes$regionType<-"fountain"

stl<-resultsByGRoverlap(fileList,nonFount,padjVal=0.05)
stl<-stl[!sapply(stl,length)<10]
nonFountRes<-do.call(rbind,lapply(stl,as.data.frame))
nonFountRes$regionType<-"between"

df<-rbind(fountRes,nonFountRes)

p<-pairwiseBoxPlotFunc(df,fountains,ymin=-1.5,ymax=1.5)
p
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"LFCsigGenesInFountains",width(fountains[1])/1000,"kb_boxplots.pdf"),plot=p, device="pdf",width=12,height=12,units="cm")


## sig genes 6kb
stl<-resultsByGRoverlap(fileList,fountains_6kb,padjVal=0.05)
stl<-stl[!sapply(stl,length)<10]
fountRes<-do.call(rbind,lapply(stl,as.data.frame))
fountRes$regionType<-"fountain"

stl<-resultsByGRoverlap(fileList,nonFount_6kb,padjVal=0.05)
stl<-stl[!sapply(stl,length)<10]
nonFountRes<-do.call(rbind,lapply(stl,as.data.frame))
nonFountRes$regionType<-"between"

df<-rbind(fountRes,nonFountRes)

p<-pairwiseBoxPlotFunc(df,fountains_6kb,ymin=-1.5,ymax=1.5)
p
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"LFCsigGenesInFountains",width(fountains_6kb[1])/1000,"kb_boxplots.pdf"),plot=p, device="pdf",width=12,height=12,units="cm")



#######################-
## get gene specific data -----
#######################-
gtf<-rtracklayer::import("/Users/semple/Documents/MeisterLab/GenomeVer/WS285/c_elegans.PRJNA13758.WS285.annotations.gtf")

numTranscripts<-as.data.frame(gtf) %>% dplyr::group_by(gene_id) %>%
  dplyr::summarise(numTranscripts=dplyr::n_distinct(transcript_id))

lengthTranscripts<-as.data.frame(gtf) %>% dplyr::group_by(gene_id) %>%
  filter(type=="transcript") %>%
  dplyr::summarise(maxTxLength=max(width),minTxLength=min(width),avrTxLength=mean(width))

numStartEnd<-as.data.frame(gtf) %>% dplyr::group_by(gene_id) %>%
  filter(type=="transcript")%>%
  dplyr::summarise(numStarts=dplyr::n_distinct(start),
                   numEnds=dplyr::n_distinct(end))

numExons<-gtf %>%  plyranges::filter(type=="exon") %>% unique() %>% plyranges::group_by(gene_id) %>%
  plyranges::summarise(numExons=plyranges::n())

geneData<-dplyr::left_join(numTranscripts,lengthTranscripts,by="gene_id")
geneData<-dplyr::left_join(geneData,numStartEnd,by="gene_id")
geneData<-dplyr::left_join(geneData,as.data.frame(numExons),by="gene_id")
geneData$gene_id<-gsub("Gene:","",geneData$gene_id)
names(geneData)[names(geneData)=="gene_id"]<-"wormbaseID"
geneData

rl<-getListOfResults(fileList,padjVal=0.05)
lapply(rl,dim)
rlg<-lapply(rl,dplyr::left_join,geneData,by="wormbaseID")
lapply(rlg,dim)
res<-do.call(rbind,rlg)

plotFeatureFunc<-function(res,gFeat,ymin=0,ymax=5){
  label_df<-res%>% dplyr::group_by(sampleName) %>% dplyr::summarise(count=n())
  p<-ggplot(res,aes(x=sampleName,y=get(gFeat))) +
    geom_violin(width=0.9,fill="lightblue",alpha=0.2,colour="lightblue")+
    ylab(gFeat) +
    ggtitle(paste0(gFeat," of significant genes"))

  if(ymax>50){
    print(paste0(ymax))
    ymax=quantile(res[,gFeat],0.99)
    p<-p+geom_boxplot(width=0.2,outlier.shape=NA,fill="lightblue")+
      geom_text(data=label_df,aes(label=count),y=ymax,size=3,colour="darkgrey") +
      coord_cartesian(ylim = c(ymin,ymax))
  } else {
    ymax=max(res[,gFeat])
    p<-p+geom_violin(width=0.9,fill="lightblue",alpha=1,colour="black",linewidth=0.2) +
      geom_text(data=label_df,aes(label=count),y=ymax,size=3,colour="darkgrey")
  }
  return(p)
}

# gFeatures<-c("numTranscripts","maxTxLength","minTxLength","avrTxLength", "numStarts",
#              "numEnds","numExons")
# plotList<-list()
# for(gFeat in gFeatures){
#   res$sampleName<-factor(res$sampleName,levels=unique(res$sampleName))
#   ymax=quantile(res[,gFeat],0.99)
#   plotList[[gFeat]]<-plotFeatureFunc(res,gFeat,ymin=0,ymax=ymax)
# }
#
# ml<-marrangeGrob(plotList,nrow=4,ncol=2)
# ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"LFCgeneFeatures.pdf"),plot=ml, device="pdf",width=19,height=29,units="cm")


##########################-
## gene features in and out of fountains
##########################-

stl<-resultsByGRoverlap(fileList,fountains)
stl<-lapply(stl,as.data.frame)
stl<-lapply(stl,dplyr::left_join,geneData,by="wormbaseID")
lapply(stl,dim)
fountRes<-do.call(rbind,stl)
fountRes$regionType<-"fountain"

stl<-resultsByGRoverlap(fileList,nonFount)
stl<-lapply(stl,as.data.frame)
stl<-lapply(stl,dplyr::left_join,geneData,by="wormbaseID")
lapply(stl,dim)
nonFountRes<-do.call(rbind,stl)
nonFountRes$regionType<-"between"

df<-rbind(fountRes,nonFountRes)


gFeatures<-c("numTranscripts","maxTxLength","avrTxLength", "numStarts",
             "numEnds","numExons")
plotList<-list()
for(gFeat in gFeatures){
  res$sampleName<-factor(res$sampleName,levels=unique(res$sampleName))
  ymax=quantile(res[,gFeat],0.999)
  ymin=0
  plotList[[gFeat]]<-pairwiseBoxPlotFunc(df[df$sampleName=="COH1cs",],fountains,yvar=gFeat, ymin=ymin,ymax=ymax,facet_by=NULL)
}


ml<-marrangeGrob(plotList,nrow=2,ncol=3)
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"LFCgeneFeatures_fountains.pdf"),plot=ml, device="pdf",width=24,height=19,units="cm")


##########################-
## gene features in and out of fountains - significant genes
##########################-

stl<-resultsByGRoverlap(fileList,fountains,padjVal=0.05)
stl<-stl[!sapply(stl,length)<10]
stl<-lapply(stl,as.data.frame)
stl<-lapply(stl,dplyr::left_join,geneData,by="wormbaseID")
lapply(stl,dim)
fountRes<-do.call(rbind,stl)
fountRes$regionType<-"fountain"

stl<-resultsByGRoverlap(fileList,nonFount,padjVal=0.05)
stl<-stl[!sapply(stl,length)<10]
stl<-lapply(stl,as.data.frame)
stl<-lapply(stl,dplyr::left_join,geneData,by="wormbaseID")
lapply(stl,dim)
nonFountRes<-do.call(rbind,stl)
nonFountRes$regionType<-"between"

df<-rbind(fountRes,nonFountRes)


gFeatures<-c("numTranscripts","maxTxLength","avrTxLength", "numStarts",
             "numEnds","numExons")
plotList<-list()
for(gFeat in gFeatures){
  res$sampleName<-factor(res$sampleName,levels=unique(res$sampleName))
  ymax=quantile(res[,gFeat],0.999)
  ymin=0
  plotList[[gFeat]]<-pairwiseBoxPlotFunc(df,fountains,yvar=gFeat, ymin=ymin,
                                         ymax=ymax,facet_by="sampleName",
                                         geneSet="significant")
}


ml<-marrangeGrob(plotList,nrow=2,ncol=3)
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"LFCsigGeneFeatures_fountains.pdf"),plot=ml, device="pdf",width=24,height=19,units="cm")

