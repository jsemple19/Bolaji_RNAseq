library(GenomicRanges)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(plyranges)
library(rtracklayer)
library(dplyr)

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
    stat_pvalue_manual(stat_df, label = "p.adj.signif", remove.bracket=F,hide.ns = T,
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


stl<-resultsByGRoverlap(fileList,fountains)
fountRes<-do.call(rbind,lapply(stl,as.data.frame))
fountRes$regionType<-"fountain"

stl<-resultsByGRoverlap(fileList,nonFount)
nonFountRes<-do.call(rbind,lapply(stl,as.data.frame))
nonFountRes$regionType<-"between"

df<-rbind(fountRes,nonFountRes)

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





######################-
## binning by distance------
######################-
binByDistance<-function(distance,maxDist=50000,binSize=2000){
  binnedDistance<-cut(distance,breaks=seq(0,maxDist,by=binSize),include.lowest=F)
  levels(binnedDistance)<-seq(binSize,maxDist,by=binSize)/1000
  levels(binnedDistance)<-c(seq(binSize,maxDist,by=binSize)/1000,paste0(">",(maxDist/1000)),0)
  binnedDistance[is.na(binnedDistance)]<-paste0(">",(maxDist/1000))
  binnedDistance[distance==0]<-0
  binnedDistance<-relevel(binnedDistance, ref="0")
  table(binnedDistance)
  return(binnedDistance)
}



#' Plot Log2 fold change of genes binned by the distance of their start from fountains
#'
#' Using a DESeq2 results table and GenomicRanges object for fountains, find
#' the distance of each gene to nearest fountain. Bin the genes by their distance
#' from the nearest fountain and plot boxplot of LFC values.
#' @param res DESeq2 results object
#' @param fountains GRanges of fountains
#' @param sampleName String with name of sample being plotted
#' @param maxDist Maximal distance away from fountain to look at (in bp)
#' @param binSize Size of bins to use (in bp)
#' @param padjVal Null by default to look at all genes, but can be set to a p value threshold
#' to filter significant genes
#' @param ymin coord_cartesian limits of the plots (minimum of y scale)
#' @param ymax coord_cartesian limits of the plots (maximum of y scale)
#' @return plot
#' @export
plotBinFountainDistanceFunc<-function(res,fountains,sampleName,maxDist,binSize,padjVal=NULL,ymin=0,ymax=1){
  salmon<-res[!is.na(res$chr),]
  salmongr<-makeGRangesFromDataFrame(salmon,keep.extra.columns = T)

  salmongr<-sort(salmongr)
  salmongr<-resize(salmongr,width=1,fix="start")

  ol<-findOverlaps(salmongr,fountains,ignore.strand=T)
  dtn<-distanceToNearest(salmongr,fountains,ignore.strand=T)
  #hist(mcols(dtn)$distance)

  salmongr$distanceToFountain<-mcols(dtn)$distance
  mcols(salmongr)<-cbind(mcols(salmongr),mcols(fountains[subjectHits(dtn),]))
  df<-data.frame(salmongr)

  #smaller bins, and zoomed in
  df$binnedDistance<-binByDistance(df$distanceToFountain,maxDist=maxDist,binSize=binSize)
  table(df$binnedDistance)
  df$upVdown<-factor(ifelse(df$log2FoldChange>0,"up","down"),levels=c("up","down"))

  if(!is.null(padjVal)){
    df<-df %>% dplyr::filter(!is.na(df$padj), df$padj<padjVal)
    table(df$binnedDistance)

    countLabels<-df %>% dplyr::group_by(binnedDistance,upVdown) %>% dplyr::summarise(count=dplyr::n())
    countLabels_a<-df %>% dplyr::group_by(binnedDistance) %>% dplyr::summarise(count=dplyr::n())
  }
  countLabels<-df %>% dplyr::group_by(binnedDistance,upVdown) %>% dplyr::summarise(count=dplyr::n())
  countLabels_a<-df %>% dplyr::group_by(binnedDistance) %>% dplyr::summarise(count=dplyr::n())

  # boxplot
  p<-ggplot(df) +
    geom_boxplot(mapping=aes(x=binnedDistance,y=abs(log2FoldChange)), fill="lightblue",
                 outlier.colour="grey",notch=F, varwidth=F) +
    theme_bw() + coord_cartesian(ylim=c(ymin,ymax)) +
    facet_grid(rows=vars(upVdown)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
    geom_text(data=countLabels,aes(x=binnedDistance,label=count),y=ymin,size=3,colour="blue") +
    xlab("Distance from fountain bins (kb)") +ggtitle(paste0(sampleName,": ",binSize/1000,"kb bins"))
  return(p)
}


subset<-c("AMA", "COH1cs", "TOP1", "TOP2", "TOP1_1h", "TOP1_2h", "TOP2_1h", "TOP2_2h")
maxDist=10000
binSize=1000
plotList<-list()
for(grp in subset){
  salmon<-readRDS(file=fileList[fileList$sampleName==grp,"filePath"])
  plotList[[grp]]<-plotBinFountainDistanceFunc(salmon,fountains,grp,maxDist=maxDist,binSize=binSize,padjVal=NULL,ymin=0,ymax=1)
}

ml<-marrangeGrob(plotList,nrow=3,ncol=1)
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"LFCDistToFountains_bin", binSize/1000, "kb_", maxDist/1000, "kb_boxplots.pdf"),plot=ml, device="pdf",width=19,height=29,units="cm")


maxDist=10000
binSize=1000
plotList<-list()
for(grp in subset){
  salmon<-readRDS(file=fileList[fileList$sampleName==grp,"filePath"])
  plotList[[grp]]<-plotBinFountainDistanceFunc(salmon,fountains,grp,maxDist=maxDist,binSize=binSize,padjVal=padjVal,ymin=0,ymax=3)
}

ml<-marrangeGrob(plotList,nrow=3,ncol=1)
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"LFCsigDistToFountains_bin", binSize/1000, "kb_", maxDist/1000, "kb_boxplots.pdf"),plot=ml, device="pdf",width=19,height=29,units="cm")


##########-
## compare start vs end overlapping with fountains
##########-

stl<-resultsByGRoverlap(fileList[fileList$sampleName=="COH1cs",],fountains)
stl[["COH1cs_start"]]<-resultsByGRoverlap(fileList[fileList$sampleName=="COH1cs",],fountains,resizeTo="start")[[1]]
stl[["COH1cs_start"]]$sampleName<-"COH1cs_start"
stl[["COH1cs_center"]]<-resultsByGRoverlap(fileList[fileList$sampleName=="COH1cs",],fountains,resizeTo="center")[[1]]
stl[["COH1cs_center"]]$sampleName<-"COH1cs_center"
stl[["COH1cs_end"]]<-resultsByGRoverlap(fileList[fileList$sampleName=="COH1cs",],fountains,resizeTo="end")[[1]]
stl[["COH1cs_end"]]$sampleName<-"COH1cs_end"

fountRes<-do.call(rbind,lapply(stl,as.data.frame))
fountRes$regionType<-"fountain"

stl<-resultsByGRoverlap(fileList[fileList$sampleName=="COH1cs",],nonFount)
stl[["COH1cs_start"]]<-resultsByGRoverlap(fileList[fileList$sampleName=="COH1cs",],nonFount,resizeTo="start")[[1]]
stl[["COH1cs_start"]]$sampleName<-"COH1cs_start"
stl[["COH1cs_center"]]<-resultsByGRoverlap(fileList[fileList$sampleName=="COH1cs",],nonFount,resizeTo="center")[[1]]
stl[["COH1cs_center"]]$sampleName<-"COH1cs_center"
stl[["COH1cs_end"]]<-resultsByGRoverlap(fileList[fileList$sampleName=="COH1cs",],nonFount,resizeTo="end")[[1]]
stl[["COH1cs_end"]]$sampleName<-"COH1cs_end"
nonFountRes<-do.call(rbind,lapply(stl,as.data.frame))
nonFountRes$regionType<-"between"

df<-rbind(fountRes,nonFountRes)

p<-pairwiseBoxPlotFunc(df,fountains)
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"LFCgenesInFountains_startVend_",width(fountains[1])/1000,"kb_boxplots.pdf"),plot=p, device="pdf",width=29,height=19,units="cm")

# significant genes
stl<-resultsByGRoverlap(fileList[fileList$sampleName=="COH1cs",],fountains,padjVal=0.05)
stl[["COH1cs_start"]]<-resultsByGRoverlap(fileList[fileList$sampleName=="COH1cs",],fountains,padjVal=0.05,resizeTo="start")[[1]]
stl[["COH1cs_start"]]$sampleName<-"COH1cs_start"
stl[["COH1cs_center"]]<-resultsByGRoverlap(fileList[fileList$sampleName=="COH1cs",],fountains,padjVal=0.05,resizeTo="center")[[1]]
stl[["COH1cs_center"]]$sampleName<-"COH1cs_center"
stl[["COH1cs_end"]]<-resultsByGRoverlap(fileList[fileList$sampleName=="COH1cs",],fountains,padjVal=0.05,resizeTo="end")[[1]]
stl[["COH1cs_end"]]$sampleName<-"COH1cs_end"

fountRes<-do.call(rbind,lapply(stl,as.data.frame))
fountRes$regionType<-"fountain"

stl<-resultsByGRoverlap(fileList[fileList$sampleName=="COH1cs",],nonFount,padjVal=0.05)
stl[["COH1cs_start"]]<-resultsByGRoverlap(fileList[fileList$sampleName=="COH1cs",],nonFount,padjVal=0.05,resizeTo="start")[[1]]
stl[["COH1cs_start"]]$sampleName<-"COH1cs_start"
stl[["COH1cs_center"]]<-resultsByGRoverlap(fileList[fileList$sampleName=="COH1cs",],nonFount,padjVal=0.05,resizeTo="center")[[1]]
stl[["COH1cs_center"]]$sampleName<-"COH1cs_center"
stl[["COH1cs_end"]]<-resultsByGRoverlap(fileList[fileList$sampleName=="COH1cs",],nonFount,padjVal=0.05,resizeTo="end")[[1]]
stl[["COH1cs_end"]]$sampleName<-"COH1cs_end"
nonFountRes<-do.call(rbind,lapply(stl,as.data.frame))
nonFountRes$regionType<-"between"

df<-rbind(fountRes,nonFountRes)
table(df$sampleName)

p<-pairwiseBoxPlotFunc(df,fountains)
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"LFCsigGenesInFountains_startVend_",width(fountains[1])/1000,"kb_boxplots.pdf"),plot=p, device="pdf",width=29,height=19,units="cm")

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

rl<-getListOfResults(fileList[-2,],padjVal=0.05)
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

gFeatures<-c("numTranscripts","maxTxLength","minTxLength","avrTxLength", "numStarts",
             "numEnds","numExons")
plotList<-list()
for(gFeat in gFeatures){
  res$sampleName<-factor(res$sampleName,levels=unique(res$sampleName))
  ymax=quantile(res[,gFeat],0.99)
  plotList[[gFeat]]<-plotFeatureFunc(res,gFeat,ymin=0,ymax=ymax)
}

ml<-marrangeGrob(plotList,nrow=4,ncol=1)
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"LFCgeneFeatures.pdf"),plot=ml, device="pdf",width=19,height=29,units="cm")

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


gFeatures<-c("numTranscripts","maxTxLength","minTxLength","avrTxLength", "numStarts",
             "numEnds","numExons")
plotList<-list()
for(gFeat in gFeatures){
  res$sampleName<-factor(res$sampleName,levels=unique(res$sampleName))
  ymax=quantile(res[,gFeat],0.99)
  ymin=0
  plotList[[gFeat]]<-pairwiseBoxPlotFunc(df[df$sampleName=="COH1cs",],fountains,yvar=gFeat, ymin=ymin,ymax=ymax,facet_by=NULL)
}


ml<-marrangeGrob(plotList,nrow=3,ncol=2)
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"LFCgeneFeatures_fountains.pdf"),plot=ml, device="pdf",width=29,height=19,units="cm")


##########################-
## gene features in and out of fountains - significant genes
##########################-

stl<-resultsByGRoverlap(fileList,fountains,padjVal=0.05)
stl<-stl[-which(lapply(stl,length)<10)]
stl<-lapply(stl,as.data.frame)
stl<-lapply(stl,dplyr::left_join,geneData,by="wormbaseID")
lapply(stl,dim)
fountRes<-do.call(rbind,stl)
fountRes$regionType<-"fountain"

stl<-resultsByGRoverlap(fileList,nonFount,padjVal=0.05)
stl<-stl[-which(lapply(stl,length)<10)]
stl<-lapply(stl,as.data.frame)
stl<-lapply(stl,dplyr::left_join,geneData,by="wormbaseID")
lapply(stl,dim)
nonFountRes<-do.call(rbind,stl)
nonFountRes$regionType<-"between"

df<-rbind(fountRes,nonFountRes)


gFeatures<-c("numTranscripts","maxTxLength","minTxLength","avrTxLength", "numStarts",
             "numEnds","numExons")
plotList<-list()
for(gFeat in gFeatures){
  res$sampleName<-factor(res$sampleName,levels=unique(res$sampleName))
  ymax=quantile(res[,gFeat],0.99)
  ymin=0
  plotList[[gFeat]]<-pairwiseBoxPlotFunc(df,fountains,yvar=gFeat, ymin=ymin,
                                         ymax=ymax,facet_by="sampleName",
                                         geneSet="significant")
}


ml<-marrangeGrob(plotList,nrow=1,ncol=1)
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"LFCsigGeneFeatures_fountains.pdf"),plot=ml, device="pdf",width=29,height=19,units="cm")
