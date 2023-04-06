#https://hbctraining.github.io/DGE_workshop_salmon/lessons/09_sleuth.html
#https://liorpachter.wordpress.com/2015/08/17/a-sleuth-for-rna-seq/
#https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html

library(wasabi)
library(sleuth)
library(annotables)
library(tidyverse)
library(AnnotationDbi)
library(ggplot2)
library(ggpubr)
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
fastqList_file <- "./fastqList_coh1.csv"
fastqList<-read.csv(fastqList_file)

grp="coh1"

source(paste0(outPath,"/functions.R"))
source(paste0(outPath,"/DESeq2_functions.R"))
source(paste0(outPath,"/Sleuth_functions.R"))

makeDirs(outPath,c("sleuth"))
clrs<-getColours()
fileNamePrefix=paste0("sleuth/")

print(paste("fastqList_file is: ",fastqList_file))
print(paste("outPath is: ",outPath))
print(paste("genomeDir is: ",genomeDir))

sf_dirs<-file.path("salmon/mRNA",fastqList$sampleName)


## convert file from salmon to kallisto format
prepare_fish_for_sleuth(sf_dirs)

# prepare sampleTable and design formula
sampleTable<-data.frame(sample=fastqList$sampleName,
                        condition=factor(fastqList$strain,levels=c("366","828")),
                        path=sf_dirs,
                        lane=fastqList$lane,
                        replicate=fastqList$repeatNum)

design<- formula("~ lane + replicate + condition")
reduced<- formula("~ lane  + replicate")

# prepare metadata
tx2gene<-getTx2Gene(genomeDir,genomeVer)
tx2gene$GENEID<-gsub("Gene:","",tx2gene$GENEID)
colnames(tx2gene)<-c("target_id","gene_id")

txMeta<-getTxMetadataGR(genomeDir,genomeVer)
#txMeta<-txMeta[txMeta$class=="protein_coding_gene"]
txMeta<-tagOscillating(txMeta)


############-
# Sleuth
############-

# read in data to sleuth
so<-sleuth_prep(sampleTable,normalize=F)

# filter input genes by counts and using metadata
counts<-sleuth_to_matrix(so,"obs_raw","est_counts")
filter<-apply(counts,1,basic_filter)

bioFilter<-txMeta$txptSeqID[txMeta$class=="protein_coding_gene" &
                              txMeta$oscillating=="no"]
txFilter<-bioFilter[which(bioFilter %in% names(filter)[filter])]

so<-sleuth_prep(sampleTable,
                full_model=design,
                target_mapping=tx2gene,
                num_cores=1,
                read_bootstrap_tpm=T,
                extra_bootstrap_summary=T,
                transformation_function=function(x) log2(x + 0.5),
                filter_target_id=txFilter)

# 18363 targets passed the filter without gene list
# 13198 targets passed the filter when removing oscillating

# fit model
so <- sleuth_fit(so)

# get names of coefficients
models(so)

# Wald test for differential expression of isoforms
oe <- sleuth_wt(so,
                which_beta = 'condition828')


p1<-plot_pc_variance(oe)
p2<-plot_loadings(oe,pc_input=1L)
p2a<-plot_loadings(oe,pc_input=2L)

p3<-plot_pca(oe,
         color_by = 'replicate',
         text_labels = TRUE)

p4<-plot_pca(oe,
         color_by = 'condition',
         text_labels = TRUE)

p4a<-plot_pca(oe, pc_x=3L, pc_y=4L,
             color_by = 'condition',
             text_labels = TRUE)

#plot_sample_heatmap(oe)

p5<-plot_group_density(oe,
                   use_filtered = FALSE,
                   units = "est_counts",
                   trans = "log",
                   grouping = "condition")

p6<-plot_group_density(oe,
                   use_filtered = TRUE,
                   units = "est_counts",
                   trans = "log",
                   grouping = "condition")


p<-ggpubr::ggarrange(p1,p2,p2a,p3,p4,p4a,p5,p6,nrow=3,ncol=3)
ggsave(paste0(outPath,"/sleuth/QC_sleuth.pdf"),p,dev="pdf",
       height=13,width=15)

# output results
res <- sleuth_results(oe,
                      test = 'condition828',
                      show_all = TRUE)

res<-left_join(res,data.frame(txMeta),by=join_by(target_id==txptSeqID))

plotTxPvalQC(res,contrastName=grp,outPath=outPath, fileNamePrefix=fileNamePrefix)


# look at skn-1
skn1<-res$target_id[res$publicID=="skn-1"]
tmp<-res[res$publicID=="skn-1",]
tmp<-tmp[complete.cases(tmp),]
tmp$txt<-round(tmp$qval,2)
tmp$txt[tmp$txt>0.05]<-"n.s."
tmp$txtPos<-0.2
tmp$label<-c("a","b")

p<-ggplot(tmp,aes(x=label,y=b)) +
  geom_bar(stat="identity",fill="lightgrey",color="black") +
  xlab("*skn-1* isoform") + ylab(label="Log<sub>2</sub>FC") +
  geom_errorbar(aes(ymin=b-se_b,ymax=b+se_b),width=0.2) +
  geom_hline(yintercept=0) + coord_cartesian(ylim=c(-0.23,0.21)) +
  geom_text(aes(y=txtPos,label=txt))
p
ggsave(paste0(outPath,"/sleuth/skn1_bar.pdf"),p,dev="pdf",
       height=3,width=1.5)

# counts1<-sleuth_to_matrix(so,"obs_norm","est_counts")
# tmp<-counts1[grep("T19E7.2",rownames(counts1))[1:2],]
# tmp
# ratios<-data.frame(tmp[1,]/tmp[2,])
# ratios$name<-gsub("_.*","",rownames(ratios))
# colnames(ratios)<-c("ratio","sample")
# ggplot(ratios,aes(x=sample,y=ratio)) + geom_jitter(width=0.1)

# plotList<-list()
# for(tid in skn1){
#   if(tid %in% res$target_id & !is.na(res[res$target_id==tid,"pval"])){
#     plotList[[tid]]<-
#       plot_bootstrap(oe,
#                target_id = tid,
#                units = "tpm",
#                color_by = "condition")
#   }
# }
#
# p<-ggpubr::ggarrange(plotlist=plotList,nrow=2,ncol=2)
# ggsave(paste0(outPath,"/sleuth/skn1.pdf"),p,dev="pdf",
#        height=13,width=15)



sigtxpts <- res %>%
  filter(qval < 0.05)

#plot_transcript_heatmap(oe,
#                        transcripts = sigtxpts$target_id[1:20])

plot_bootstrap(oe,
               target_id = sigtxpts$target_id[1],
               units = "est_counts",
               color_by = "condition")
#so$sample_to_covariates

plot_mean_var(oe)
plot_qq(oe,test="condition828",test_type="wt")
plot_ma(oe,test="condition828",test_type="wt")
#sleuth_live(oe)

write.csv(res,paste0(outPath,"/sleuth/coh1cs_DTE.csv"))
gr<-makeGRangesFromDataFrame(res,keep.extra.column=T)
seqinfo(gr)<-seqinfo(Celegans)
gr<-sort(gr)
saveRDS(gr,paste0(outPath,"/sleuth/coh1cs_DTE.RDS"))

#######################-
# find DTE ------
#######################-
so<-readRDS(paste0(outPath,"/sleuth/coh1cs_DTE.RDS"))

upGenes<-data.frame(so) %>% dplyr::filter(!is.na(qval),qval<0.05,b>0) %>% dplyr::select(wormbaseID) %>% unique()

uptx<-data.frame(so) %>% dplyr::filter(wormbaseID %in% upGenes$wormbaseID) %>%
  dplyr::group_by(wormbaseID) %>%
  dplyr::mutate(count=n(),notNA=which(!is.na(qval))[1]) %>%
  dplyr::filter(count>1) %>% mutate(isoRatio=b/b[notNA]) %>%
  dplyr::filter(!is.na(isoRatio),abs(isoRatio)>1)

print(uptx,width=Inf)

diffuptx<-unique(uptx$wormbaseID)
so$diffuptx<-F
so$diffuptx[so$wormbaseID %in% diffuptx]<-T


fountains<-readRDS("/Users/semple/Documents/MeisterLab/otherPeopleProjects/fountains/detected_fountains_equalQ.RDS")
fountains$fountainName<-paste0("fount",1:length(fountains))
fountains<-resize(fountains, width=6000,fix="center")

nonFount<-gaps(fountains)
nonFount<-resize(nonFount,width=6000,fix="center")

so$olFount<-F
ol<-findOverlaps(so[so$diffuptx],fountains,ignore.strand=T)
so$olFount[queryHits(ol)]<-T
table(so$olFount,so$diffuptx)

so$olGap<-F
ol<-findOverlaps(so[so$diffuptx],nonFount,ignore.strand=T)
so$olGap[queryHits(ol)]<-T
table(so$olGap,so$diffuptx)

library(ggbio)
library(GENOVA)
library(patchwork)
# get hic data
tev<-load_contacts(signal_path="/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/2021_HiCworked/366_cis100bins_1000.cool", sample_name= "TEVonly", centromeres=F, balancing = F, verbose = T, colour="black")
coh1<-load_contacts(signal_path="/Users/semple/Documents/MeisterLab/otherPeopleProjects/Moushumi/2021_HiCworked/828_cis100bins_1000.cool", sample_name= "coh1", centromeres=F, balancing = F, verbose = T, colour="black")
# get fountains
fountains<-readRDS("/Users/semple/Documents/MeisterLab/otherPeopleProjects/fountains/detected_fountains_equalQ.RDS")
fountains$fountainName<-paste0("fount",1:length(fountains))
# get daugherty enhancers
daugherty<-readRDS("/Users/semple/Documents/MeisterLab/otherPeopleProjects/fountains/publicData/daugherty2017_L3enhancers_ce11.rds")
activedaugherty<-daugherty[daugherty$L3_chromHMMState=="L3_activeEnhancer"]
repdaugherty<-daugherty[daugherty$L3_chromHMMState=="L3_repressedEnhancer" | daugherty$L3_chromHMMState=="L3_H3K27me3Repressed"]
# get jaenes enhancers
jaenes<-import("/Users/semple/Documents/MeisterLab/otherPeopleProjects/fountains/publicData/Jaenes2018_enhancers_ce11_stages_L3chromHMM.bed")
activejaenes<-jaenes[jaenes$name=="Active enhancer"]
repjaenes<-jaenes[jaenes$name=="Repressed enhancer" | jaenes$name=="H3K27me3 repressed"]

# get transcripts
if(!file.exists(paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer,
                       "_ce11.annotations.sqlite"))){
  makeTxDbsqlite_ce11(genomeDir,genomeVer)
}
txdb<-loadDb(paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer,
                    "_ce11.annotations.sqlite"))

dir.create(paste0(outPath,"/sleuth/sleuthDTE"))

theme_bed = function(){
  theme(title=element_text(size=8),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank())
}


plotList=list()
for (wbid in unique(so$wormbaseID[so$diffuptx==T])) {
#wbid="WBGene00020131"
  gr<-GenomicRanges::reduce(so[so$wormbaseID==wbid],drop.empty.ranges=T)
  #f<-nearest(gr,fountains,ignore.strand=T)
  #gr<-reduce(fountains[f],gr)
  gr<-trim(resize(gr,width=100000,fix="center"))
  p.bg <- autoplot(Celegans, which = gr) + xlim(gr)
  p.txdb<-autoplot(txdb,which=gr,label.size=1.5,fill="grey80") + xlim(gr)
  p.fount <- autoplot(subsetByOverlaps(fountains,gr), geom="rect",aes(fill=j.score),legend=F) +
    xlim(gr) + ggtitle("Detected fountains") + theme_bed()
  p.activeD<-autoplot(subsetByOverlaps(activedaugherty,gr), geom="rect",fill="darkgreen") +
    xlim(gr) + ggtitle("Active enhancers (Daugherty et al.)") + theme_bed()
  p.activeJ<-autoplot(subsetByOverlaps(activejaenes,gr), geom="rect",fill="darkgreen") +
    xlim(gr) + ggtitle("Active enhancers (Jaenes et al.)") #+ theme_bed()
  p.repD<-autoplot(subsetByOverlaps(repdaugherty,gr), aes(fill=L3_chromHMMState),geom="rect",fill="darkgray") +
    xlim(gr) + ggtitle("Repressed enhancers (Daugherty et al.)") + theme_bed()
  p.repJ<-autoplot(subsetByOverlaps(repjaenes,gr), aes(fill=name),geom="rect",fill="darkgray") +
    xlim(gr)+ ggtitle("Repressed enhancers (Jaenes et al.)") + theme_bed()


  p.hic<-pyramid(exp=tev,chrom=as.character(seqnames(gr)),start=start(gr),end=end(gr),crop_y=100000) +
    ggtitle("TEVonly") + theme(legend.key.size = unit(0.2,'cm'))
  p.hic1<-pyramid(exp=coh1,chrom=as.character(seqnames(gr)),start=start(gr),end=end(gr),crop_y=100000) +
    ggtitle("COH-1cs") + theme(legend.key.size = unit(0.2,'cm'))
  p.hicdiff<-pyramid_difference(coh1,tev,chrom=as.character(seqnames(gr)),start=start(gr),end=end(gr),crop_y=100000) +
    ggtitle("COH-1cs - TEV") + theme(legend.position=c(0.9,0.8))
# look at transcripts
  pptxdb<-autoplot(subsetByOverlaps(so,gr),aes(fill=b)) +xlim(gr)+
    scale_fill_gradient2(na.value="grey90") + theme(legend.key.size = unit(0.2,'cm'))+
    labs(fill="LFC")

  p1<-(p.hic+p.hic1)/p.hicdiff/p.txdb@ggplot/pptxdb@ggplot/p.fount@ggplot/plot_spacer()/p.activeD@ggplot/plot_spacer()/p.activeJ@ggplot/plot_spacer()/p.repD@ggplot/plot_spacer()/p.repJ@ggplot/plot_spacer()/p.bg@ggplot +
    plot_layout(heights=c(10,25,10,5,0.5,-2,0.5,-2,0.5,-2,0.5,-2,0.5,-2,0.1)) +
    plot_annotation(title=paste0(wbid,": ",so$publicID[so$wormbaseID==wbid]))
  ggsave(paste0(outPath,"/sleuth/sleuthDTE/hic_",wbid,".pdf"),p1,device="pdf",height=29,width=19,units="cm")


  tmp<-data.frame(subsetByOverlaps(so[so$wormbaseID==wbid],gr))
  tmp$txt<-paste0("p=",round(tmp$qval,2))
  tmp$txt[tmp$qval>0.05]<-"n.s."

  plotList[[wbid]]<-ggplot(tmp,aes(x=target_id,y=b,fill=b)) +
    geom_bar(stat="identity",color="black") +
    xlab(paste0(tmp$wormbaseID, " isoforms")) + ylab(label="Log<sub>2</sub>FC") +
    geom_errorbar(aes(ymin=b-se_b,ymax=b+se_b),width=0.2) +
    geom_hline(yintercept=0) + scale_fill_gradient2() +
    geom_text(aes(y=0,label=txt))
}





# TODO: DTU analysis
#https://bioconductor.org/packages/release/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html#getting-help
#https://bioconductor.org/packages/release/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html

