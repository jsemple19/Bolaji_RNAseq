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

prepare_fish_for_sleuth(sf_dirs)

sampleTable<-data.frame(sample=fastqList$sampleName,
                        condition=factor(fastqList$strain,levels=c("366","828")),
                        path=sf_dirs,
                        lane=fastqList$lane,
                        replicate=fastqList$repeatNum)

design<- formula("~ lane + replicate + condition")
reduced<- formula("~ lane  + replicate")

tx2gene<-getTx2Gene(genomeDir,genomeVer)
tx2gene$GENEID<-gsub("Gene:","",tx2gene$GENEID)
colnames(tx2gene)<-c("target_id","gene_id")

txMeta<-getTxMetadataGR(genomeDir,genomeVer)
#txMeta<-txMeta[txMeta$class=="protein_coding_gene"]
txMeta<-tagOscillating(txMeta)

bioFilter<-txMeta$txptSeqID[txMeta$class=="protein_coding_gene" &
                             txMeta$oscillating=="no"]


so<-sleuth_prep(sampleTable,normalize=F)
counts<-sleuth_to_matrix(so,"obs_raw","est_counts")
filter<-apply(counts,1,basic_filter)
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

so <- sleuth_fit(so)

models(so)

# Wald test for differential expression of isoforms
oe <- sleuth_wt(so,
                which_beta = 'condition828')


p1<-plot_pc_variance(oe)
p2<-plot_loadings(oe)
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
plot_qq(oe,,test="condition828",test_type="wt")
plot_ma(oe,test="condition828",test_type="wt")
#sleuth_live(oe)

write.csv(res,paste0(outPath,"/sleuth/coh1cs_DTE.csv"))
gr<-makeGRangesFromDataFrame(res,keep.extra.column=T)
saveRDS(gr,paste0(outPath,"/sleuth/coh1cs_DTE.RDS"))

# TODO: DTU analysis
#https://bioconductor.org/packages/release/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html#getting-help
#https://bioconductor.org/packages/release/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html
