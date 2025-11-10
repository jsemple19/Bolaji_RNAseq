#https://github.com/bartongroup/RATS
#http://127.0.0.1:30134/session/Rvig.121c731ae36d.html
#https://f1000research.com/articles/8-213/v1

library("data.table")
library("matrixStats")
library("ggplot2")
library("rhdf5")
library("rtracklayer")
library("GenomicRanges")
library("ggbio")
library("rats")
library("AnnotationDbi")
library("parallel")
library("dplyr")

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
source(paste0(outPath,"/RATs_functions.R"))


clrs<-getColours()
fileNamePrefix=paste0("RATs/all_")
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

# load a txdb of wormbase data and create a tx2gene object
txdb<-loadDb(paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer,
                    ".annotations.sqlite"))
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
tx2gene$TXNAME<-gsub("Transcript:","",tx2gene$TXNAME)
tx2gene$GENEID<-gsub("Gene:","",tx2gene$GENEID)
colnames(tx2gene)<-c("target_id","parent_id")

##############-
# RATs analysis
##############-

sf_dirs<-file.path("salmon/mRNA",fastqList$sampleName)



## convert file from salmon to kallisto format
#prepare_fish_for_sleuth(sf_dirs)

# prepare sampleTable and design formula
sampleTable<-data.frame(sample=fastqList$sampleName,
                        condition=factor(fastqList$strain,levels=c("366","828")),
                        path=sf_dirs,
                        lane=fastqList$lane,
                        replicate=fastqList$repeatNum,
                        libSize=getLibSizes(sf_dirs)$size/1e6)

design<- formula("~ lane + replicate + condition")
reduced<- formula("~ lane  + replicate")

head(txMeta)
numTxpts<-data.frame(txMeta) %>% dplyr::group_by(wormbaseID) %>%
  dplyr::mutate(numTxpts=n()) %>% select(numTxpts)
if(filterOscillating){
  idx<-tx2gene$target_id %in% txMeta$txptSeqID[txMeta$class=="protein_coding_gene" & txMeta$oscillating=="no"]
} else {
  idx<-tx2gene$target_id %in% txMeta$txptSeqID[txMeta$class=="protein_coding_gene"]
}

rats<-fish4rodents(A_paths=sampleTable$path[sampleTable$condition=="366"],
             B_paths=sampleTable$path[sampleTable$condition=="828"],
             annot=tx2gene[idx,])

# these are TPM (abundance) that are imported
ratsdtu <- call_DTU(annot= tx2gene[idx,], boot_data_A=rats[[1]],
                  boot_data_B=rats[[2]], verbose= T,
                  scaling= 1,#sampleTable$libSize,
                  dprop_thresh = 0.1,
                  abund_thresh = 5,
                  name_A= "TEV_only", name_B= "COH-1cs",
                  varname= "strain",
                  description="Coh-1 cleavage analysis using RATs.",
                  threads=2)


saveRDS(ratsdtu,file=paste0(outPath,"/",fileNamePrefix,"RATsDTU_coh-1.RDS"))
ratsdtu<-readRDS(file=paste0(outPath,"/",fileNamePrefix,"RATsDTU_coh-1.RDS"))
print( names(ratsdtu) )
print( names(ratsdtu$Parameters) )
print( names(ratsdtu$Genes) )
print( names(ratsdtu$Transcripts) )
print( names(ratsdtu$Abundances) )
print( head(ratsdtu$Abundances[[1]]) )


print( dtu_summary(ratsdtu) )
ids <- get_dtu_ids(ratsdtu)
print( names(ids) )
print(lapply(ids,head))

### Summary of isoform switching
print( dtu_switch_summary(ratsdtu) )
ids <- get_switch_ids(ratsdtu)
print( names(ids) )

### Summary of DTU plurality
# how many isoforms are affected per gene.
print( dtu_plurality_summary(ratsdtu) )
# The gene IDs displaying isoform switching.
ids <- get_plurality_ids(ratsdtu)

## Grouping by condition (DEFAULT)
gene_id="WBGene00007810"
gene_id="WBGene00004804" # skn-1
gene_id="WBGene00000591" # coh-1

plot_gene(ratsdtu, gene_id, style="bycondition")

# Grouping by isoform:
plot_gene(ratsdtu, gene_id, style="byisoform")


#models<-annot2models(paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer,
#                     "_ce11.annotations.gtf"),threads=2)
#saveRDS(models,paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer,
#                            "_ce11.annotations.RDS"))

models<-readRDS(paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer,
                      "_ce11.annotations.RDS"))

tx<-unlist(models[[gene_id]])
tx$target_id=names(tx)
lfc<-ratsdtu[["Transcripts"]][ratsdtu[["Transcripts"]]$parent_id==gene_id,]

tx$score<-left_join(data.frame(tx),lfc)$log2FC
txl<-split(tx,tx$target_id)


autoplot(txl,aes(fill=score))

lapply(models[[gene_id]],length)

# Proportion change VS transcript-level significance. Each point is a transcript
plot_overview(ratsdtu, type="tvolcano")

# This can also be plotted for genes, by using the largest isoform effect size as proxy.
plot_overview(ratsdtu , type="gvolcano")


# Distribution of proportion change.
plot_overview(ratsdtu, type="dprop")

# Distribution of largest isoform proportion change per gene.
plot_overview(ratsdtu, type="maxdprop")

# Proportion change VS transcript-level significance. Each point is a transcript
plot_overview(ratsdtu, type="fcvolcano")

# This can also be plotted for genes, by using the largest isoform effect size as proxy.
plot_overview(ratsdtu, type="fcVSdprop")

## Diagnostic plots
# Pairwise Pearson's correlations among samples.
plot_diagnostics(ratsdtu, type='cormat') # Default type.


# Start the interactive volcano plot.
#plot_shiny_volcano(mydtu)

tpm366<-ratsdtu$Abundances$condA[parent_id==gene_id,1:6]
tpm828<-ratsdtu$Abundances$condB[parent_id==gene_id,1:6]

ratio366<-tpm366[1,]/tpm366[2,]
ratio828<-tpm828[1,]/tpm828[2,]

df<-data.frame(condition=rep(c("TEVonly","COH-1cs"),each=6),
           skn1aVskn1b=c(t(ratio366),t(ratio828)))

ggplot(df,aes(x=condition,y=skn1aVskn1b)) + geom_boxplot()+ geom_jitter()



