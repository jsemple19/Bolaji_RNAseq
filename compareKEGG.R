library(dplyr)
library(tibble)
library(clusterProfiler)
library(pathview)
library(enrichplot)


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

scriptName="kegg"
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

### upregulated genes
sigUp<-data.frame(Wormbase.ID=getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                                                  namePadjCol="padj",
                                                  nameLfcCol="log2FoldChange",
                                                  direction="gt",
                                                  chr="all", nameChrCol="chr")$wormbaseID)

write.table(sigUp,file="sigUp_DEseq2.csv",row.names=F,quote=F,col.names=F)
sigDown<-data.frame(Wormbase.ID=getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                                                    namePadjCol="padj",
                                                    nameLfcCol="log2FoldChange",
                                                    direction="lt",
                                                    chr="all", nameChrCol="chr")$wormbaseID)
write.table(sigDown,file="sigDown_DEseq2.csv",row.names=F,quote=F,col.names=F)

