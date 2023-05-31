#https://ycl6.gitbook.io/guide-to-rna-seq-analysis/differential-expression-analysis/differential-transcript-usage/dtu-using-drimseq
#http://bioconductor.org/packages/release/bioc/vignettes/DRIMSeq/inst/doc/DRIMSeq.pdf
#https://bioconductor.org/packages/release/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html#ref-Nowicka2016DRIMSeq

library(data.table)
library(GenomicFeatures)
library(tximport)
library(dplyr)
library(reshape2)
library(ggplot2)
library(DRIMSeq)
library(stageR)
library(ggbeeswarm)
library(annotables)
library(tidyverse)
library(AnnotationDbi)
library(ggplot2)
library(ggpubr)
library(BSgenome.Celegans.UCSC.ce11)
library(ggbio)
library(patchwork)
library(gridExtra)
library(DEXSeq)

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
#source(paste0(outPath,"/Sleuth_functions.R"))


clrs<-getColours()
fileNamePrefix=paste0("DEXseq/")
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


################-
## DRIMseq analysis ------
################-
quantFiles<-file.path("salmon/mRNA",fastqList$sampleName,"quant.sf")
names(quantFiles)<-fastqList$sampleName

# load a txdb of wormbase data and create a tx2gene object
txdb<-loadDb(paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer,
                    ".annotations.sqlite"))
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
tx2gene$TXNAME<-gsub("Transcript:","",tx2gene$TXNAME)
tx2gene$GENEID<-gsub("Gene:","",tx2gene$GENEID)

txi = tximport(quantFiles, type = "salmon", tx2gene = tx2gene,
               txIn = TRUE, txOut = TRUE, countsFromAbundance = "dtuScaledTPM")

## cbuild data frame from counts
cts = txi$counts
cts = cts[rowSums(cts) > 0,]
range(colSums(cts)/1e6)


txdf.sub = tx2gene[match(rownames(cts),tx2gene$TXNAME),]
counts = data.frame(gene_id = txdf.sub$GENEID, feature_id = txdf.sub$TXNAME, cts)
dim(counts)
counts[1:5,]


if(filterOscillating){
  idx<-row.names(counts) %in% txMeta$txptSeqID[txMeta$oscillating=="no"]
  counts<-counts[idx,]
}



# prepare sampleTable and design formula
sampleTable<-data.frame(sample_id=paste0("X",fastqList$sampleName),
                        group=factor(fastqList$strain,levels=c("366","828")),
                        path=quantFiles,
                        lane=fastqList$lane,
                        replicate=fastqList$repeatNum)

# create a dmDSdata object d
d = dmDSdata(counts = counts, samples = sampleTable)
table(table(counts(d)$gene_id))
#10247 genes have only one transcript
#    1     2     3     4     5     6     7     8     9    10    11    12    13
#10247  2194   687   294   146    68    41    37    20    10    13     8     4
#14    15    16    17    22
#3     1     1     1     1

# filtering
n = nrow(sampleTable)
n.small = min(table(sampleTable$group))

d = dmFilter(d,
             min_samps_feature_expr = n.small, min_feature_expr = 5,
             min_samps_feature_prop = n.small, min_feature_prop = 0.05,
             min_samps_gene_expr = n, min_gene_expr = 5)

table(table(counts(d)$gene_id))
#   2    3    4    5    6    7    8
#1381  369  111   47   13    4    2
#~2000 genes left after filtering

