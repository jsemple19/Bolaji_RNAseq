library(rtracklayer)
#library(BSgenome.Celegans.UCSC.ce11)

args = commandArgs(trailingOnly=TRUE)
wbVer<-args[1]
genomeDir<-args[2]
url="https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.chrom.sizes"

setwd(genomeDir)
download.file(url,basename(url))

chrSizes<-read.delim(basename(url),header=F)

#chrSizes<-chrSizes[order(chrSizes$V1)[c(1:4,6,7,5)],]

chrSizes$V1<-gsub("chr","",chrSizes$V1)
chrSizes$V1<-gsub("M","MtDNA",chrSizes$V1)

write.table(chrSizes,paste0(wbVer,".chrom.sizes"),row.names=F,col.names=F,quote=F,sep="\t")

