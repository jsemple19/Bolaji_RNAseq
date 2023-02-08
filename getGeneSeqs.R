library(rtracklayer)
library(Biostrings)
library(BSgenome)

args<-commandArgs(trailingOnly=TRUE)
genomeVer=args[1]
genomeDir=args[2]


gnm<-readDNAStringSet(paste0(genomeDir,"/c_elegans.PRJNA13758.",genomeVer,".genomic.fa"))
gtf<-import(paste0(genomeDir,"/c_elegans.PRJNA13758.",genomeVer,".annotations.gtf"))

gns<-unlist(reduce(split(gtf,gtf$gene_id)))
gtfa<-getSeq(gns,gnm)

writeXStringSet(gtfa,paste0(genomeDir,"/c_elegans.PRJNA13758.",genomeVer,".genes.fa.gz"),compress=T)
