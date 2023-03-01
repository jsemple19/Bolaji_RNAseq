library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(eulerr)

outPath="."
makeDirs(path=outPath,dirNameList=c("publicData"))
###############-
#  Meeuse et al 2020 - Oscillating genes -----------
###############-
# Developmental function and state transitions of a gene expression oscillator in Caenorhabditis elegan Meeuse...Grosshans Mol Syst Biol (2020)
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7370751/


meeuseURL<-"https://www.embopress.org/action/downloadSupplement?doi=10.15252%2Fmsb.20209498&file=msb209498-sup-0003-DatasetEV1.xlsx"
meeuseFileName<-"MSB-16-e9498-s003.xlsx"

if(!file.exists(paste0(outPath,"/publicData/oscillatingGenes_Meesue.tsv"))){
  download.file(url=meeuseURL,
                destfile=paste0(outPath,"/publicData/",meeuseFileName))

  meeuse<-readxl::read_excel(paste0(outPath,"/publicData/",meeuseFileName),
                             col_types=c(rep("text",3),rep("numeric",2),"text"))
  oscillating<-meeuse[meeuse$Class=="Osc",]
  names(oscillating)[1:3]<-c("wormbaseID","publicID","sequenceID")

  write.table(oscillating,file=paste0(outPath,"/publicData/oscillatingGenes_Meeuse.tsv"),
              row.names=F, col.names=T,quote=F,sep="\t")

  file.remove(paste0(outPath,"/publicData/",meeuseFileName))
}

###############-
#  Latorre et al 2015 - Oscillating genes ------------------
###############-
# The DREAM complex promotes gene body H2A.Z for target repression Latorre..Ahringer (2015)
# https://pubmed.ncbi.nlm.nih.gov/25737279/
latorreURL<-"http://genesdev.cshlp.org/content/suppl/2015/03/03/29.5.495.DC1/Supplemental_TableS7.xlsx"
latorreFileName<-"Supplemental_TableS7.xlsx"

if(! file.exists(paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv"))){
  download.file(url=latorreURL,destfile=paste0(outPath,"/publicData/",latorreFileName))

  latorre<-readxl::read_excel(paste0(outPath,"/publicData/",latorreFileName),
                              col_names=F)
  colnames(latorre)<-"Osc_Latorre2015"

  #metadata<-readRDS(paste0(outPath,"/wbGeneGR_WS275.rds"))
  sum(latorre$Osc_Latorre2015 %in% metadata$sequenceID)
  #3235
  length(latorre$Osc_Latorre2015)
  #3269

  idx<-match(latorre$Osc_Latorre2015,metadata$sequenceID)
  latorre$wormbaseID<-metadata$wormbaseID[idx]
  latorre<-latorre[!is.na(latorre$wormbaseID),]

  write.table(latorre,file=paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv"),
              row.names=F,col.names=T,quote=F,sep="\t")

  file.remove(paste0(outPath,"/publicData/",latorreFileName))
}

latorre<-read.delim(file=paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv"),
                    header=T,sep="\t")
osc<-read.delim(paste0(outPath,"/publicData/oscillatingGenes_Meeuse.tsv"))
sum(latorre$Osc_Latorre2015 %in% osc$SequenceName)
#2473
dim(latorre)
#3269
dim(osc)
#3739

OSC<-list(Meeuse=osc$wormbaseID,Latorre=latorre$wormbaseID)
length(unique(unlist(OSC))) #4522
fit<-euler(OSC)
#pdf(paste0(outPath, "/publicData/venn_MeeuseVsLatorre.pdf"),width=5, height=10,
#    paper="a4")
p1<-plot(fit, quantities=list(type="counts"),
         main=list(label=paste0("Meeuse(2020) vs Latorre(2015) oscillating genes"),
                   fontsize=8, y=0.7))
print(p1)
#dev.off()
