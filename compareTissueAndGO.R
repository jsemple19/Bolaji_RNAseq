library(GenomicRanges)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(plyranges)
library(rtracklayer)
library(dplyr)
library(readxl)
library(wormcat)
library(xlsx)
library(reticulate)

genomeVer="WS285"
outPath="."
scriptPath="."
source(paste0(scriptPath,"/functions.R"))

theme_set(
  theme_bw()+
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
)

scriptName="compareTissueAndGO"
fileNamePrefix=paste0(scriptName,"/")
makeDirs(outPath,dirNameList=c(paste0(c("wormcat/"),scriptName)))

padjVal=0.05
lfcVal=0
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


############################-
## Wormcat-----
############################-
## need to process by dropping excel file into http://www.wormcat.com
## since r function does not seem to work.

if(!dir.exists(paste0(outPath,"/wormcat"))) {
  dir.create(paste0(outPath,"/wormcat"),
             recursive=T)
}

wormcatDataURL="http://www.wormcat.com/static/download/whole_genome_nov-16-2019.csv"
wormcatDataURL="https://raw.githubusercontent.com/dphiggs01/Wormcat/master/inst/extdata/whole_genome_v2_nov-11-2021.csv"
wormcatData=paste0(getwd(),"/publicData/",basename(wormcatDataURL))
if(!file.exists(wormcatData)){
  download.file(url=wormcatDataURL,destfile=wormcatData)
}

grp=fileList$sampleName[1]

for(grp in fileList$sampleName){
  salmon<-readRDS(fileList[fileList$sampleName==grp,"filePath"])

  ### upregulated genes
  sigTable<-data.frame(Wormbase.ID=getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                                                       namePadjCol="padj",
                                                       nameLfcCol="log2FoldChange",
                                                       direction="gt",
                                                       chr="all", nameChrCol="chr")$wormbaseID)
  write.table(sigTable,paste0(outPath,"/wormcat/",grp,
                              "_up_wormcat.csv"),
              row.names=F,quote=F,col.names=T)
  worm_cat_fun(file_to_process=paste0(outPath,"/wormcat/",grp,
                                      "_up_wormcat.csv"),
               title=paste(grp,"up"),
               output_dir=paste0(outPath,"/wormcat/",grp,"_up"),
               annotation_file=wormcatData,
               input_type="Wormbase.ID", rm_dir=FALSE,
               zip_files=FALSE)

  ### down regulated genes
  sigTable<-data.frame(Wormbase.ID=getSignificantGenes(salmon, padj=padjVal, lfc= -lfcVal,
                                                       namePadjCol="padj",
                                                       nameLfcCol="log2FoldChange",
                                                       direction="lt",
                                                       chr="all", nameChrCol="chr")$wormbaseID)
  write.table(sigTable,paste0(outPath,"/wormcat/",grp,"_down_wormcat.csv"),
              row.names=F,quote=F,col.names=T)
  worm_cat_fun( file_to_process=paste0(outPath,"/wormcat/",  grp,
                                       "_down_wormcat.csv"),
                title=paste(grp,"down"),
                output_dir=paste0(outPath,"/wormcat/",grp,"_down"),
                annotation_file=wormcatData,
                input_type="Wormbase.ID", rm_dir=FALSE,
                zip_files=FALSE)
}




########-
## TEA - tissue enrichment analysis-----
########-
# https://www.micropublication.org/media/2018/03/microPublication.biology-10.17912-W25Q2N.pdf

if(!dir.exists(paste0(outPath,"/tissue/tea/"))) {
  dir.create(paste0(outPath,"/tissue/tea/"),
             recursive=T)
}


# fetch dictionaries
# anatomy
anaURL=paste0("http://caltech.wormbase.org/TissueEnrichmentAnalysis/DICTs/anatomy_dict_95_33_", genomeVer, ".csv")
anaDict=paste0(outPath, "/publicData/",basename(anaURL))
if(!file.exists(anaDict)){
  download.file(url=anaURL, destfile=anaDict)
}

# go
goURL=paste0("http://caltech.wormbase.org/TissueEnrichmentAnalysis/DICTs/go_dict_95_33_", genomeVer, ".csv")
goDict=paste0(outPath, "/publicData/",basename(goURL))
if(!file.exists(goDict)){
  download.file(url=goURL, destfile=goDict)
}

# phenotype
pheURL=paste0("http://caltech.wormbase.org/TissueEnrichmentAnalysis/DICTs/phenotype_dict_95_50_", genomeVer, ".csv")
pheDict=paste0(outPath,"/publicData/",basename(pheURL))
if(!file.exists(pheDict)){
  download.file(url=pheURL, destfile=pheDict)
}


condaActivate<-gsub("conda$","activate",conda_binary(conda = "auto"))
sink(file=paste0(outPath,"/runTea.sh"),append=FALSE, type="output")
cat("#! /bin/bash\n")
cat(paste0("source ",condaActivate, " tea\n"))
cat(paste0("cd ./tissue/tea/","\n"))
#cat("cd ./tissue/tea\n")
sink()


## upregulated genes
sigTables<-list()
for(grp in fileList$sampleName){
  salmon<-readRDS(fileList[fileList$sampleName==grp,"filePath"])

  sigTables[[paste0(grp)]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="gt",
                        chr="all", nameChrCol="chr"))
}

sigGenes<-lapply(sigTables, "[", ,"wormbaseID")
lapply(sigGenes,length)
sigGenes<-lapply(sigGenes,na.omit)

for(grp in fileList$sampleName){
  write.table(sigGenes[[grp]], file=paste0(outPath,"/tissue/tea/",grp,"_upGenes_WBID.txt"), quote=F, row.names=F,col.names=F)
  sink(file=paste0(outPath,"/runTea.sh"),append=TRUE, type="output")
  cat(paste0("tea -d ../../",anaDict," -q 0.05 -s ",grp,"_upGenes_WBID.txt ",grp,"_up_tissue tissue\n"))
  cat(paste0("tea -d ../../",pheDict," -q 0.05 -s ",grp,"_upGenes_WBID.txt ",grp,"_up_phe phenotype\n"))
  cat(paste0("tea -d ../../",goDict," -q 0.05 -s ",grp,"_upGenes_WBID.txt ",grp,"_up_go go\n"))
  sink()
}


## downregulated genes
sigTables<-list()
for(grp in fileList$sampleName){
  salmon<-readRDS(fileList[fileList$sampleName==grp,"filePath"])

  sigTables[[paste0(grp)]]<-as.data.frame(
    getSignificantGenes(salmon, padj=padjVal, lfc=lfcVal,
                        namePadjCol="padj",
                        nameLfcCol="log2FoldChange",
                        direction="lt",
                        chr="all", nameChrCol="chr"))
}

sigGenes<-lapply(sigTables, "[", ,"wormbaseID")
lapply(sigGenes,length)
sigGenes<-lapply(sigGenes,na.omit)

for(grp in fileList$sampleName){
  write.table(sigGenes[[grp]], file=paste0(outPath,"/tissue/tea/",grp,"_downGenes_WBID.txt"), quote=F, row.names=F,col.names=F)
  sink(file=paste0(outPath,"/runTea.sh"),append=TRUE, type="output")
  cat(paste0("tea -d ../../",anaDict," -q 0.05 -s ",grp,"_downGenes_WBID.txt ",grp,"_down_tissue tissue\n"))
  cat(paste0("tea -d ../../",pheDict," -q 0.05 -s ",grp,"_downGenes_WBID.txt ",grp,"_down_phe phenotype\n"))
  cat(paste0("tea -d ../../",goDict," -q 0.05 -s ",grp,"_downGenes_WBID.txt ",grp,"_down_go go\n"))
  sink()
}

sink(file=paste0(outPath,"/runTea.sh"),append=TRUE, type="output")
cat("cd ../../\n")
sink()

system(paste0("grep -v '^>' ",outPath,"/runTea.sh > ",outPath,"/runTea1.sh"))
file.remove(paste0(outPath,"/runTea.sh"))
system(paste0("chmod +x ",outPath,"/runTea1.sh"))
system(paste0(outPath,"/runTea1.sh"),wait=F)


## combine images
rmarkdown::render(input="AllPlots_WORMCAT.Rmd", output_format="pdf_document",
                  output_file=paste0(outPath,"/wormcat/","allPlots_WORMCAT.pdf"))

rmarkdown::render(input="AllPlots_TEA.Rmd",output_format="pdf_document",
                  output_file=paste0(outPath,"/tissue/tea/","allPlots_TEA.pdf"))



###### compile results in excel-----
file.remove("allTEAresults.xlsx")
categories=c("go","tissue","phe")
for(grp in fileList$sampleName){
  go<-NULL
  for(category in categories){
    if(file.exists(paste0(outPath,"/tissue/tea/",grp,"_up_",category,".tsv"))){
      up<-read.delim(paste0(outPath,"/tissue/tea/",grp,"_up_",
                          category,".tsv"))
      if(nrow(up)>0){
        up$SMC<-grp
        up$category<-category
        up$UpOrDownRegulated<-"up"
        if(is.null(go)){
          go<-up
        } else {
          go<-rbind(go,up)
        }
      }
    }
    if(file.exists(paste0(outPath,"/tissue/tea/",grp,"_down_",category,".tsv"))){
      down<-read.delim(paste0(outPath,"/tissue/tea/",grp,"_down_",
                            category,".tsv"))
      if(nrow(down)>0){
        down$SMC<-grp
        down$category<-category
        down$UpOrDownRegulated<-"down"
        if(is.null(go)){
          go<-down
        } else {
          go<-rbind(go,down)
        }
      }
    }
  }
  write.xlsx(go,file="allTEAresults.xlsx",sheetName=grp,
             append=T,row.names=F)
}



