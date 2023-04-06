

#' Create GRanges of gene names
#'
#' Create a GRanges object for protein coding genes from wormbase annotation
#' txdb sqlite file. Requires txdb and geneIDs file from wormbase, both
#' saved into genomeDir. Metadata object uses ucsc format seqnames.
#' @param genomeDir Directory where txdb sqlite and geneIDs files are stored.
#' @param genomeVer Version of wormbase to use, e.g. "WS285"
#' @return Metadata object
#' @export
getTxMetadataGR<-function(genomeDir,genomeVer,outPath=".",remake=F){
  if(!file.exists(paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer, ".annotations.sqlite"))){
    makeTxDbsqlite(genomeDir, genomeVer)
  }
  # load a txdb of wormbase data and create a tx2gene object
  if(!file.exists(paste0(outPath,"/ce11TxptGR_",genomeVer,".rds"))){
    txdb<-loadDb(paste0(genomeDir, "/c_elegans.PRJNA13758.", genomeVer,
                        ".annotations.sqlite"))
    Txpts<-transcripts(txdb,columns=c("tx_id","tx_name","gene_id"))
    print(paste0("Txpts object length: ", length(Txpts)))

    #geneGR<-unlist(range(TxptByGene))
    #mcols(geneGR)$wormbase<-names(geneGR)

    txptdf<-as.data.frame(Txpts)
    colnames(txptdf)<-c("seqnames","start","end","width","strand","txptID","txptSeqID","wormbaseID")

    txptdf$wormbaseID<-gsub("^Gene:","",txptdf$wormbaseID)
    txptdf$txptSeqID<-gsub("^Transcript:","",txptdf$txptSeqID)

    geneIDs<-read.delim(gzfile(paste0(genomeDir,"/c_elegans.PRJNA13758.",genomeVer,".geneIDs.txt.gz")), header=FALSE, sep=",")
    colnames(geneIDs) <- c("species.ID","WormBase.Gene.ID","Public.Name","Sequence.Name","gene.anotation.status","class")
    # can be ignored: david
    #david<-read.delim("/Users/semple/Documents/MeisterLab/GenomeVer/annotations/david_wbid2entrez_WS278.txt")

    ### combine wormbase gene IDs with gene location ###
    metadata<-left_join(txptdf,geneIDs, genedf,by=c("wormbaseID"="WormBase.Gene.ID"))%>%
      dplyr::select(seqnames, start, end, width, strand, txptID, txptSeqID, wormbaseID, Public.Name, Sequence.Name, class) %>%
      collect %>% GenomicRanges::GRanges()

    names(mcols(metadata))<-c( "txptID", "txptSeqID", "wormbaseID", "publicID","sequenceID","class")

    table(metadata$class)

    seqlevelsStyle(metadata)<-"ucsc"
    seqinfo(metadata)<-seqinfo(BSgenome.Celegans.UCSC.ce11::Celegans)
    metadata<-sort(metadata)

    # prtncd<-metadata[metadata$class=="protein_coding_gene"]
    # prtncd$name<-prtncd$wormbaseID
    # export(prtncd,"proteinCodingGenes.bed")
    #
    # ncd<-metadata[metadata$class!="protein_coding_gene"]
    # ncd$name<-ncd$class
    # export(ncd,"nonCodingGenes.bed")

    saveRDS(metadata,paste0(outPath,"/ce11TxptGR_",genomeVer,".rds"))
  } else {
    metadata<-readRDS(paste0(outPath,"/ce11TxptGR_",genomeVer,".rds"))
  }
  return(metadata)
}


#' Plot QC of pvalues and filtering threshold
#'
#' Plot histogram of pvalues and ajusted pvalues. Expect to see uniform distribution
#' and elevated counts in smallest bin.
#' Also plot criteria used to choose independent filtering threshold.
#' @param resLFC DESeqResults object with shrunken LFC
#' @param contrastName String indicating the name of the contrast used to make res
#' @param outPath Path to working directory where results will be written
#' @param fileNamePrefix Prefix to add to the name of the output file.
#' @return Plots written to pdf file
#' @export
plotTxPvalQC<-function(resLFC,contrastName,outPath=".", fileNamePrefix="sleuth_"){
  # histograms of pvals
  pdf(file=paste0(outPath,"/",fileNamePrefix,contrastName,
                  "_QCdge_pvalHist_filtThresh.pdf"), width=5,height=11,paper="a4")
  par(mfrow=c(3,1))
  if("oscillating" %in% names(resLFC)){
    oscFac<-factor(resLFC$oscillating,levels=c("no","Meeuse","Latorre","Meeuse&Latorre"))
    pvalLevels<-cut(resLFC$pval,seq(0,1,0.1),include.lowest=T)
    barplot(table(oscFac, pvalLevels), legend.text=levels(oscFac),
            main=paste0(contrastName," pvalue"))
    pvalLevels<-cut(resLFC$qval,seq(0,1,0.1),include.lowest=T)
    barplot(table(oscFac, pvalLevels), legend.text=levels(oscFac),
            main=paste0(contrastName," adjusted pvalue"))
  } else {
    hist(resLFC$pval,breaks=9,main=paste0(contrastName," pvalue"))
    hist(resLFC$qval,breaks=9, main=paste0(contrastName," adjusted pvalue"))
  }
  par(mfrow=c(1,1))
  dev.off()
}





#' Remove oscillating genes
#'
#' Remove genes that oscillate during larval development. Based on two studies
#' from Meeuse et al. (2020) and Laorre (2015).
#' @param dds DESeq2 object
#' @param remove One of c("either","both","Meeuse","Latorre"). If "either"
#' removes genes present in either one of these studies. "both" removes genes
#' present in both studies. And "Meeuse" and "Latorre" only remove genes from
#' that one particular study
#' @return DESeq2 object with rows corresponding to specificed genes removed
#' @export
removeOscillating<-function(dds,remove="either"){
  if(!file.exists(paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv")) |
     !file.exists(paste0(outPath,"/publicData/oscillatingGenes_Meeuse.tsv"))){
    source("./processPublished.R")
  }
  latorre<-read.delim(file=paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv"),
                      header=T,sep="\t")
  meeuse<-read.delim(paste0(outPath,"/publicData/oscillatingGenes_Meeuse.tsv"))
  if(remove=="either"){
    idlist=unique(c(latorre$wormbaseID,meeuse$wormbaseID))
  }
  if(remove=="both"){
    idx<-which(meeuse$wormbaseID %in% latorre$wormbaseID)
    idlist<-meeuse$wormbaseID[idx]
  }
  if(remove=="Meeuse"){
    idlist<-meeuse$wormbaseID
  }
  if(remove=="Latorre"){
    idlist<-latorre$wormbaseID
  }
  idx<-rowData(dds)$gene %in% idlist
  dds<-dds[!idx,]
  print(paste0("Removing ",sum(idx)," oscillating genes from ",remove," dataset"))
  return(dds)
}

#' Remove oscillating genes
#'
#' Remove genes that oscillate during larval development. Based on two studies
#' from Meeuse et al. (2020) and Laorre (2015).
#' @param metadata GRanges object
#' @return GRanges object with genes tagged by oscillating info
#' @export
tagOscillating<-function(metadata){
  if(!file.exists(paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv")) |
     !file.exists(paste0(outPath,"/publicData/oscillatingGenes_Meeuse.tsv"))){
    source("./processPublished.R")
  }
  latorre<-read.delim(file=paste0(outPath,"/publicData/oscillatingGenes_latorre.tsv"),
                      header=T,sep="\t")
  meeuse<-read.delim(paste0(outPath,"/publicData/oscillatingGenes_Meeuse.tsv"))
  idlist=unique(c(latorre$wormbaseID,meeuse$wormbaseID))
  df<-data.frame(wormbaseID=unique(c(latorre$wormbaseID,meeuse$wormbaseID)),
                 oscillating=NA)
  idx<-which(df$wormbaseID %in% latorre$wormbaseID)
  df$oscillating[idx]<-"Latorre"
  idx<-which(df$wormbaseID %in% meeuse$wormbaseID)
  df$oscillating[idx]<-"Meeuse"
  idx<-which(df$wormbaseID %in% meeuse$wormbaseID & df$wormbaseID %in% latorre$wormbaseID)
  df$oscillating[idx]<-"Meeuse&Latorre"
  idx1<-which(metadata$wormbaseID %in% df$wormbaseID)
  idx2<-match(metadata$wormbaseID[idx1], df$wormbaseID)
  metadata$oscillating<-"no"
  metadata$oscillating[idx1]<-df$oscillating[idx2]
  return(metadata)
}
