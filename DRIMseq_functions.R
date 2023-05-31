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




#' Post-hoc filtering on the standard deviation in proportions
#' https://ycl6.gitbook.io/guide-to-rna-seq-analysis/differential-expression-analysis/differential-transcript-usage/dtu-using-drimseq
smallProportionSD <- function(d, filter = 0.1) {
  # Generate count table
  cts = as.matrix(subset(counts(d), select = -c(gene_id, feature_id)))
  # Summarise count total per gene
  gene.cts = rowsum(cts, counts(d)$gene_id)
  # Use total count per gene as count per transcript
  total.cts = gene.cts[match(counts(d)$gene_id, rownames(gene.cts)),]
  # Calculate proportion of transcript in gene
  props = cts/total.cts
  rownames(props) = rownames(total.cts)

  # Calculate standard deviation
  propSD = sqrt(rowVars(props))
  # Check if standard deviation of per-sample proportions is < 0.1
  propSD < filter
}



# Plot expression function from https://ycl6.gitbook.io/guide-to-rna-seq-analysis/differential-expression-analysis/differential-transcript-usage/dtu-using-drimseq
plotExpression <- function(expData = NULL, geneID = NULL, samps = NULL, isProportion = FALSE){
  colnames(expData)[1:2] = c("gid","tid")
  sub = subset(expData, gid == geneID)
  publicID<-unique(sub$publicID)
  colidx<-colnames(sub) %in% sampleTable$sample_id
  sub<-cbind(sub[,c(1,2)],sub[,colidx])
  sub = reshape2::melt(sub, id = c("gid", "tid"))
  sub = merge(samps, sub, by.x = "sample_id", by.y = "variable")
  if(!isProportion) {
    sub$value = log(sub$value)
  }

  clrs = c("dodgerblue3", "maroon2",  "forestgreen", "darkorange1", "blueviolet", "firebrick2",
           "deepskyblue", "orchid2", "chartreuse3", "gold", "slateblue1", "tomato" , "blue", "magenta", "green3",
           "yellow", "purple3", "red" ,"darkslategray1", "lightpink1", "lightgreen", "khaki1", "plum3", "salmon")

  p = ggplot(sub, aes(tid, value, color = group, fill = group)) +
    geom_boxplot(alpha = 0.4, outlier.shape = NA, width = 0.8, lwd = 0.5) +
    stat_summary(fun = mean, geom = "point", color = "black", shape = 5, size = 3, position=position_dodge(width = 0.8)) +
    scale_color_manual(values = clrs) + scale_fill_manual(values = clrs) +
    geom_quasirandom(color = "black", size = 1, dodge.width = 0.8) + theme_bw() +
    ggtitle(paste0(geneID,": ",publicID)) + xlab("Transcripts")

  if(!isProportion) {
    p = p + ylab("log(Expression)")
  } else {
    p = p + ylab("Proportions")
  }
  p
}
