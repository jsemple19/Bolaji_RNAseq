library(jsonlite)

outPath="."
dioptVer="v9"
sourceOrg<-"6239"
destOrg<-"9606"
filter<-"none"
entrez<-read.delim(file="wormbaseIDtoEntrezID.tsv")
colnames(entrez)<-c("wormbaseID","entrezID","Species","Gene.Name")

orthologs=NULL
for(e in 16333:nrow(entrez)) {
  url=paste0("https://www.flyrnai.org/tools/diopt/web/diopt_api/",dioptVer,
           "/get_orthologs_from_entrez/",sourceOrg,"/",entrez$entrezID[e],"/",
           destOrg,"/",filter)
  js<-fromJSON(url,simplifyVector=F)
  js
  queryID<-js$search_details$search_gene_entrez

  if(length(js$results)>0) {
    numTargets<-length(js$results[[1]])
    targets<-names(js$results[[1]])
    for(i in 1:numTargets) {
        numFeat<-length(js$results[[1]][targets[[i]]][[1]])
        method<-data.frame(method=paste(unlist(js$results[[1]][targets[[i]]][[1]][[1]]),collapse=","))
        feat<-data.frame(js$results[[1]][targets[[i]]][[1]][2:numFeat])
        linedf<-cbind(entrez[e,c("wormbaseID","entrezID")],
                      feat,method)
        if(is.null(orthologs)){
          orthologs<-linedf
        } else {
          orthologs<-rbind(orthologs,linedf)
        }
    }
  # } else { # NA record for genes without orthologs
  #   linedf<-data.frame(wormbaseID= entrez[e,c("wormbaseID")],
  #                      entrezID=entrez[e,c("entrezID")],
  #                      best_score=NA, score=NA,max_score=NA, best_score_rev=NA,
  #                      best_score_count=NA, confidence=NA, mist_ppi=NA,
  #                      mist_genetic=NA,geneid=NA, species_id=NA,symbol=NA,
  #                      species_specific_geneid=NA, species_specific_geneid_type=NA,
  #                      count=NA, method=NA)
  #   if(is.null(orthologs)){
  #     orthologs<-linedf
  #   } else {
  #     orthologs<-rbind(orthologs,linedf)
  #   }
  }
  if(e%%100==0){
    print(e)
    print(dim(orthologs))
  }
}

orthologs1<-orthologs
write.table(orthologs,paste0(outPath,"/publicData/DIOPTorthologs_worm-human_20230615.tsv"),
            sep="\t",quote=F,row.names=F,col.names=T)

#how many worm genes have orthologs
length(unique(orthologs$wormbaseID[!is.na(orthologs$geneid)])) # 11821
length(unique(orthologs$symbol[!is.na(orthologs$geneid)]))# 15749
table(orthologs$confidence)
#    high      low moderate
#7948   181946    36664
table(orthologs$best_score)
#    No    Yes
#190294  36264
table(orthologs$best_score_rev)
#No    Yes
#192607  33951
table(orthologs$best_score_count)
#0      1      2
#167315  48271  10972
orthologs[orthologs$best_score_count==1,] # not clear what this count means, as in some cases tehre are >2 best scores

orthologs[orthologs$wormbaseID=="WBGene00008700",]
