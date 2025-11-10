library(rtracklayer)
library(zoo)
library(BSgenome.Celegans.UCSC.ce11)

bw<-import("/Users/semple/Documents/MeisterLab/Datasets/ChIP_seq_data_GEO/modEncode_SMC/SDQ0809_COH1_fem2_AD_ChIP_Rep1.bigwig")

seqlevels(bw)
seqlevels(bw)<-seqlevels(Celegans)

cov<-coverage(bw,weight = "score")

tileWidth=100
winWidth=10000
tiles<-unlist(tileGenome(seqlengths(Celegans),tilewidth=tileWidth))
seqinfo(tiles)<-seqinfo(Celegans)

tile1000<-trim(resize(tiles,width=winWidth,fix="center"))

tile1000<-binnedAverage(tile1000,cov,"COH1_r1")

tiles$score<-tile1000$COH1_r1
export.bw(tiles,paste0("coh1_r1_tile",tileWidth,"_smootheWin",winWidth,".bw"))


#https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data/54507329#54507329
ThresholdingAlgo <- function(y,lag,threshold,influence) {
  signals <- rep(0,length(y))
  filteredY <- y[0:lag]
  avgFilter <- NULL
  stdFilter <- NULL
  avgFilter[lag] <- mean(y[0:lag], na.rm=TRUE)
  stdFilter[lag] <- sd(y[0:lag], na.rm=TRUE)
  for (i in (lag+1):length(y)){
    if (abs(y[i]-avgFilter[i-1]) > threshold*stdFilter[i-1]) {
      if (y[i] > avgFilter[i-1]) {
        signals[i] <- 1;
      } else {
        signals[i] <- -1;
      }
      filteredY[i] <- influence*y[i]+(1-influence)*filteredY[i-1]
    } else {
      signals[i] <- 0
      filteredY[i] <- y[i]
    }
    avgFilter[i] <- mean(filteredY[(i-lag):i], na.rm=TRUE)
    stdFilter[i] <- sd(filteredY[(i-lag):i], na.rm=TRUE)
  }
  return(list("signals"=signals,"avgFilter"=avgFilter,"stdFilter"=stdFilter))
}

lag=50 # amount of historical data to take into account
thresh=2.5 # the z score at which a peak is called (difference from moving average)
infl=0  # influence of new signals on moving
hist(tiles$score)
abline(v=mean(tiles$score),col="red")
abline(v=c(mean(tiles$score)+sd(tiles$score)*thresh),col="red",lty=2)
pks<-ThresholdingAlgo(tiles$score,lag=lag,threshold=thresh,influence=infl)
idx<-tiles$score>pks$avgFilter+pks$stdFilter
idx[is.na(idx)]<-F
pkgr<-reduce(tiles[idx])
export(pkgr,paste0("threshPk_lag",lag,"_th",thresh,"_infl",infl,"_tol100w10000.bed"))

#library(PeakSegDP)
#cDPA(tiles$score, weight = rep(1, length(tiles$score)), maxSegments)
#requires poisson data. not clear what maxSegments should be

