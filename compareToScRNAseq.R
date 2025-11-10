library(monocle3)
library(ggplot2)
library(cowplot)
library(readxl)

outPath="."
scriptPath="."
source(paste0(scriptPath,"/functions.R"))

theme_set(
  theme_bw()+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.title.y=ggtext::element_markdown())
)

scriptName="compareToScRNAseq"
fileNamePrefix=paste0(scriptName,"/")
makeDirs(outPath,dirNameList=c(paste0(c("plots/"),scriptName)))

padjVal=0.05
lfcVal=0


#cds_baseline_post_sub.rds
#https://zenodo.org/record/7958249/files/cds_baseline_post_sub.rds?download=1
cds<-readRDS("/Users/semple/Documents/MeisterLab/Datasets/Ghaddar_2023_scRNAseq_yAd/cds_baseline_post_sub.rds")
sigUp<-read.csv("/Users/semple/Documents/MeisterLab/otherPeopleProjects/Bolaji/BolajiRNAseq_20211216/sigUp_DEseq2.csv",header=F)
sigDown<-read.csv("/Users/semple/Documents/MeisterLab/otherPeopleProjects/Bolaji/BolajiRNAseq_20211216/sigDown_DEseq2.csv",header=F)
sigUp$upVdown<-"COH-1_up"
sigDown$upVdown<-"COH-1_down"

smcTbl<-data.frame(gene=c("WBGene00000591","WBGene00004737"),kleisin=c("coh-1","scc-1"))

sig<-rbind(sigUp,sigDown)

#colData(cds)$cell_type_group<-factor(colData(cds)$cell_type_group)
#colData(cds)$cluster<-factor(colData(cds)$cell_type_group)
colData(cds)

p1<-plot_cells(cds,color_cells_by="cell_type_group", cell_size=1,
               label_cell_groups = T,group_cells_by="partition",
               group_label_size=1.5) +
  ggtitle("Tissue type")


p2<-plot_cells(cds,color_cells_by="assigned_cell_type", cell_size=1,
               label_cell_groups = T,group_cells_by="cluster",
               group_label_size=1.2) +
  ggtitle("Cell type")

p3<-plot_cells(cds,genes = sig, label_cell_groups = F, cell_size=1) +
  #scale_color_gradient(low="#FEFEE3",high="#8B1A1A")
  scale_color_gradient(low="white",high="black")

p4<-plot_cells(cds,genes = smcTbl$gene[smcTbl$kleisin=="coh-1"], label_cell_groups = F, cell_size=0.7) +
  #scale_color_gradient(low="#FEFEE3",high="#8B1A1A")
  scale_color_gradient(low="white",high="black")

p5<-plot_cells(cds,genes = smcTbl$gene[smcTbl$kleisin=="scc-1"], label_cell_groups = F, cell_size=0.7) +
  #scale_color_gradient(low="#FEFEE3",high="#8B1A1A")
  scale_color_gradient(low="white",high="black")

p<-plot_grid(plot_grid(p1,p2,ncol=2,labels=c("A ","B ")),p3,nrow=2,labels=c("","C "))
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"COH1cs_upVdown_scRNAseq_Ghaddar2023r1.png"),p,height=21,width=29,units="cm",device="png")

p<-plot_grid(plot_grid(p1,p2,ncol=2,labels=c("A ","B ")),
             plot_grid(p4,p5,ncol=2,labels=c("C ","C ")),
             nrow=2)
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"COH1vSCC1_scRNAseq_Ghaddar2023r1.png"),p,height=21,width=29,units="cm",device="png")




# tissue types - heatmap -----
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)),
                                cell_group=colData(cds)$cell_type_group)

agg_mat <- aggregate_gene_expression(cds, sig, cell_group_df)
#row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
#colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))

p<-pheatmap::pheatmap(t(agg_mat), cluster_rows=T, cluster_cols=F,
                   scale="column", clustering_method="ward.D2",
                   fontsize=10, angle_col=45)
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"COH1cs_upVdown_scRNAseq_heatmapTissues_Ghaddar2023.png"),p,height=11,width=9,units="cm",device="png")

# Closer look at Neurons -----
## neuron types ----
neurons_cds <- cds[,grepl("Neurons", colData(cds)$cell_type_group, ignore.case=TRUE)]
cell_group_df <- tibble::tibble(cell=row.names(colData(neurons_cds)),
                                cell_group=colData(neurons_cds)$assigned_cell_type)

agg_mat <- aggregate_gene_expression(neurons_cds, sig, cell_group_df)
#row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
#colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))

p<-pheatmap::pheatmap(t(agg_mat), cluster_rows=T, cluster_cols=F,
                      scale="column", clustering_method="ward.D2",
                      fontsize=9, angle_col=45)
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"COH1cs_upVdown_scRNAseq_heatmapNeurons_Ghaddar2023.png"),p,height=29,width=9,units="cm",device="png")

## neural groups ----
library(readxl)
library(dplyr)
nt<-read_excel("/Users/semple/Documents/MeisterLab/Datasets/wormatlas/Ce_NTtables_Loer&Rand2022.xlsx",sheet="All neurons, sorted by name",skip=7)
options(dplyr.print_max = 1e9)
#narows<-which(is.na(nt$`Neuron class`))
notnarows<-which(!is.na(nt$`Neuron class`))
#nt$NeuronClass<-nt$`Neuron class`
# for(narow in narows){
#   closest<-max(notnarows[notnarows<narow])
#   nt$NeuronClass[narow]<-nt$`Neuron class`[closest]
# }

nt<-nt[notnarows,]
tail(nt,20) %>% print(width=Inf)

nt$`h / m`[is.na(nt$`h / m`)]<-"h&m"
nt$progenitor<-gsub("([[:space:]]|\\.).*$","",nt$`WA - Lineage`)
nt$finalAxis<-gsub(".*?(a|p|l|r|d|v)$","\\1",gsub("(\\,|\\().*$","",nt$`WA - Lineage`))


### Lineage timings ----
# get timings from http://wormweb.org/celllineage#c=Eprp&z=0.37
timings=data.frame(stage=c("fertilisation","hatch","L1molt","L2molt",
                           "L3molt","L4molt"),
                   time=c(0,800,1760,2300,2840,3500))

ww<-read.delim("/Users/semple/Documents/MeisterLab/Datasets/CelegansLineage/lineage_wormweb.tsv",
                sep="\t",header=F)
names(ww)<-c("did","id","name","levelTime","totalTime","deathTime","type","children","childBirth")
head(ww)
ww$birthStage<-NA
ww$birthStage[ww$totalTime<timings$time[timings$stage=="hatch"]]<-"embryo"
ww$birthStage[ww$totalTime<timings$time[timings$stage=="L1molt"] &
                ww$totalTime>=timings$time[timings$stage=="hatch"]]<-"L1"
ww$birthStage[ww$totalTime<timings$time[timings$stage=="L2molt"] &
                ww$totalTime>=timings$time[timings$stage=="L1molt"]]<-"L2"
ww$birthStage[ww$totalTime<timings$time[timings$stage=="L3molt"] &
                ww$totalTime>=timings$time[timings$stage=="L2molt"]]<-"L3"
ww$birthStage[ww$totalTime<timings$time[timings$stage=="L4molt"] &
                ww$totalTime>=timings$time[timings$stage=="L3molt"]]<-"L4"
table(ww$birthStage)

nt<-left_join(nt,ww,by=join_by("Neuron"=="name"))


cell_group_df1<-left_join(cell_group_df,nt,by=join_by("cell_group"=="Neuron class"),keep=F)
mixedType<-is.na(cell_group_df1$N)

cell_group_df_doublets<-tibble(cell_group_df) %>% tidyr::extract(cell_group, c("A", "B"), "([[:alnum:]]+)_([[:alnum:]]+)") %>%
  dplyr::filter(!is.na(A)) %>% tidyr::pivot_longer(cols=c("A","B"),names_to="mixedCells",values_to="cell_group") %>%
  dplyr::left_join(nt,by=join_by("cell_group"=="Neuron class"),keep=F)%>%
  dplyr::filter(!is.na(Neuron)) %>% select(-mixedCells)

cell_group_df_triplets<-tibble(cell_group_df) %>% tidyr::extract(cell_group, c("A", "B","C"), "([[:alnum:]]+)_([[:alnum:]]+)_([[:alnum:]]+)") %>%
  dplyr::filter(!is.na(A)) %>% tidyr::pivot_longer(cols=c("A","B","C"),names_to="mixedCells",values_to="cell_group") %>%
  dplyr::left_join(nt,by=join_by("cell_group"=="Neuron class"),keep=F) %>%
  dplyr::filter(!is.na(Neuron)) %>% select(-mixedCells)

cell_group_df2<-rbind(cell_group_df1[!mixedType,],cell_group_df_doublets,cell_group_df_triplets)

# mixed neuron types
agg_mat <- aggregate_gene_expression(neurons_cds, sig, cell_group_df1[mixedType,c("cell","cell_group")])

p<-pheatmap::pheatmap(t(agg_mat), cluster_rows=T, cluster_cols=F,
                      scale="column", clustering_method="ward.D2",
                      fontsize=9, angle_col=45)
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"COH1cs_upVdown_scRNAseq_heatmapMixedNeurons_Ghaddar2023.png"),p,height=29,width=9,units="cm",device="png")

# neurontransmitter
agg_mat <- aggregate_gene_expression(neurons_cds, sig, cell_group_df1[!mixedType,c("cell","Neurotransmitter 1")])

p<-pheatmap::pheatmap(t(agg_mat), cluster_rows=T, cluster_cols=F,
                      scale="column", clustering_method="ward.D2",
                      fontsize=9, angle_col=45)
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"COH1cs_upVdown_scRNAseq_heatmapNeuroTransmitter_Ghaddar2023.png"),p,height=29,width=9,units="cm",device="png")
table(nt$`Neurotransmitter 1`)

# neurontransmitter
agg_mat <- aggregate_gene_expression(neurons_cds, sig, cell_group_df2[,c("cell","Neurotransmitter 1")])

p<-pheatmap::pheatmap(t(agg_mat), cluster_rows=T, cluster_cols=F,
                      scale="column", clustering_method="ward.D2",
                      fontsize=9, angle_col=45)
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"COH1cs_upVdown_scRNAseq_heatmapNeuroTransmitter_withAmbiguous_Ghaddar2023.png"),p,height=29,width=9,units="cm",device="png")
table(nt$`Neurotransmitter 1`)

# soma location
agg_mat <- aggregate_gene_expression(neurons_cds, sig, cell_group_df1[!mixedType,c("cell","Soma Location")])

p<-pheatmap::pheatmap(t(agg_mat), cluster_rows=T, cluster_cols=F,
                      scale="column", clustering_method="ward.D2",
                      fontsize=9, angle_col=45)
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"COH1cs_upVdown_scRNAseq_heatmapLocation_Ghaddar2023.png"),p,height=29,width=9,units="cm",device="png")
table(nt$`Soma Location`)

# herm vs male
cell_group_df1$`h / m`[is.na(cell_group_df1$`h / m`)]<-"h&m"
agg_mat <- aggregate_gene_expression(neurons_cds, sig, cell_group_df1[!mixedType,c("cell","h / m")])
p<-pheatmap::pheatmap(t(agg_mat), cluster_rows=T, cluster_cols=F,
                      scale="column", clustering_method="ward.D2",
                      fontsize=9, angle_col=45)
ggsave(filename=paste0(outPath, "/plots/",fileNamePrefix,"COH1cs_upVdown_scRNAseq_hermVmale_Ghaddar2023.png"),p,height=29,width=9,units="cm",device="png")
table(nt$`h / m`)




##########-
## compare to postembryonic neuronal development genes------
##########-
# Temporal transitions in the post-mitotic nervous system of Caenorhabditis elegans. HaoSheng Sun & Oliver Hobert  https://doi.org/10.1038/s41586-021-04071-4
# supplementary table 3
devNeur<-read_excel("~/Downloads/41586_2021_4071_MOESM7_ESM.xlsx",
                    col_types=c("text",rep("numeric",12)))
dim(devNeur)
head(devNeur)
devNeur<-devNeur[-1,]
names(devNeur)[1]<-"publicID"
names(devNeur)<-gsub(" ","",names(devNeur))
names(devNeur)[grep("\\.\\.\\.",names(devNeur))-1]<-paste0("lfc_",names(devNeur)[grep("\\.\\.\\.",names(devNeur))-1])
names(devNeur)[grep("\\.\\.\\.",names(devNeur))]<-gsub("lfc_","padj_",names(devNeur)[grep("\\.\\.\\.",names(devNeur))-1])
devNeur

salmon<-readRDS("./rds/coh1_noOsc/coh1_noOsc_COH1vsTEVonly_DESeq2_fullResults.rds")
salmon$publicID[salmon$publicID==""]<-salmon$sequenceID[salmon$publicID==""]

devNeur1<-dplyr::inner_join(devNeur,data.frame(salmon),by=c("publicID"))

colIdx<-c(grep("lfc_",names(devNeur1)), grep("log2FoldChange",names(devNeur1)))

GGally::ggpairs(devNeur1[colIdx]) + geom_vline(xintercept=0,type=2,alpha=0.3) +
  geom_hline(yintercept=0,type=2,alpha=0.3)

# distance to fountains ----
devNeurSig<-read_excel("~/Downloads/41586_2021_4071_MOESM8_ESM.xlsx")
dim(devNeurSig)
head(devNeurSig)
names(devNeurSig)[1]<-"sequenceID"

devNeurNS<- devNeur[!(devNeur$sequenceID %in% devNeurSig$sequenceID),]

metadata<-readRDS("ce11GeneGR_WS285.rds")

sigGR<-metadata[metadata$sequenceID %in% devNeurSig$sequenceID |
                  metadata$publicID %in% devNeurSig$sequenceID] # only 2612 from 2639

nsGR<-metadata[metadata$sequenceID %in% devNeurNS$sequenceID |
                  metadata$publicID %in% devNeurNS$sequenceID] # only 5281 from 5335

fountains<-readRDS("/Users/semple/Documents/MeisterLab/otherPeopleProjects/fountains/detected_fountains_equalQ.RDS")
fountains$fountainName<-paste0("fount",1:length(fountains))

dfsig<-data.frame(type="sig",distance=mcols(distanceToNearest(sigGR,fountains,ignore.strand=T))$distance)

dfns<-data.frame(type="ns",distance=mcols(distanceToNearest(nsGR,fountains,ignore.strand=T))$distance)

ggplot(rbind(dfsig,dfns),aes(x=type,y=log10(distance+1))) + geom_violin(fill="grey")
wilcox.test(log10(dfsig$distance+1),log10(dfns$distance+1))

ggplot(rbind(dfsig,dfns),aes(color=type,x=log10(distance+1))) + geom_density(linewidth=2)

ggplot(rbind(dfsig,dfns),aes(color=type,x=distance)) + geom_density(linewidth=2)



###-
# check specific genes-----
###-
metadata<-readRDS("/Users/semple/Documents/MeisterLab/otherPeopleProjects/Bolaji/BolajiRNAseq_20211216/ce11GeneGR_WS285.rds")
geneList<-c("C39H7.4","lag-1","cpi-1","srsx-16","skn-1")
matchName<-function(g){
  newg<-ifelse(!is.na(match(g, metadata$publicID)),match(g, metadata$publicID), match(g, metadata$sequenceID))
  return(metadata$wormbaseID[newg])
}
geneListwb<-sapply(geneList, matchName)

ASIandGut<-cds[,(grepl("Intestine", colData(cds)$cell_type_group, ignore.case=TRUE)| grepl("ASI",colData(cds)$assigned_cell_type,ignore.case=TRUE))]

ASIandGutsubset<-ASIandGut[rowData(ASIandGut)$id %in% geneListwb,]

p<-plot_genes_violin(ASIandGutsubset, group_cells_by="assigned_cell_type", ncol=5) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
p+ggtitle("Ghaddar(2023)  scRNAseq")




### neuronal connectivity:
### https://wormwiring.org
### neuronal adjacency:
### https://www.cell.com/current-biology/fulltext/S0960-9822(23)00553-5
