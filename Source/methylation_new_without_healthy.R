library(TCGAbiolinks)
BiocManager::install("sesameData")
BiocManager::install("methylKit")
BiocManager::install("sesame")

BiocManager::install("DMRcate")
BiocManager::install("genomatio")

setwd("T:/Dokumente")
test <- read.csv(file='HNSC_aliquot.tsv', sep = '\t', header = TRUE)
test2 <- read.csv(file='LUSC_aliquot.tsv', sep = '\t', header = TRUE)

#larynx
query.met.gbm <- GDCquery(project = "TCGA-HNSC",#legacy=TRUE,
                          data.category = "DNA Methylation",data.type="Methylation Beta Value",
                          platform = "Illumina Human Methylation 450",
                          barcode = test$aliquot_submitter_id)

GDCdownload(query.met.gbm)

met.gbm.450 <- GDCprepare(query = query.met.gbm,save = TRUE,
                          save.filename = "gbmDNAmet450k.rda" ,
                          summarizedExperiment = TRUE)

#lung
query.met.lgg <- GDCquery(project = "TCGA-LUSC",#legacy=TRUE,
                          data.category = "DNA Methylation",data.type="Methylation Beta Value",
                          platform = "Illumina Human Methylation 450",
                          barcode = test2$aliquot_submitter_id)


GDCdownload(query.met.lgg)
met.lgg.450 <- GDCprepare(query = query.met.lgg,save = TRUE,
                          save.filename = "lggDNAmet450k.rda",summarizedExperiment = TRUE)

met.lgg.450 <-get(load("lggDNAmet450k.rda"))
met.gbm.450 <-get(load("gbmDNAmet450k.rda"))

#bind rows

#colnames(colData(met.gbm.450))

library(SummarizedExperiment)
idx <- match(colnames(colData(met.gbm.450)), colnames(colData(met.lgg.450)))
colData(met.lgg.450)<-colData(met.lgg.450)[na.omit(idx)]
#colnames(colData(met.lgg.450))<-colnames(colData(met.lgg.450))[na.omit(idx)]
idx <- match(colnames(colData(met.lgg.450)), colnames(colData(met.gbm.450)))
colData(met.gbm.450)<-colData(met.gbm.450)[na.omit(idx)]
#colnames(colData(met.gbm.450))<-colnames(colData(met.gbm.450))[na.omit(idx)]
#dim(met.lgg.450) 
#met.lgg.450

#met.gbm.450<-met.gbm.450[,met.gbm.450@colData$shortLetterCode=="TP"]
#met.lgg.450<-met.lgg.450[,met.lgg.450@colData$shortLetterCode=="TP"]

met.gbm.lgg <- SummarizedExperiment::cbind(met.lgg.450, met.gbm.450)
#met.gbm.lgg@colData$name

#


#––––––––––––––––––––––––––––
# Obtaining DNA methylation
#––––––––––––––––––––––––––––
#library(TCGAbiolinks)
#library(stringr)
# Samples
#matched_met_exp <- function(project, n = NULL){
  # get primary solid tumor samples: DNA methylation
 # message("Download DNA methylation information")
 # met450k <- GDCquery(project = project,
    #                  data.category = "DNA methylation",
   #                   platform = "Illumina Human Methylation 450",
     #                 legacy = TRUE,
     #                 sample.type = c("Primary solid Tumor"))
 # met450k.tp <-  met450k$results[[1]]$cases
  
  # get primary solid tumor samples: RNAseq
#  message("Download gene expression information")
#  exp <- GDCquery(project = project,
        #          data.category = "Gene expression",
        #          data.type = "Gene expression quantification",platform = "Illumina HiSeq",
         #         file.type =  "results",sample.type = c("Primary solid Tumor"),
         #         legacy = TRUE)
  
#  exp.tp <-  exp$results[[1]]$cases
#  print(exp.tp[1:10])
  # Get patients with samples in both platforms
  #patients <- unique(substr(exp.tp,1,15)[substr(exp.tp,1,12) %in% substr(met450k.tp,1,12)])
#  if(!is.null(n)) patients <- patients[1:n] # get only n samples
#  return(patients)
#}
#lgg.samples <- matched_met_exp("TCGA–LGG", n = 10)
#gbm.samples <- matched_met_exp("TCGA–GBM", n = 10)
#samples <- c(lgg.samples,gbm.samples)


#–––––––––––––––––––––––––––––––––––
# 1 – Methylation
# ––––––––––––––––––––––––––––––––––
# For methylation it is quicker in this case to download the tar.gz file
# and get the samples we want instead of downloading files by files
#query.lgg <- GDCquery(project = "TCGA-LGG",data.category = "DNA methylation",
#platform = "Illumina Human Methylation 450",
#legacy = TRUE, barcode = lgg.samples)
#GDCdownload(query.lgg)
#met.lgg <-GDCprepare(query.lgg, save = FALSE)

#query.gbm <- GDCquery(project = "TCGA-GBM", data.category = "DNA methylation",
#platform = "Illumina Human Methylation 450",
#legacy = TRUE, barcode = gbm.samples)
#GDCdownload(query.gbm)
#met.gbm <- GDCprepare(query.gbm, save = FALSE)
#met <- SummarizedExperiment::cbind(met.lgg, met.gbm)

#––––––––––––––––––––––––––––
# Mean methylation
#––––––––––––––––––––––––––––
# Plot a barplot for the groups in the disease column in the
# summarizedExperiment object

# remove probes with NA (similar to na.omit)
#met.gbm.lgg 
met.gbm.lgg <-subset(met.gbm.lgg ,subset = (rowSums(is.na(assay(met.gbm.lgg))) == 0))

# remove probes in chromossomes X, Y and NA
met.gbm.lgg <- subset(met.gbm.lgg ,subset = !as.character(seqnames(met.gbm.lgg)) %in% c("chrNA","chrX","chrY"))

TCGAvisualize_meanMethylation(met,groupCol = "name",
                              group.legend  = "Groups",
                              filename = "mean_lgg_gbm.png",
                              print.pvalue = TRUE)

# Becareful! Depending on the number of probes and samples this function might take some days.
# To make this example faster we used only the chromosome 9
# This should take some minutes
#met.chr9 <- subset(met,subset = as.character(seqnames(met)) %in% c("chr9"))
#met.chr9 <- met
#cancer.met <- subset(met.chr9,subset = (rowSums(is.na(assay(met.chr9))) == 0)) 
#cancer.met <- subset(met,subset = (rowSums(is.na(assay(met))) == 0)) 
met.chr9 <- TCGAanalyze_DMC( met.gbm.lgg,
                            groupCol = "name", # a column in the colData matrix
                            group1 = "Head and Neck Squamous Cell Carcinoma", # a type of the disease type column
                            group2= "Lung Squamous Cell Carcinoma", # a type of the disease column
                            p.cut = 10^-2,diffmean.cut = 0.2,legend = "State",
                            plot.filename = "LGG_GBM_metvolcano_without_healthy.png",
                            cores = 1 # if set to 1 there will be a progress bar
)






colData(cancer.met)$smoker<-colData(cancer.met)$pack_years_smoked>100
met.chr9 <- TCGAanalyze_DMC(cancer.met,
                            groupCol = "smoker", # a column in the colData matrix
                            group1 = "TRUE", # a type of the disease type column
                            group2= "FALSE", # a type of the disease column
                            p.cut = 10^-2,diffmean.cut = 0.25,legend = "State",
                            plot.filename = "LGG_GBM_metvolcano_smoking.png",
                            cores = 1 # if set to 1 there will be a progress bar
)
nrow(met.chr9)
nrow(ll)
ll<-(met.chr9[met.chr9$p.value.adj.Head.and.Neck.Squamous.Cell.Carcinoma.Lung.Squamous.Cell.Carcinoma<0.01,])
nrow(met.chr9[met.chr9$p.value.adj.Head.and.Neck.Squamous.Cell.Carcinoma.Lung.Squamous.Cell.Carcinoma<0.05,][1])

ll
#annotate gene
ll$gene<-met.gbm.lgg@rowRanges$gene[match(rownames(ll), rownames(met.gbm.lgg))]
ll
write.csv(ll,"cpg_significant_w_gene2.csv")







########################bis hier hin
library(minfi)
library(DMRcate)
colData(cancer.met)$years_smoked[is.na(colData(cancer.met)$years_smoked)]<-0
mdsPlot(assay(cancer.met),
        # sampNames =colData(cancer.met)$barcode, 
        sampGroups =colData(cancer.met)$shortLetterCode ,legendPos = c(9,10.11,12))

mdsPlot(assay(cancer.met),
        # sampNames =colData(cancer.met)$barcode, 
        sampGroups =colData(cancer.met)$years_smoked ,legendPos = c(9,10.11,12))
PCASamples(met.chr9, screeplot=FALSE, adj.lim=c(0.0004,0.1), scale=TRUE,
           center=TRUE,comp=c(1,2),transpose=TRUE,sd.filter=TRUE,
           sd.threshold=0.5,filterByQuantile=TRUE,obj.return=FALSE,chunk.size)

cancer.met@rowRanges
aggregateAcrossFeatures(assay(cancer.met), cancer.met@colData$name)
geneBody <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
geneExons <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene,by="gene")
myGeneCounts <- summarizeOverlaps(geneExons,cancer.met@rowRanges,mode="Union", singleEnd=FALSE, fragments=TRUE, )
Age.math <- cancer.met %>% group_by(gene_HGNC) %>% summarize(
  Mean = mean(Age),
  Max = max(Age),
  Min = min(Age),
  sd = sd(Age))

#––––––––––––––––––––––––––
# DNA methylation heatmap
#–––––––––––––––––––––––––
library(ComplexHeatmap)
clin.gbm <- GDCquery_clinic("TCGA-GBM", "Clinical")
clin.lgg <- GDCquery_clinic("TCGA-LGG", "Clinical")
clinical <- plyr::rbind.fill(clin.lgg, clin.gbm)

# get the probes that are Hypermethylated or Hypomethylated
# met is the same object of the section 'DNA methylation analysis'
sig.met <- met.chr9[values(met.chr9)[,"status.Glioblastoma.Multiforme.Brain.Lower.Grade.Glioma"] %in%
                      c("Hypermethylated","Hypomethylated"),]

# To speed up the example, we will select no more than 100 probes
nb.probes <- ifelse(nrow(sig.met) > 100, 100, nrow(sig.met)) # If there is more than 100 get 100
sig.met.100 <- sig.met[1:nb.probes,]

# top annotation, which sampples are LGG and GBM
# We will add clinical data as annotation of the samples
# we will sort the clinical data to have the same order of the DNA methylation matrix
clinical.order <- clinical[match(substr(colnames(sig.met.100),1,12),clinical$bcr_patient_barcode),]
ta = HeatmapAnnotation(df = clinical.order[,c("disease","gender","vital_status","race")],
                       col = list(disease = c("LGG" = "grey", "GBM" = "black"),
                                  gender = c("male"="blue","female"="pink")))
# row annotation: add the status for LGG in relation to GBM
# For exmaple: status.gbm.lgg Hypomethyated means that the
# mean DNA methylation of probes for lgg are hypomethylated
# compared to GBM ones.
ra = rowAnnotation(df = values(sig.met.100)["status.Glioblastoma.Multiforme.Brain.Lower.Grade.Glioma"],
                   col = list("status.Glioblastoma.Multiforme.Brain.Lower.Grade.Glioma" = 
                                c("Hypomethylated" = "orange",
                                  "Hypermethylated" = "darkgreen")),
                   width = unit(1, "cm"))

heatmap <- Heatmap(assay(sig.met.100), name = "DNA methylation",
                   col = matlab::jet.colors(200),
                   show_row_names = F,
                   cluster_rows = T,
                   cluster_columns = F,
                   show_column_names = F,
                   bottom_annotation = ta,
                   column_title = "DNA methylation") 
# Save to pdf
pdf("heatmap.pdf",width = 10, height = 8)
draw(heatmap, annotation_legend_side =  "bottom")
dev.off()

ww<-as.data.frame(assay(cancer.met))

met.chr9<-as.data.frame(read.csv(sep = ",",header=TRUE,"DMR_results_name_Head.and.Neck.Squamous.Cell.Carcinoma_Lung.Squamous.Cell.Carcinoma_pcut_0.01_meancut_0.25_without_healthy.csv"))
write.csv(met.chr9,"met.csv", row.names = FALSE)

dataDEGsFiltLevel<-as.data.frame(read.csv(sep = ",",header=TRUE,"dataDEGsFiltLevel.csv"))


#------------------- Starburst plot ------------------------------
starburst <-TCGAvisualize_starburst(met.chr9,    # DNA methylation with results
                                    dataDEGs,    # DEG results
                                    group1 = "HNSC",
                                    group2 = "LUSC",
                                    filename = "starburst.png",
                                    met.p.cut = 10^-2,
                                    exp.p.cut = 10^-2,
                                    genome="hg38",
                                    met.platform="Illumina Human Methylation 450",
                                    diffmean.cut = 0.25,
                                    logFC.cut = 1,width = 15,height = 10,
                                    names = TRUE)


-----------------------
  
  
  1 library(summarizedExperiment)
2 # get expression matrix
3 data <– assay(exp.gbm.lgg)
4
5 # get genes information
6 genes.info <– rowRanges(exp.gbm.lgg)
7
8 # get sample information
9 sample.info <– colData(exp.gbm.lgg)