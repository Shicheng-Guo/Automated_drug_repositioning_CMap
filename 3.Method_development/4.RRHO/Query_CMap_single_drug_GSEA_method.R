##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                                  QUERY CMAP DATABASE - RRHO                                                            #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################

## AUTHOR: HAMEL PATEL
## EMAIL: Hamel.patel@kcl.ac.uk
## DATE: 21/11/2016

##### DESCRIPTION OF ANALYSIS ####
## MCF7 cell line across 3 chips has been processed and DE analysis performed.
## WORKING WITH FULL PROBE DATASETS - 12443
## THIS SCRIPT WILL RUN Rank-Rank Hypergeometric test USING USER DEFINED DEG LIST IN CMAP
#####

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### LIBRARY #####

library(RRHO)
library(foreach)
library(doParallel)

##### SET DIRECTORIES ####

work_dir="/media/hamel/Workspace/Dropbox/Projects/Brain_expression/7.Query_CMap/7.4.single_drug_GSEA_method"

setwd(work_dir)

# Cmap directory
HG_U133A_MCF7_dir<-"/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/CMAP/HG_U133A_MCF7_processing/DE_results"
HT_HG_U133A_MCF7_dir<-"/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/CMAP/HT_HG_U133A_MCF7_processing/DE_results"
HT_HG_U133A_EA_MCF7_dir<-"/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/CMAP/HT_HG_U133A_EA_MCF7_processing/DE_results"

# AD DEG dir
DEG_dir<-"/media/hamel/Workspace/Dropbox/Projects/Brain_expression/6.DE_Meta_Analysis/"

##### LOAD CMAP DATA ####

setwd(HG_U133A_MCF7_dir)
load("full_DE_results.Rdata")
HG_U133A_MCF7<-full_DE_results

setwd(HT_HG_U133A_MCF7_dir)
load("full_DE_results.Rdata")
HT_HG_U133A_MCF7<-full_DE_results

setwd(HT_HG_U133A_EA_MCF7_dir)
load("full_DE_results.Rdata")
HT_HG_U133A_EA_MCF7<-full_DE_results

rm(full_DE_results)

##### CHECK CMAP DATA #####

dim(HT_HG_U133A_EA_MCF7[[1]])
dim(HT_HG_U133A_MCF7[[1]])
dim(HG_U133A_MCF7[[1]])

names(HT_HG_U133A_EA_MCF7)
names(HT_HG_U133A_MCF7)
names(HG_U133A_MCF7)

##### LOAD AD DEG RESULTS ######

setwd(DEG_dir)

Cerebellum_DEG<-read.table("Cerebellum_results/Cerebellum_AW.OC_DEG.txt", head=T, as.is=T)

Frontal_lobe_DEG<-read.table("Frontal_Lobe_results/Frontal_lobe_AW.OC_DEG.txt", head=T, as.is=T)

Parietal_lobe_DEG<-read.table("Parietal_Lobe_results/Parietal_lobe_AW.OC_DEG.txt", head=T, as.is=T)

Temporal_lobe_DEG<-read.table("Temporal_Lobe_results/Temporal_lobe_AW.OC_DEG.txt", head=T, as.is=T)

##### CHECK AD DEG RESULTS ######

head(Cerebellum_DEG)
head(Frontal_lobe_DEG)
head(Parietal_lobe_DEG)
head(Temporal_lobe_DEG)

##### RANK QUERY SIGNATURE #####

#create ranking function for query signature

rank_data<-function(dataset){
  #susbet query to only genes available in CMap
  dataset<-subset(dataset, (dataset$Entrez_Gene_ID %in% rownames(full_DE_results[[1]]))==T)
  #convert up/down to +1 and -1
  dataset[dataset$Regulation=="up",2]<-as.numeric(1)
  dataset[dataset$Regulation=="down",2]<-as.numeric(-1)
  # change column to numeric
  dataset$Regulation<-as.numeric(dataset$Regulation)
  #subset to up/down
  datasetup<-subset(dataset, Regulation=="1")
  datasetdown<-subset(dataset, Regulation=="-1")
  #rank - most sig up genes ranked highest
  datasetup$rank<-rank(1-datasetup$AW_FDR_adjusted_p_val)
  #rank - most sig down genes ranked lowest (using negative ranking)
  datasetdown$rank<--rank(1-datasetdown$AW_FDR_adjusted_p_val)
  # merge up/down list
  datasetmerge<-rbind(datasetup, datasetdown)
  # assign new rank - highest rank is most sig up gene, lowest rank is sig down gene, middle is where genes not changed much
  datasetmerge$true_rank<-rank(datasetmerge$rank, ties.method = "min")
  #rearrgae by ture rank
  datasetmerge<-datasetmerge[order(-datasetmerge$true_rank),]
  #remove unwanted columns
  datasetmerge$Regulation<-NULL
  datasetmerge$AW_FDR_adjusted_p_val<-NULL
  datasetmerge$rank<-NULL
  #change colname
  colnames(datasetmerge)[2]<-"Rank"
  #asign rownames
  rownames(datasetmerge)<-1:nrow(datasetmerge)
  return(datasetmerge)
}

#apply function
Temporal_lobe_DEG_ranked<-rank_data(Temporal_lobe_DEG)
head(Temporal_lobe_DEG_ranked)
dim(Temporal_lobe_DEG)
summary(Temporal_lobe_DEG_ranked$Rank)


Parietal_lobe_DEG_ranked<-rank_data(Parietal_lobe_DEG)
head(Parietal_lobe_DEG_ranked)
dim(Parietal_lobe_DEG)
summary(Parietal_lobe_DEG_ranked$Rank)


Frontal_lobe_DEG_ranked<-rank_data(Frontal_lobe_DEG)
head(Frontal_lobe_DEG_ranked)
dim(Frontal_lobe_DEG)
summary(Frontal_lobe_DEG_ranked$Rank)


Cerebellum_DEG_ranked<-rank_data(Cerebellum_DEG)
head(Cerebellum_DEG_ranked)
dim(Cerebellum_DEG)
summary(Cerebellum_DEG_ranked$Rank)

##### CREATE CLUSTER FOR PARALLEL PROCESSING #####

#no_cores <- detectCores() - 2
registerDoParallel(15)

##### FRONTAL LOBE RRHO RESULTS #####

#run RRHO in parallel
Frontal_lobe_drug_results<-foreach (x=1:length(names(full_DE_results)),
         .verbose=T,
         .packages="RRHO", 
         .combine=rbind) %dopar% {
  dataset<-subset(full_DE_results[[x]], rownames(full_DE_results[[1]]) %in% Frontal_lobe_DEG_ranked$Entrez_Gene_ID)
  # keep only logFC and P.Value columns
  dataset<-dataset[grep("logFC|adj.P.Val", colnames(dataset))]
  #subset to up/down
  datasetup<-subset(dataset, logFC>0)
  datasetdown<-subset(dataset, logFC<0)
  # reverse rank to query - give up regulated genes negative value and down regulated positive rank
  datasetdown$rank<-rank(1-datasetdown$adj.P.Val)
  #rank - most sig down genes ranked lowest (using negative ranking)
  datasetup$rank<--rank(1-datasetup$adj.P.Val)
  # merge up/down list
  datasetmerge<-rbind(datasetup, datasetdown)
  # assign new rank - highest rank is most sig up gene, lowest rank is sig down gene, middle is where genes not changed much
  datasetmerge$true_rank<-rank(datasetmerge$rank, ties.method = "min")
  #rearrgae by ture rank
  datasetmerge<-datasetmerge[order(-datasetmerge$true_rank),]
  #remove unwanted columns
  datasetmerge$logFC<-NULL
  datasetmerge$adj.P.Val<-NULL
  datasetmerge$rank<-NULL
  #add Entrez_id
  datasetmerge$Entrez_Gene_ID<-rownames(datasetmerge)
  #rearrange columns
  datasetmerge<-datasetmerge[c(2,1)]
  #change colname
  colnames(datasetmerge)[2]<-"Rank"
  #asign rownames
  rownames(datasetmerge)<-1:nrow(datasetmerge)
  #RRHO object
  RRHO.test <-  RRHO(datasetmerge, Frontal_lobe_DEG_ranked, BY=TRUE, alternative='enrichment')
  RRHO_results<-pvalRRHO(RRHO.test, 1000)
  return(RRHO_results$pval)
  gc()
}


# create results table- and arrange by p_val
Frontal_lobe_drug_results_table<-as.data.frame(unlist(Frontal_lobe_drug_results))
Frontal_lobe_drug_results_table$Drug_names<-names(full_DE_results)
colnames(Frontal_lobe_drug_results_table)[1]<-"P_Val"
Frontal_lobe_drug_results_table<-Frontal_lobe_drug_results_table[order(Frontal_lobe_drug_results_table$P_Val),]
rownames(Frontal_lobe_drug_results_table)<-1:nrow(Frontal_lobe_drug_results_table)
Frontal_lobe_drug_results_table<-Frontal_lobe_drug_results_table[c(2,1)]

#check results
head(Frontal_lobe_drug_results_table)
Frontal_lobe_drug_results_table[grep("memantine|galantamine|trichostatin", Frontal_lobe_drug_results_table$Drug_names),]

#write results
setwd(work_dir)
write.table(Frontal_lobe_drug_results_table, "Frontal_lobe_RRHO_results_table.txt", quote=F, sep=",")

##### TEMPORAL LOBE RRHO RESULTS #####

#run RRHO in parallel
Temporal_lobe_drug_results<-foreach (x=1:length(names(full_DE_results)),
                                     .verbose=T,
                                     .packages="RRHO", 
                                     .combine=rbind) %dopar% {
                                       dataset<-subset(full_DE_results[[x]], rownames(full_DE_results[[1]]) %in% Temporal_lobe_DEG_ranked$Entrez_Gene_ID)
                                       # keep only logFC and P.Value columns
                                       dataset<-dataset[grep("logFC|adj.P.Val", colnames(dataset))]
                                       #subset to up/down
                                       datasetup<-subset(dataset, logFC>0)
                                       datasetdown<-subset(dataset, logFC<0)
                                       # reverse rank to query - give up regulated genes negative value and down regulated positive rank
                                       datasetdown$rank<-rank(1-datasetdown$adj.P.Val)
                                       #rank - most sig down genes ranked lowest (using negative ranking)
                                       datasetup$rank<--rank(1-datasetup$adj.P.Val)
                                       # merge up/down list
                                       datasetmerge<-rbind(datasetup, datasetdown)
                                       # assign new rank - highest rank is most sig up gene, lowest rank is sig down gene, middle is where genes not changed much
                                       datasetmerge$true_rank<-rank(datasetmerge$rank, ties.method = "min")
                                       #rearrgae by ture rank
                                       datasetmerge<-datasetmerge[order(-datasetmerge$true_rank),]
                                       #remove unwanted columns
                                       datasetmerge$logFC<-NULL
                                       datasetmerge$adj.P.Val<-NULL
                                       datasetmerge$rank<-NULL
                                       #add Entrez_id
                                       datasetmerge$Entrez_Gene_ID<-rownames(datasetmerge)
                                       #rearrange columns
                                       datasetmerge<-datasetmerge[c(2,1)]
                                       #change colname
                                       colnames(datasetmerge)[2]<-"Rank"
                                       #asign rownames
                                       rownames(datasetmerge)<-1:nrow(datasetmerge)
                                       #RRHO object
                                       RRHO.test <-  RRHO(datasetmerge, Temporal_lobe_DEG_ranked, BY=TRUE, alternative='enrichment')
                                       RRHO_results<-pvalRRHO(RRHO.test, 1000)
                                       return(RRHO_results$pval)
                                       gc()
                                     }


# create results table- and arrange by p_val
Temporal_lobe_drug_results_table<-as.data.frame(unlist(Temporal_lobe_drug_results))
Temporal_lobe_drug_results_table$Drug_names<-names(full_DE_results)
colnames(Temporal_lobe_drug_results_table)[1]<-"P_Val"
Temporal_lobe_drug_results_table<-Temporal_lobe_drug_results_table[order(Temporal_lobe_drug_results_table$P_Val),]
rownames(Temporal_lobe_drug_results_table)<-1:nrow(Temporal_lobe_drug_results_table)
Temporal_lobe_drug_results_table<-Temporal_lobe_drug_results_table[c(2,1)]

#check results
head(Temporal_lobe_drug_results_table)
Temporal_lobe_drug_results_table[grep("memantine|galantamine|trichostatin", Temporal_lobe_drug_results_table$Drug_names),]

#write results
setwd(work_dir)
write.table(Temporal_lobe_drug_results_table, "Temporal_lobe_RRHO_results_table.txt", quote=F, sep=",")

##### PARIETAL LOBE RRHO RESULTS #####

#run RRHO in parallel
Parietal_lobe_drug_results<-foreach (x=1:length(names(full_DE_results)),
                                     .verbose=T,
                                     .packages="RRHO", 
                                     .combine=rbind) %dopar% {
                                       dataset<-subset(full_DE_results[[x]], rownames(full_DE_results[[1]]) %in% Parietal_lobe_DEG_ranked$Entrez_Gene_ID)
                                       # keep only logFC and P.Value columns
                                       dataset<-dataset[grep("logFC|adj.P.Val", colnames(dataset))]
                                       #subset to up/down
                                       datasetup<-subset(dataset, logFC>0)
                                       datasetdown<-subset(dataset, logFC<0)
                                       # reverse rank to query - give up regulated genes negative value and down regulated positive rank
                                       datasetdown$rank<-rank(1-datasetdown$adj.P.Val)
                                       #rank - most sig down genes ranked lowest (using negative ranking)
                                       datasetup$rank<--rank(1-datasetup$adj.P.Val)
                                       # merge up/down list
                                       datasetmerge<-rbind(datasetup, datasetdown)
                                       # assign new rank - highest rank is most sig up gene, lowest rank is sig down gene, middle is where genes not changed much
                                       datasetmerge$true_rank<-rank(datasetmerge$rank, ties.method = "min")
                                       #rearrgae by ture rank
                                       datasetmerge<-datasetmerge[order(-datasetmerge$true_rank),]
                                       #remove unwanted columns
                                       datasetmerge$logFC<-NULL
                                       datasetmerge$adj.P.Val<-NULL
                                       datasetmerge$rank<-NULL
                                       #add Entrez_id
                                       datasetmerge$Entrez_Gene_ID<-rownames(datasetmerge)
                                       #rearrange columns
                                       datasetmerge<-datasetmerge[c(2,1)]
                                       #change colname
                                       colnames(datasetmerge)[2]<-"Rank"
                                       #asign rownames
                                       rownames(datasetmerge)<-1:nrow(datasetmerge)
                                       #RRHO object
                                       RRHO.test <-  RRHO(datasetmerge, Parietal_lobe_DEG_ranked, BY=TRUE, alternative='enrichment')
                                       RRHO_results<-pvalRRHO(RRHO.test, 1000)
                                       return(RRHO_results$pval)
                                       gc()
                                     }


# create results table- and arrange by p_val
Parietal_lobe_drug_results_table<-as.data.frame(unlist(Parietal_lobe_drug_results))
Parietal_lobe_drug_results_table$Drug_names<-names(full_DE_results)
colnames(Parietal_lobe_drug_results_table)[1]<-"P_Val"
Parietal_lobe_drug_results_table<-Parietal_lobe_drug_results_table[order(Parietal_lobe_drug_results_table$P_Val),]
rownames(Parietal_lobe_drug_results_table)<-1:nrow(Parietal_lobe_drug_results_table)
Parietal_lobe_drug_results_table<-Parietal_lobe_drug_results_table[c(2,1)]

#check results
head(Parietal_lobe_drug_results_table)
Parietal_lobe_drug_results_table[grep("memantine|galantamine|trichostatin", Parietal_lobe_drug_results_table$Drug_names),]

#write results
setwd(work_dir)
write.table(Parietal_lobe_drug_results_table, "Parietal_lobe_RRHO_results_table.txt", quote=F, sep=",")

##### CEREBELLUM LOBE RRHO RESULTS #####

#run RRHO in parallel
Cerebellum_drug_results<-foreach (x=1:length(names(full_DE_results)),
                                  .verbose=T,
                                  .packages="RRHO", 
                                  .combine=rbind) %dopar% {
                                    dataset<-subset(full_DE_results[[x]], rownames(full_DE_results[[1]]) %in% Cerebellum_DEG_ranked$Entrez_Gene_ID)
                                    # keep only logFC and P.Value columns
                                    dataset<-dataset[grep("logFC|adj.P.Val", colnames(dataset))]
                                    #subset to up/down
                                    datasetup<-subset(dataset, logFC>0)
                                    datasetdown<-subset(dataset, logFC<0)
                                    # reverse rank to query - give up regulated genes negative value and down regulated positive rank
                                    datasetdown$rank<-rank(1-datasetdown$adj.P.Val)
                                    #rank - most sig down genes ranked lowest (using negative ranking)
                                    datasetup$rank<--rank(1-datasetup$adj.P.Val)
                                    # merge up/down list
                                    datasetmerge<-rbind(datasetup, datasetdown)
                                    # assign new rank - highest rank is most sig up gene, lowest rank is sig down gene, middle is where genes not changed much
                                    datasetmerge$true_rank<-rank(datasetmerge$rank, ties.method = "min")
                                    #rearrgae by ture rank
                                    datasetmerge<-datasetmerge[order(-datasetmerge$true_rank),]
                                    #remove unwanted columns
                                    datasetmerge$logFC<-NULL
                                    datasetmerge$adj.P.Val<-NULL
                                    datasetmerge$rank<-NULL
                                    #add Entrez_id
                                    datasetmerge$Entrez_Gene_ID<-rownames(datasetmerge)
                                    #rearrange columns
                                    datasetmerge<-datasetmerge[c(2,1)]
                                    #change colname
                                    colnames(datasetmerge)[2]<-"Rank"
                                    #asign rownames
                                    rownames(datasetmerge)<-1:nrow(datasetmerge)
                                    #RRHO object
                                    RRHO.test <-  RRHO(datasetmerge, Cerebellum_DEG_ranked, BY=TRUE, alternative='enrichment')
                                    RRHO_results<-pvalRRHO(RRHO.test, 1000)
                                    return(RRHO_results$pval)
                                    gc()
                                  }


# create results table- and arrange by p_val
Cerebellum_drug_results_table<-as.data.frame(unlist(Cerebellum_drug_results))
Cerebellum_drug_results_table$Drug_names<-names(full_DE_results)
colnames(Cerebellum_drug_results_table)[1]<-"P_Val"
Cerebellum_drug_results_table<-Cerebellum_drug_results_table[order(Cerebellum_drug_results_table$P_Val),]
rownames(Cerebellum_drug_results_table)<-1:nrow(Cerebellum_drug_results_table)
Cerebellum_drug_results_table<-Cerebellum_drug_results_table[c(2,1)]

#check results
head(Cerebellum_drug_results_table)
Cerebellum_drug_results_table[grep("memantine|galantamine|trichostatin", Cerebellum_drug_results_table$Drug_names),]

#write results
setwd(work_dir)
write.table(Cerebellum_drug_results_table, "Cerebellum_RRHO_results_table.txt", quote=F, sep=",")

##### CLOSE CLUTER #####

stopImplicitCluster()

##### SAVE #####

setwd(work_dir)

save.image("CMap_query_RRHO.Rdata")
