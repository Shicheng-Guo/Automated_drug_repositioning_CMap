##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                                         QUERY CMAP DATABASE                                                            #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################

## AUTHOR: HAMEL PATEL
## EMAIL: Hamel.patel@kcl.ac.uk
## DATE: 21/11/2016

##### DESCRIPTION OF ANALYSIS ####
## MCF7 cell line across 3 chips has been processed and DE analysis performed.
## WORKING WITH FULL PROBE DATASETS - 12500
## THIS SCRIPT WILL CREATE DATABASE
#####

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### LIBRARY #####

library(WGCNA)

##### SET DIRECTORIES ####

work_dir="/media/hamel/Workspace/Dropbox/Projects/CMap2/3.Method_development/3.single_drug_anti_correlation/"

setwd(work_dir)

# Cmap directory
HT_HG_U133A_MCF7_dir<-"/media/hamel/Workspace/Dropbox/Projects/CMap2/1.Data/3.HT_HG_U133A_MCF7_processing/DE_results/"

##### LOAD CMAP DATA ####

setwd(HT_HG_U133A_MCF7_dir)
load("full_DE_results.Rdata")
HT_HG_U133A_MCF7<-full_DE_results

rm(full_DE_results)

##### CHECK CMAP DATA #####

dim(HT_HG_U133A_MCF7[[1]])
names(HT_HG_U133A_MCF7)

##### ADD +1 and -1 to CMAP data #####

#subset CMAP data to where at least 1 gene is DE

full_DE_results<-HT_HG_U133A_MCF7

head(full_DE_results[[1]])

length(names(full_DE_results))

for (x in 1:length(names(full_DE_results))) {
  #create FC column
  full_DE_results[[x]]$FC<-full_DE_results[[x]]$logFC
  #change all positive values to 1
  full_DE_results[[x]][full_DE_results[[x]]$FC>0,9]<-1
  #change all negative values to -1
  full_DE_results[[x]][full_DE_results[[x]]$FC<0,9]<--1
}

##### CREATE DATAFRAME OF ALL DRUG FC ######

#empty dataframe

CMap_FC<-as.data.frame(matrix(ncol=1, nrow=nrow(full_DE_results[[1]])))

#assign rownames
rownames(CMap_FC)<-rownames(full_DE_results[[1]])

for (x in 1:length(names(full_DE_results))) {
  #extract FC column
  CMap_FC<-merge(CMap_FC,full_DE_results[[x]][9], by="row.names")
  #remove row.names
  rownames(CMap_FC)<-CMap_FC$Row.names
  CMap_FC$Row.names<-NULL
  #assign drug name
  colnames(CMap_FC)[grep("FC", colnames(CMap_FC))]<-as.character(names(full_DE_results)[x])
}

head(CMap_FC)[1:5]

#remove 1st column
CMap_FC$V1<-NULL

head(CMap_FC)[1:5]

(names(full_DE_results)[1])
(names(full_DE_results)[2])
(names(full_DE_results)[3])
(names(full_DE_results)[4])

#check colnames same

any(colnames(CMap_FC)==names(full_DE_results))==F


##### SAVE #####

setwd(work_dir)

write.table(CMap_FC, file="cmap_query_database_method2.txt")

