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
## WORKING WITH FULL PROBE DATASETS - 12501
## THIS SCRIPT WILL CREATE QUERY DATABASE
#####

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### LIBRARY #####

library(WGCNA)

##### SET DIRECTORIES ####

work_dir="/media/hamel/Workspace/Dropbox/Projects/CMap2/3.Method_development/1.single_drug_bi_directional_enrichment/"

setwd(work_dir)

# Cmap directory
HT_HG_U133A_MCF7_dir<-"/media/hamel/Workspace/Dropbox/Projects/CMap2/1.Data/3.HT_HG_U133A_MCF7_processing/DE_results/"

# AD DEG dir
#DEG_dir<-"/media/hamel/Workspace/Dropbox/Projects/Brain_expression/6.Meta_Analysis/AW_results"

##### LOAD CMAP DATA ####

setwd(HT_HG_U133A_MCF7_dir)
load("full_DE_results.Rdata")
HT_HG_U133A_MCF7<-full_DE_results

##### CHECK CMAP DATA #####

dim(HT_HG_U133A_MCF7[[1]])

names(HT_HG_U133A_MCF7)

# ##### TEST DRUG ####
# 
# setwd(work_dir)
# 
# # test no drug - trichostatin
# head(names(HT_HG_U133A_MCF7))
# 
# test_drug<-HT_HG_U133A_MCF7[grep("trichostatin A_1e-07",names(HT_HG_U133A_MCF7))]
# 
# head(test_drug[[1]])
# names(test_drug)
# 
# # extract full list of genes
# 
# test_drug_full_gene_list<-rownames(test_drug[[1]])
# test_drug_sig_up<-rownames(subset(test_drug[[1]], adj.P.Val<=0.05 & logFC > 0))
# test_drug_sig_down<-rownames(subset(test_drug[[1]], adj.P.Val<=0.05 & logFC < 0))
# 
# # make copy of gene list to convert into background or sig
# 
# categories_up<-test_drug_full_gene_list
# 
# for (x in 1:length(categories_up)){
#   if (categories_up[x]%in%test_drug_sig_up==T) {
#     categories_up[x]<-"sig"
#   }
#   else{
#     categories_up[x]<-"background"
#   }
# }
# 
# 
# # same for down genes
# categories_down<-test_drug_full_gene_list
# 
# for (x in 1:length(categories_down)){
#   if (categories_down[x]%in%test_drug_sig_down==T) {
#     categories_down[x]<-"sig"
#   }
#   else{
#     categories_down[x]<-"background"
#   }
# }
# 
# setwd(work_dir)
# 
# write(c("TESTLIST1",as.character(Temporal_lobe_AW.OC_DEG_down$Entrez_Gene_ID), sep="\n"),"TESTLIST1.txt")
# write(c("TESTLIST2",as.character(Temporal_lobe_AW.OC_DEG_up$Entrez_Gene_ID),sep="\n"),"TESTLIST2.txt")
# 
# 
# 
# #test enrichment
# testResults_down = userListEnrichment(
#   # full gene list CMap - per drug
#   geneR=test_drug_full_gene_list,
#   # assign significant and background genes
#   labelR=categories_up,
#   # query llist - brain region genes
#   fnIn="TESTLIST1.txt")
# 
# # To see a list of all significant enrichments,
# testResults_down$sigOverlaps
# testResults_down$pValues
# testResults_down$ovGenes
# 
# # reverse test
# 
# testResults_up = userListEnrichment(
#   # full gene list CMap - per drug
#   geneR=test_drug_full_gene_list,
#   # assign significant and background genes
#   labelR=categories_down,
#   # query llist - brain region genes
#   fnIn="TESTLIST2.txt")
# 
# # To see a list of all significant enrichments,
# testResults_up$sigOverlaps
# testResults_up$pValues
# testResults_up$ovGenes
# 
##### SUBSET CMAP DATA #####

#subset CMAP data to where at least 1 gene is DE

full_DE_results<-HT_HG_U133A_MCF7

#count number of sig adj p value for each drug
#get drug names
number_of_DEG_in_full_data<-as.data.frame(names(full_DE_results))
#convert to character format
number_of_DEG_in_full_data$`names(full_DE_results)`<-as.character(number_of_DEG_in_full_data$`names(full_DE_results)`)
# rename column
colnames(number_of_DEG_in_full_data)<-"Drug"
number_of_DEG_in_full_data$sig_DEG_count<-as.numeric("0")
head(number_of_DEG_in_full_data)

#add drug count
for (x in 1:dim(number_of_DEG_in_full_data)[1]){
  number_of_DEG_in_full_data[x,2]<-dim(subset(full_DE_results[[x]], adj.P.Val<=0.05))[1]
}


number_of_DEG_in_full_data<-as.data.frame(number_of_DEG_in_full_data[order(number_of_DEG_in_full_data$sig_DEG_count, decreasing = T),])

rownames(number_of_DEG_in_full_data)<-1:dim(number_of_DEG_in_full_data)[1]

head(number_of_DEG_in_full_data)
tail(number_of_DEG_in_full_data)

#list drugs with 0 DE
drug_zero<-(subset(number_of_DEG_in_full_data, sig_DEG_count ==0))$Drug
drug_zero

#subset CMap data

subset_CMap<-subset(HT_HG_U133A_MCF7, !(names(HT_HG_U133A_MCF7)%in%drug_zero))

#### SAVE QUERY CMAP DATABASE ######
setwd(work_dir)

save(subset_CMap, file="cmap_database_for_bi_directional_enrichment_search.Rdata")

##### SAVE #####

setwd(work_dir)

save.image("CMap_query.Rdata")
