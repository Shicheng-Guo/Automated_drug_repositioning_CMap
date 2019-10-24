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
## THIS SCRIPT WILL CREATE DATABASE
##  
#####

##### library #####

library(foreach)
library(doParallel)

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### SET DIRECTORIES ####

work_dir="/media/hamel/Workspace/Dropbox/Projects/CMap2/3.Method_development/2.combined_drug_bi_directional_enrichment/"

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

#subset CMap data to exclude drugs with no affect

subset_CMap<-subset(HT_HG_U133A_MCF7, !(names(HT_HG_U133A_MCF7)%in%drug_zero))

##### CREATE CLUSTER #####
#ceate cluster for multi-threading
#Calculate the number of cores
no_cores <- detectCores() - 6
registerDoParallel(no_cores)

#for windows
#cores<- detectCores(logical=T)
#cl<-makeCluster(cores)
#cl<-makeCluster(22) # each core requires 8gb ram?
#registerDoParallel(cl)


##### CREATE DRUG COMBINATIONS DATABASE #####

head(subset_CMap[[1]])
length(names(subset_CMap))

# all possible combinations of drugs - each possible combination is represented by 2 rows.
all_drug_combinations<-combn(names(subset_CMap),2)
length(all_drug_combinations)
class(all_drug_combinations)
head(all_drug_combinations)[,1:5]

# create empty list
Drug_combination_database<-list()

Drug_combination_database<-foreach(x=1:(dim(all_drug_combinations)[2])) %dopar% {
  # extract drug combination positions
  drug1_position<-which(names(subset_CMap) == all_drug_combinations[1,x])
  drug2_position<-which(names(subset_CMap) == all_drug_combinations[2,x])
  # extract sig genes for each drug
  drug1_sig_up<-rownames(subset(subset_CMap[[drug1_position]], adj.P.Val<=0.05 & logFC > 0))
  drug1_sig_down<-rownames(subset(subset_CMap[[drug1_position]], adj.P.Val<=0.05 & logFC < 0))
  drug2_sig_up<-rownames(subset(subset_CMap[[drug2_position]], adj.P.Val<=0.05 & logFC > 0))
  drug2_sig_down<-rownames(subset(subset_CMap[[drug2_position]], adj.P.Val<=0.05 & logFC < 0))
  # remove genes which are reversed by drug i.e up-regulated by 1st drug but down regulated by 2nd drug
  drug1_sig_up2<-subset(drug1_sig_up, !(drug1_sig_up%in%drug2_sig_down))
  drug1_sig_down2<-subset(drug1_sig_down, !(drug1_sig_down%in%drug2_sig_up))
  drug2_sig_up2<-subset(drug2_sig_up, !(drug2_sig_up%in%drug1_sig_down))
  drug2_sig_down2<-subset(drug2_sig_down, !(drug2_sig_down%in%drug1_sig_up))
  # add "UP" + "DOWN naming
  drug1_sig_up2<-paste(drug1_sig_up2, "DOWN", sep="_")
  drug1_sig_down2<-paste(drug1_sig_down2, "UP", sep="_")
  drug2_sig_up2<-paste(drug2_sig_up2, "DOWN", sep="_")
  drug2_sig_down2<-paste(drug2_sig_down2, "UP", sep="_")
  # combine genes
  combined_genes<-c(drug1_sig_up2, drug1_sig_down2, drug2_sig_up2, drug2_sig_down2)
  return(combined_genes)
}

#add drug name

for (x in 1:dim(all_drug_combinations)[2]) {
  names(Drug_combination_database)[x]<-paste(all_drug_combinations[1,x], all_drug_combinations[2,x], sep="_+_")
}

head(names(Drug_combination_database))
length(names(Drug_combination_database))

# save Drug combination database

setwd(work_dir)

save(Drug_combination_database, file="Drug_combination_database.Rdata")

##### STOP CLUSTER #####

stopCluster(cl)

##### SAVE #####

#setwd(work_dir)
#save.image("Query_CMap_bi-directional_with_drug_combination.Rdata")
