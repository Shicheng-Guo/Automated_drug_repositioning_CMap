##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                               LIST OF DRUGS TO REMOVE FROM CMAP                                                        #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################

## AUTHOR: HAMEL PATEL
## EMAIL: Hamel.patel@kcl.ac.uk
## DATE: 07/02/2017

##### DESCRIPTION OF ANALYSIS ####
## CMAP CONTAINS COMPOUNDS WITHRAWN BY FDA FOR VARIOUS REASONS INCLUDING SEVERE SIDE-EFFECTS.
## THIS WILL CREATE A LISAT TO REMOVE
## WITHDARN LIST DOWNLOADED FROM DRUG BANK ON 07/02/2017
#####

##### LIBRARY #####

library(reshape2)

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### SET DIRECTORIES ####

work_dir="/media/hamel/Workspace/Dropbox/Projects/CMap2/1.Data/5.Drug_bank_withdrawn_list/"

setwd(work_dir)

##### LOAD LIST OF DRUGS IN CMAP #####

cmap_pheno<-read.table("/media/hamel/Workspace/Dropbox/Projects/CMap2/1.Data/3.HT_HG_U133A_MCF7_processing/Clean_data/subset_pheno_info.txt",
                       head=T,
                       sep=",", 
                       as.is=T)


head(cmap_pheno)
dim(cmap_pheno)

##### LOAD WITDRAWN DRUG LIST ######

setwd(work_dir)

withdrawn_drugs<-read.csv("Drug_Bank_Withrawn_list.csv",
                          head=T,
                          as.is=T)

head(withdrawn_drugs)

##### CREATE LIST OF DRUGS TO REMOVE FROM CMAP ######

# create cmap drug name list
drugs_to_remove<-cmap_pheno[grep(paste(withdrawn_drugs$Name, collapse = "|"), cmap_pheno$drug_conc, ignore.case=TRUE),]
drugs_to_remove<-drugs_to_remove$drug_conc

drugs_to_remove

##### WRITE OUT FILE OF DRUGS TO REMOVE ######

setwd(work)
write(drugs_to_remove, file="CMap_drug_withdrawn.txt")
