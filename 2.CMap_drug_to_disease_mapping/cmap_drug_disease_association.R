##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                                  CMAP - LINK DRUGS TO DISEASES                                                         #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################

## AUTHOR: HAMEL PATEL
## EMAIL: Hamel.patel@kcl.ac.uk
## DATE: 04/01/2017

##### DESCRIPTION OF ANALYSIS ####
## MAPPING DRUGS IN CMAP VERSION 2 TO DISEASES THAT THE DRUG IS CURRENTLY PRESCRIBED FOR
## CMAP DRUG NAMES SUBSET TO UNIQUE (1310 DRUGS)
## DRUG-DISEASE DATABASE (31614 DRUG-DISEASE) DOWNLOADED FROM THERAPEUTIC TARGETS AND DISEASE (TTD) ON 03/01/2017
## 
#####

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### LIBRARY #####

library(splitstackshape)
library(ggplot2)

##### SET DIRECTORIES ####

work_dir="/media/hamel/Workspace/Dropbox/Projects/CMap2/2.CMap_drug_disease_association"

setwd(work_dir)

##### LOAD DATA #####

#load cmap drug list
#cmap_drugs<-read.csv("cmap_unique_drug_names.csv", header=F)
cmap_drugs<-read.csv("/media/hamel/Workspace/Dropbox/Projects/CMap2/1.Data/3.HT_HG_U133A_MCF7_processing/Clean_data/subset_pheno_info.txt")
# change colname to Drug
names(cmap_drugs)
names(cmap_drugs)[3]<-"Drug"

head(cmap_drugs)

#load TTD database

TTD<-read.csv("drug-disease_TTD2016.txt", skip = 12, header=T, sep="\t")
head(TTD)
dim(TTD)

##### MERGE CMAP AND TTD DATABASE #####

# keep only LNM and Indication columns - rename
TTD<-TTD[c(2,3)]
names(TTD)<-c("Drug", "Diseases")
head(TTD)

# change Drug column in both cmap and TTD from uppercase to lower case
TTD$Drug<-tolower(TTD$Drug)
cmap_drugs$Drug<-tolower(cmap_drugs$Drug)

#merge

cmap_TTD<-merge(cmap_drugs, TTD, by="Drug")

# keep drug name + disease

cmap_TTD<-cmap_TTD[c(1,11)]

head(cmap_TTD)
dim(cmap_TTD)

cmap_TTD

grep("memantine", cmap_TTD$Drug)

# reshape - collpase multiple disesse to unique rows

cmap_TTD<-cSplit(cmap_TTD, "Diseases", ";", "long")

head(cmap_TTD)
dim(cmap_TTD)
length(unique(sort(cmap_TTD$Drug)))
length(unique(sort(cmap_TTD$Diseases)))

# remove duplicate rows
cmap_TTD<-cmap_TTD[!duplicated(cmap_TTD),]

# create table to count disease

disease_count<-as.data.frame(table(cmap_TTD$Diseases))
#reverse order
disease_count<-disease_count[order(-disease_count$Freq),]

head(disease_count, 30)

#colnames

colnames(disease_count)<-c("Disease", "Number_of_Drugs")

disease_count$Number_of_Drugs<-as.numeric(disease_count$Number_of_Drugs)

#subset to more than 1 drug

disease_count_subset<-subset(disease_count, disease_count$Number_of_Drugs >1)
dim(disease_count_subset)

#plot 
ggplot(data=disease_count_subset, aes(x=reorder(Disease, Number_of_Drugs), y=Number_of_Drugs)) +
  geom_bar(position="dodge",stat="identity") + 
  coord_flip() +
  ggtitle("Treatable diseases by cmap drugs")

#plot

png("Treatable diseases by cmap drugs.png", width = 1000, height =1000, units="px")
ggplot(data=disease_count_subset, aes(x=reorder(Disease, Number_of_Drugs), y=Number_of_Drugs)) +
  geom_bar(position="dodge",stat="identity") + 
  coord_flip() +
  ggtitle("Treatable diseases by cmap drugs") +
  labs(x = "Number of CMap drugs", y = "Disease")
dev.off()

#  geom_bar(stat='identity') +
#  coord_flip()
## AD check

subset(disease_count, Disease=="Alzheimer's disease")

subset(cmap_TTD, Diseases=="Alzheimer's disease")

##### AFTER QC IN CMAP #####

# mapping to HT HG U133A

data_CMap<-read.csv("/media/hamel/Workspace/Dropbox/Projects/CMap2/1.Data/3.HT_HG_U133A_MCF7_processing/DE_results/DEG_count.txt")

head(data_CMap)

colnames(data_CMap)[1]<-c("drug_conc")

library(stringr)
data_CMap$Drug<-print(str_split_fixed(data_CMap$drug_conc, "_", 2)[,1], quote=F)

head(data_CMap)

dim(data_CMap)

# TDD needs to be further split

library(splitstackshape)

dim(TTD)
head(TTD)
tail(TTD)

TTD$Diseases[grep(";", TTD$Diseases)]

TTD2<-cSplit(TTD, "Diseases", direction="long", sep=";")

dim(TTD2)
head(TTD2)
tail(TTD2)

TTD2$Diseases[grep(";", TTD2$Diseases)]

cmap_TTD2<-merge(data_CMap, TTD2, by="Drug")

dim(cmap_TTD2)
length(unique(cmap_TTD2$Diseases))
length(unique(cmap_TTD2$Drug))


# create table to count disease

disease_count2<-as.data.frame(table(cmap_TTD2$Diseases))
#reverse order
disease_count2<-disease_count2[order(-disease_count2$Freq),]

head(disease_count2, 30)

#colnames

colnames(disease_count2)<-c("Disease", "Number_of_Drugs")

disease_count2$Number_of_Drugs<-as.numeric(disease_count2$Number_of_Drugs)

#subset to more than 1 drug

disease_count2_subset<-subset(disease_count2, disease_count2$Number_of_Drugs >5)
dim(disease_count2_subset)

#plot 
ggplot(data=disease_count2_subset, aes(x=reorder(Disease, Number_of_Drugs), y=Number_of_Drugs)) +
  geom_bar(position="dodge",stat="identity") + 
  coord_flip() +
  ggtitle("Treatable diseases by cmap drugs")

#plot

png("Treatable diseases by cmap drugs after QC.png", width = 1000, height =1000, units="px")
ggplot(data=disease_count2_subset, aes(x=reorder(Disease, Number_of_Drugs), y=Number_of_Drugs)) +
  geom_bar(position="dodge",stat="identity") + 
  coord_flip() +
  ggtitle("Treatable diseases by CMap compounds") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y = "Number of CMap compounds", x = "Disease")
dev.off()

#  geom_bar(stat='identity') +
#  coord_flip()
## AD check

subset(disease_count2, Disease=="Alzheimer's disease")

subset(cmap_TTD2, Diseases=="Alzheimer's disease")


##### WRITE ######

setwd(work_dir)
write.table(cmap_TTD2, file="cmap_drug_disease_mapping.txt", sep="\t", row.names=F, quote=F)
#write.table(disease_count, file="number_of_cmap_drugs_per_disease.txt", sep="\t", row.names=F, quote=F)
write.table(subset(disease_count2, disease_count2$Number_of_Drugs>0), file="number_of_cmap_drugs_per_disease_after_QC.txt", sep="\t", row.names=F, quote=F)

###### SAVE IMAGE #####

setwd(work_dir)
#save.image("cmap_drug_disease_association.Rdata")
