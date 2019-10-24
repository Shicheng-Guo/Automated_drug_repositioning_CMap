
##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                                       CMap build 2.0 PROCESSING - PART 1                                               #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################

# MICROARRAY PLATFORM - Affymetrix
# EXPRESSION CHIP - 
# NUMBER OF SAMPLES - 
# TISSUE - 
#
# NOTES
# using gcRMA instead of MAs5 to background correct data
 # only use MCF7

## AUTHOR: HAMEL PATEL
## EMAIL: Hamel.patel@kcl.ac.uk
## DATE: 20/09/2016

##### DESCRIPTION OF ANALYSIS ####
## CMAP build 2 downloaded from from http://portals.broadinstitute.org/cmap/ on 4/08/2016.
## This script will pre-process the raw data - separate date by chip + tissue and background correct + normalise
#####

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### SET DIRECTORIES ####

data_dir="/media/hamel/Workspace/Dropbox/Projects/CMap/1.Data/CMAP/Raw_Data/"

work_dir="/media/hamel/Workspace/Dropbox/Projects/CMap2/1.Data/"

unzipped_dir="/media/hamel/Workspace/Dropbox/Projects/CMap/1.Data/CMAP/Raw_Data/unzipped/"

setwd(data_dir)

dir.create(paste(work_dir,"background_corrected_and_normalised", sep="/"))
background_corrected_dir=paste(work_dir,"background_corrected_and_normalised", sep="/")

dir.create(paste(work_dir,"normalised", sep="/"))
normalised_dir=paste(work_dir,"normalised", sep="/")

dir.create(paste(work_dir,"Raw_data_R_object", sep="/"))
Raw_data_R_object_dir=paste(work_dir,"Raw_data_R_object", sep="/")

##### LIBRARIES #####

#library(R.utils)
library(affy)
library(lumi)
library(simpleaffy)
library(stringr)
library(gcrma)
library(hthgu133aprobe)
#library(hgu133aprobe)

##### READ IN CMAP DATA INFO FILE ####

# file downloaded from cmap - 
setwd(work_dir)

cmap_data_info<-read.csv("/media/hamel/Workspace/Dropbox/Projects/CMap/1.Data/CMAP/cmap_instances_02.csv",
                         head=T,
                         fill=T,
                         check.names = F)

head(cmap_data_info[1:5])
tail(cmap_data_info)
colnames(cmap_data_info)

dim(cmap_data_info)

#keep 1st 6100 rows
cmap_data_info<-cmap_data_info[1:6100,]

dim(cmap_data_info)

head(cmap_data_info[1:5])

#number of chips
table(cmap_data_info$array3)

# random "'" found in perturbation scan id and vehicle_scan_id4 - remove "'" from original file
grep("'", cmap_data_info$perturbation_scan_id)
grep("'", cmap_data_info$vehicle_scan_id4)

#convert columns to character
cmap_data_info$perturbation_scan_id<-as.character(cmap_data_info$perturbation_scan_id)
cmap_data_info$vehicle_scan_id4<-as.character(cmap_data_info$vehicle_scan_id4)

for (x in 1:dim(cmap_data_info)[1]){
  if (substring(cmap_data_info$perturbation_scan_id[x],1,1)=="'"){
    cmap_data_info$perturbation_scan_id[x]<-substring(cmap_data_info$perturbation_scan_id[x],2)
  }
}

for (x in 1:dim(cmap_data_info)[1]){
  if (substring(cmap_data_info$vehicle_scan_id4[x],1,1)=="'"){
    cmap_data_info$vehicle_scan_id4[x]<-substring(cmap_data_info$vehicle_scan_id4[x],2)
  }
}

grep("'", cmap_data_info$perturbation_scan_id)
grep("'", cmap_data_info$vehicle_scan_id4)

#create a column with cell IDs - append scan ID to ".CELL"
cmap_data_info$cell_id<-paste(cmap_data_info$perturbation_scan_id,".CEL", sep="")

#separete cell ID by chip
colnames(cmap_data_info)

HG_U133A<-as.vector(cmap_data_info[cmap_data_info$array3=="HG-U133A", 16])
length(HG_U133A)

HT_HG_U133A<-cmap_data_info[cmap_data_info$array3=="HT_HG-U133A", 16]
length(HT_HG_U133A)

HT_HG_U133A_EA<-cmap_data_info[cmap_data_info$array3=="HT_HG-U133A_EA", 16]
length(HT_HG_U133A_EA)

# remove characters from batch_id column
cmap_data_info$batch_id<-as.character(cmap_data_info$batch_id)
grep("a|b", cmap_data_info$batch_id)

# loop through batch_id column - if last character == a | b  - remove
for (x in 1:dim(cmap_data_info)[1]){
  # check if last character is a or b
  if (substring(cmap_data_info$batch_id[x],nchar(cmap_data_info$batch_id[x]),)=="a" || substring(cmap_data_info$batch_id[x],nchar(cmap_data_info$batch_id[x]),)=="b" ){
    # remove last charcter if above is true
    cmap_data_info$batch_id[x]<-substring(cmap_data_info$batch_id[x],1,nchar(cmap_data_info$batch_id[x])-1)
  }
}

grep("a|b", cmap_data_info$batch_id)

##### CREATE CONTROL INFO FILE #####

# make copy
cmap_data_info_control<-cmap_data_info

# count number of perioids in "vehicle_scan_id4 column" 

cmap_data_info_control$count<-str_count(cmap_data_info_control$vehicle_scan_id4, fixed('.'))
head(cmap_data_info_control)
tail(cmap_data_info_control)
table(cmap_data_info_control$count)

# split file by count. <3 for those that need .CEL attached to end of name, and that above 3 for those that require aditional cel files

cmap_data_info_control_1<-subset(cmap_data_info_control, count<3)
cmap_data_info_control_5<-subset(cmap_data_info_control, count==5)
cmap_data_info_control_6<-subset(cmap_data_info_control, count==6)

#drop levels
cmap_data_info_control_5<-droplevels(cmap_data_info_control_5)
cmap_data_info_control_6<-droplevels(cmap_data_info_control_6)

table(cmap_data_info_control_5$vehicle_scan_id4)
table(cmap_data_info_control_6$vehicle_scan_id4)

dim(cmap_data_info_control_1)
dim(cmap_data_info_control_5)
dim(cmap_data_info_control_6)

# for cmap_data_info_control_1 - add .cell 
cmap_data_info_control_1$cell_id<-paste(cmap_data_info_control_1$vehicle_scan_id4,".CEL", sep="")

# for cmap_data_info_control_5 - split by "." and keep id
for (x in 1:dim(cmap_data_info_control_5)[1]){
  cmap_data_info_control_5$cell_id[x]<-unlist(strsplit(as.character(cmap_data_info_control_5$cell_id[x]), "[.]"))[1]
}

# for cmap_data_info_control_6 - split by "." and keep id
for (x in 1:dim(cmap_data_info_control_6)[1]){
  cmap_data_info_control_6$cell_id[x]<-unlist(strsplit(as.character(cmap_data_info_control_6$cell_id[x]), "[.]"))[1]
}

head(cmap_data_info_control_5)
head(cmap_data_info_control_6)

# make 5 cpoies
cmap_data_info_control_5_1<-cmap_data_info_control_5
cmap_data_info_control_5_2<-cmap_data_info_control_5
cmap_data_info_control_5_3<-cmap_data_info_control_5
cmap_data_info_control_5_4<-cmap_data_info_control_5
cmap_data_info_control_5_5<-cmap_data_info_control_5

# add vehicle_scan_id4 to cell_id
for (x in 1:dim(cmap_data_info_control_5_1)[1]){
  cmap_data_info_control_5_1$cell_id[x]<-paste(cmap_data_info_control_5_1$cell_id[x], unlist(strsplit(as.character(cmap_data_info_control_5_1$vehicle_scan_id4[x]), "[.]"))[2], sep=".")
}

for (x in 1:dim(cmap_data_info_control_5_2)[1]){
  cmap_data_info_control_5_2$cell_id[x]<-paste(cmap_data_info_control_5_2$cell_id[x], unlist(strsplit(as.character(cmap_data_info_control_5_2$vehicle_scan_id4[x]), "[.]"))[3], sep=".")
}

for (x in 1:dim(cmap_data_info_control_5_3)[1]){
  cmap_data_info_control_5_3$cell_id[x]<-paste(cmap_data_info_control_5_3$cell_id[x], unlist(strsplit(as.character(cmap_data_info_control_5_3$vehicle_scan_id4[x]), "[.]"))[4], sep=".")
}

for (x in 1:dim(cmap_data_info_control_5_4)[1]){
  cmap_data_info_control_5_4$cell_id[x]<-paste(cmap_data_info_control_5_4$cell_id[x], unlist(strsplit(as.character(cmap_data_info_control_5_4$vehicle_scan_id4[x]), "[.]"))[5], sep=".")
}

for (x in 1:dim(cmap_data_info_control_5_5)[1]){
  cmap_data_info_control_5_5$cell_id[x]<-paste(cmap_data_info_control_5_5$cell_id[x], unlist(strsplit(as.character(cmap_data_info_control_5_5$vehicle_scan_id4[x]), "[.]"))[6], sep=".")
}

#add .CEL
cmap_data_info_control_5_1$cell_id<-paste(cmap_data_info_control_5_1$cell_id,".CEL", sep="")
cmap_data_info_control_5_2$cell_id<-paste(cmap_data_info_control_5_2$cell_id,".CEL", sep="")
cmap_data_info_control_5_3$cell_id<-paste(cmap_data_info_control_5_3$cell_id,".CEL", sep="")
cmap_data_info_control_5_4$cell_id<-paste(cmap_data_info_control_5_4$cell_id,".CEL", sep="")
cmap_data_info_control_5_5$cell_id<-paste(cmap_data_info_control_5_5$cell_id,".CEL", sep="")

head(cmap_data_info_control_5_1)
head(cmap_data_info_control_5_2)
head(cmap_data_info_control_5_3)
head(cmap_data_info_control_5_4)
head(cmap_data_info_control_5_5)

# make 6 cpoies
cmap_data_info_control_6_1<-cmap_data_info_control_6
cmap_data_info_control_6_2<-cmap_data_info_control_6
cmap_data_info_control_6_3<-cmap_data_info_control_6
cmap_data_info_control_6_4<-cmap_data_info_control_6
cmap_data_info_control_6_5<-cmap_data_info_control_6
cmap_data_info_control_6_6<-cmap_data_info_control_6

# add vehicle_scan_id4 to cell_id
for (x in 1:dim(cmap_data_info_control_6_1)[1]){
  cmap_data_info_control_6_1$cell_id[x]<-paste(cmap_data_info_control_6_1$cell_id[x], unlist(strsplit(as.character(cmap_data_info_control_6_1$vehicle_scan_id4[x]), "[.]"))[2], sep=".")
}

for (x in 1:dim(cmap_data_info_control_6_2)[1]){
  cmap_data_info_control_6_2$cell_id[x]<-paste(cmap_data_info_control_6_2$cell_id[x], unlist(strsplit(as.character(cmap_data_info_control_6_2$vehicle_scan_id4[x]), "[.]"))[3], sep=".")
}

for (x in 1:dim(cmap_data_info_control_6_3)[1]){
  cmap_data_info_control_6_3$cell_id[x]<-paste(cmap_data_info_control_6_3$cell_id[x], unlist(strsplit(as.character(cmap_data_info_control_6_3$vehicle_scan_id4[x]), "[.]"))[4], sep=".")
}

for (x in 1:dim(cmap_data_info_control_6_4)[1]){
  cmap_data_info_control_6_4$cell_id[x]<-paste(cmap_data_info_control_6_4$cell_id[x], unlist(strsplit(as.character(cmap_data_info_control_6_4$vehicle_scan_id4[x]), "[.]"))[5], sep=".")
}

for (x in 1:dim(cmap_data_info_control_6_5)[1]){
  cmap_data_info_control_6_5$cell_id[x]<-paste(cmap_data_info_control_6_5$cell_id[x], unlist(strsplit(as.character(cmap_data_info_control_6_5$vehicle_scan_id4[x]), "[.]"))[6], sep=".")
}

for (x in 1:dim(cmap_data_info_control_6_6)[1]){
  cmap_data_info_control_6_6$cell_id[x]<-paste(cmap_data_info_control_6_6$cell_id[x], unlist(strsplit(as.character(cmap_data_info_control_6_6$vehicle_scan_id4[x]), "[.]"))[7], sep=".")
}

#add .CEL
cmap_data_info_control_6_1$cell_id<-paste(cmap_data_info_control_6_1$cell_id,".CEL", sep="")
cmap_data_info_control_6_2$cell_id<-paste(cmap_data_info_control_6_2$cell_id,".CEL", sep="")
cmap_data_info_control_6_3$cell_id<-paste(cmap_data_info_control_6_3$cell_id,".CEL", sep="")
cmap_data_info_control_6_4$cell_id<-paste(cmap_data_info_control_6_4$cell_id,".CEL", sep="")
cmap_data_info_control_6_5$cell_id<-paste(cmap_data_info_control_6_5$cell_id,".CEL", sep="")
cmap_data_info_control_6_6$cell_id<-paste(cmap_data_info_control_6_6$cell_id,".CEL", sep="")

head(cmap_data_info_control_6_1)
head(cmap_data_info_control_6_2)
head(cmap_data_info_control_6_3)
head(cmap_data_info_control_6_4)
head(cmap_data_info_control_6_5)
head(cmap_data_info_control_6_6)

# merge all control files back togther
cmap_data_info_control_merged<-rbind(cmap_data_info_control_1,
                                     cmap_data_info_control_5_1,
                                     cmap_data_info_control_5_2,
                                     cmap_data_info_control_5_3,
                                     cmap_data_info_control_5_4,
                                     cmap_data_info_control_5_5,
                                     cmap_data_info_control_6_1,
                                     cmap_data_info_control_6_2,
                                     cmap_data_info_control_6_3,
                                     cmap_data_info_control_6_4,
                                     cmap_data_info_control_6_5,
                                     cmap_data_info_control_6_6)

# subset unique chip and cell id- keep only 2 columns

colnames(cmap_data_info_control_merged)
#cmap_data_info_control_merged_unique<-cmap_data_info_control_merged[c(8,7,16)]
cmap_data_info_control_merged_unique<-cmap_data_info_control_merged[c(2, 7, 8, 11, 16)]
head(cmap_data_info_control_merged_unique)
cmap_data_info_control_merged_unique<-unique(cmap_data_info_control_merged_unique)
head(cmap_data_info_control_merged_unique)
dim(cmap_data_info_control_merged_unique)
table(cmap_data_info_control_merged_unique$array3)
table(cmap_data_info_control_merged_unique$cell2, cmap_data_info_control_merged_unique$array3)

# total number of files = 6100 case + 956 controls

# add control cel names to case cel names

HG_U133A<-c(HG_U133A, as.vector(cmap_data_info_control_merged_unique[cmap_data_info_control_merged_unique$array3=="HG-U133A", 3]))
HT_HG_U133A<-c(HT_HG_U133A, as.vector(cmap_data_info_control_merged_unique[cmap_data_info_control_merged_unique$array3=="HT_HG-U133A", 3]))
HT_HG_U133A_EA<-c(HT_HG_U133A_EA, as.vector(cmap_data_info_control_merged_unique[cmap_data_info_control_merged_unique$array3=="HT_HG-U133A_EA", 3]))

length(HG_U133A)
length(HT_HG_U133A)
length(HT_HG_U133A_EA)

grep("'", HG_U133A)
grep("'", HT_HG_U133A)
grep("'", HT_HG_U133A_EA)

##### SAVE INFO FILE #####

#save info file
setwd(work_dir)
write.table(cmap_data_info, file="cmap_instances_only_clean.txt", sep="\t")
write.table(cmap_data_info_control_merged, file="cmap_controls_only_clean.txt", sep="\t")

#save full control file
setwd(work_dir)
write.table(cmap_data_info_control_merged, file="cmap_controls_only_full.txt", sep="\t")

##### CREATE LIST OF .CEL PER CHIP + TISSUE ######

# split chip by cell line
table(subset(cmap_data_info, array3=="HG-U133A")$cell2)
table(subset(cmap_data_info, array3=="HT_HG-U133A")$cell2)
table(subset(cmap_data_info, array3=="HT_HG-U133A_EA")$cell2)
table(subset(cmap_data_info_control_merged_unique, array3=="HG-U133A")$cell2)
table(subset(cmap_data_info_control_merged_unique, array3=="HT_HG-U133A")$cell2)
table(subset(cmap_data_info_control_merged_unique, array3=="HT_HG-U133A_EA")$cell2)

# extract cell_id from individual cell lines and chip from both drug and control compound
HG_U133A_HL60<-c(subset(cmap_data_info, array3=="HG-U133A" & cell2=="HL60")$cell_id,subset(cmap_data_info_control_merged_unique, array3=="HG-U133A" & cell2=="HL60")$cell_id)
HG_U133A_MCF7<-c(subset(cmap_data_info, array3=="HG-U133A" & cell2=="MCF7")$cell_id,subset(cmap_data_info_control_merged_unique, array3=="HG-U133A" & cell2=="MCF7")$cell_id)
HG_U133A_PC3<-c(subset(cmap_data_info, array3=="HG-U133A" & cell2=="PC3")$cell_id,subset(cmap_data_info_control_merged_unique, array3=="HG-U133A" & cell2=="PC3")$cell_id)
HG_U133A_SKMEL5<-c(subset(cmap_data_info, array3=="HG-U133A" & cell2=="SKMEL5")$cell_id,subset(cmap_data_info_control_merged_unique, array3=="HG-U133A" & cell2=="SKMEL5")$cell_id)
HG_U133A_ssMCF7<-c(subset(cmap_data_info, array3=="HG-U133A" & cell2=="ssMCF7")$cell_id,subset(cmap_data_info_control_merged_unique, array3=="HG-U133A" & cell2=="ssMCF7")$cell_id)

HT_HG_U133A_HL60<-c(subset(cmap_data_info, array3=="HT_HG-U133A" & cell2=="HL60")$cell_id,subset(cmap_data_info_control_merged_unique, array3=="HT_HG-U133A" & cell2=="HL60")$cell_id)
HT_HG_U133A_MCF7<-c(subset(cmap_data_info, array3=="HT_HG-U133A" & cell2=="MCF7")$cell_id,subset(cmap_data_info_control_merged_unique, array3=="HT_HG-U133A" & cell2=="MCF7")$cell_id)
HT_HG_U133A_PC3<-c(subset(cmap_data_info, array3=="HT_HG-U133A" & cell2=="PC3")$cell_id,subset(cmap_data_info_control_merged_unique, array3=="HT_HG-U133A" & cell2=="PC3")$cell_id)
HT_HG_U133A_SKMEL5<-c(subset(cmap_data_info, array3=="HT_HG-U133A" & cell2=="SKMEL5")$cell_id,subset(cmap_data_info_control_merged_unique, array3=="HT_HG-U133A" & cell2=="SKMEL5")$cell_id)
HT_HG_U133A_ssMCF7<-c(subset(cmap_data_info, array3=="HT_HG-U133A" & cell2=="ssMCF7")$cell_id,subset(cmap_data_info_control_merged_unique, array3=="HT_HG-U133A" & cell2=="ssMCF7")$cell_id)

HT_HG_U133A_EA_HL60<-c(subset(cmap_data_info, array3=="HT_HG-U133A_EA" & cell2=="HL60")$cell_id,subset(cmap_data_info_control_merged_unique, array3=="HT_HG-U133A_EA" & cell2=="HL60")$cell_id)
HT_HG_U133A_EA_MCF7<-c(subset(cmap_data_info, array3=="HT_HG-U133A_EA" & cell2=="MCF7")$cell_id,subset(cmap_data_info_control_merged_unique, array3=="HT_HG-U133A_EA" & cell2=="MCF7")$cell_id)
HT_HG_U133A_EA_PC3<-c(subset(cmap_data_info, array3=="HT_HG-U133A_EA" & cell2=="PC3")$cell_id,subset(cmap_data_info_control_merged_unique, array3=="HT_HG-U133A_EA" & cell2=="PC3")$cell_id)
HT_HG_U133A_EA_SKMEL5<-c(subset(cmap_data_info, array3=="HT_HG-U133A_EA" & cell2=="SKMEL5")$cell_id,subset(cmap_data_info_control_merged_unique, array3=="HT_HG-U133A_EA" & cell2=="SKMEL5")$cell_id)
HT_HG_U133A_EA_ssMCF7<-c(subset(cmap_data_info, array3=="HT_HG-U133A_EA" & cell2=="ssMCF7")$cell_id,subset(cmap_data_info_control_merged_unique, array3=="HT_HG-U133A_EA" & cell2=="ssMCF7")$cell_id)

length(HG_U133A_HL60)
length(HG_U133A_MCF7)
length(HG_U133A_PC3)
length(HG_U133A_SKMEL5)
length(HG_U133A_ssMCF7)

length(HT_HG_U133A_HL60)
length(HT_HG_U133A_MCF7)
length(HT_HG_U133A_PC3)
length(HT_HG_U133A_SKMEL5) # none
length(HT_HG_U133A_ssMCF7) # none

length(HT_HG_U133A_EA_HL60) #none
length(HT_HG_U133A_EA_MCF7)
length(HT_HG_U133A_EA_PC3) # none
length(HT_HG_U133A_EA_SKMEL5) # none
length(HT_HG_U133A_EA_ssMCF7) #none

# check it adds upto 7056 - total number of .CEL files

length(HG_U133A_HL60)+
length(HG_U133A_MCF7)+
length(HG_U133A_PC3)+
length(HG_U133A_SKMEL5)+
length(HG_U133A_ssMCF7)+
length(HT_HG_U133A_HL60)+
length(HT_HG_U133A_MCF7)+
length(HT_HG_U133A_PC3)+
length(HT_HG_U133A_EA_MCF7)

##### READ IN RAW CMAP DATA #####

#read in data by chip+tissue
setwd(unzipped_dir)

# only reading MCF7

HG_U133A_MCF7_raw<-ReadAffy(filenames = HG_U133A_MCF7)
HT_HG_U133A_MCF7_raw<-ReadAffy(filenames = HT_HG_U133A_MCF7)
HT_HG_U133A_EA_MCF7_raw<-ReadAffy(filenames = HT_HG_U133A_EA_MCF7)

HG_U133A_MCF7_raw
HT_HG_U133A_MCF7_raw
HT_HG_U133A_EA_MCF7_raw

# change annotation


# save raw data as object
setwd(Raw_data_R_object_dir)

save(HG_U133A_MCF7_raw, file="HG_U133A_MCF7.Rdata")
save(HT_HG_U133A_MCF7_raw, file="HT_HG_U133A_MCF7.Rdata")
save(HT_HG_U133A_EA_MCF7_raw, file="HT_HG_U133A_EA_MCF7.Rdata")

##### PRE-PROCESS - BACKGROUND CORRECT + NORMALISE #####

setwd(background_corrected_dir)

HG_U133A_MCF7_bc<-gcrma(HG_U133A_MCF7_raw)
HT_HG_U133A_MCF7_bc<-gcrma(HT_HG_U133A_MCF7_raw)
HT_HG_U133A_EA_MCF7_bc<-gcrma(HT_HG_U133A_EA_MCF7_raw)

#save bc + normalise - does quantile normalisation

save(HG_U133A_MCF7_bc, file="HG_U133A_MCF7_background_corrected.Rdata")
save(HT_HG_U133A_MCF7_bc, file="HT_HG_U133A_MCF7_background_corrected.Rdata")
save(HT_HG_U133A_EA_MCF7_bc, file="HT_HG_U133A_EA_MCF7_background_corrected.Rdata")

# ##### PRE-PROCESS - NORMALISE #####
# 
# # partially done on rosalind due to lack of computing power (HT_HG_U133A chip require ~200 gb ram).
# 
# HG_U133A_MCF7_normalised<-rsn(log2(exprs(HG_U133A_MCF7_bc)))
# HHT_HG_U133A_MCF7_normalised<-rsn(log2(exprs(HT_HG_U133A_MCF7_bc)))
# HT_HG_U133A_EA_MCF7_normalised<-rsn(log2(exprs(HT_HG_U133A_EA_MCF7_bc)))

##### SAVE NORMALISED TABLES #####

# setwd(normalised_dir)
# 
# save(HG_U133A_HL60_normalised, file="HG_U133A_HL60_normalised.Rdata")
# save(HG_U133A_MCF7_normalised, file="HG_U133A_MCF7_normalised.Rdata")
# save(HG_U133A_PC3_normalised, file="HG_U133A_PC3_normalised.Rdata")
# save(HG_U133A_SKMEL5_normalised, file="HG_U133A_SKMEL5_normalised.Rdata")
# save(HG_U133A_ssMCF7_normalised, file="HG_U133A_ssMCF7_normalised.Rdata")
# 
# save(HT_HG_U133A_HL60_normalised, file="HT_HG_U133A_HL60_normalised.Rdata")
# save(HT_HG_U133A_MCF7_normalised, file="HT_HG_U133A_MCF7_normalised.Rdata")
# save(HT_HG_U133A_PC3_normalised, file="HT_HG_U133A_PC3_normalised.Rdata")
# 
# save(HT_HG_U133A_EA_MCF7_normalised, file="HT_HG_U133A_EA_MCF7_normalised.Rdata")
# 
# 
# 
