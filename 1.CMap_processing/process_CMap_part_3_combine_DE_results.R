##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                          CMap build 2.0 PROCESSING - PART 3 -                                                          #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################

## AUTHOR: HAMEL PATEL
## EMAIL: Hamel.patel@kcl.ac.uk
## DATE: 17/10/2016

##### DESCRIPTION OF ANALYSIS ####
## MCF7 cell line across 3 chips has been processed and DE analysis performed.
## This script will combine information across 3 chips to create a lookup reference database
## WORKING WITH FULL PROBE DATASETS - 12443
#####

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### SET DIRECTORIES ####

work_dir="/media/hamel/Workspace/Dropbox/Projects/CMap2/1.Data/5.Explore_DE_Results/"

setwd(work_dir)

HG_U133A_MCF7_dir<-"/media/hamel/Workspace/Dropbox/Projects/CMap2/1.Data/2.HG_U133A_MCF7_processing/DE_results"
HT_HG_U133A_MCF7_dir<-"/media/hamel/Workspace/Dropbox/Projects/CMap2/1.Data/3.HT_HG_U133A_MCF7_processing/DE_results/"

# create dir to plot overlap of same drug in different chip
dir.create(paste(work_dir, "Volcano_plots"))

##### LIBRARY #####
library(calibrate)

##### LOAD DATA ####

setwd(HG_U133A_MCF7_dir)
load("full_DE_results.Rdata")
HG_U133A_MCF7<-full_DE_results

setwd(HT_HG_U133A_MCF7_dir)
load("full_DE_results.Rdata")
HT_HG_U133A_MCF7<-full_DE_results

rm(full_DE_results)

##### CHECK DATA #####

dim(HT_HG_U133A_MCF7[[1]])
dim(HG_U133A_MCF7[[1]])

names(HT_HG_U133A_MCF7)
names(HG_U133A_MCF7)

##### COMPARE DUPLICATED DRUGS #####

duplicate_drugs<-subset(names(HT_HG_U133A_MCF7), names(HT_HG_U133A_MCF7) %in% names(HG_U133A_MCF7))

duplicate_drugs

head(HT_HG_U133A_MCF7[[grep("alpha-estradiol_1e-08", names(HT_HG_U133A_MCF7))]])
head(HG_U133A_MCF7[[grep("alpha-estradiol_1e-08", names(HG_U133A_MCF7))]])

##### VOLCANO PLOTS ######

# function
plot_volcano<-function(dataset, drug, adj.pval_cutt_off, logFc_cutoff, label ){
  # grep drug
  dataset_subset<-as.data.frame(dataset[(grep(drug, names(dataset)))])
  #add rownames to column
  dataset_subset$Entrez<-rownames(dataset_subset)
  #add colnames
  colnames(dataset_subset)<-c("logFC", "CI.L", "CI.R", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
  # plot
  with(dataset_subset, plot(logFC, -log10(P.Value), pch=20, main=label, xlim=c(-6,6)))
  # Add colored points: red if adj.P.Val<0adj.pval_cutt_off, orange of loglogFc_cutoffFC>1, green if both)
  with(subset(dataset_subset, adj.P.Val<adj.pval_cutt_off ), points(logFC, -log10(P.Value), pch=20, col="red"))
  with(subset(dataset_subset, abs(logFC)>logFc_cutoff), points(logFC, -log10(P.Value), pch=20, col="orange"))
  with(subset(dataset_subset, adj.P.Val<adj.pval_cutt_off & abs(logFC)>logFc_cutoff), points(logFC, -log10(P.Value), pch=20, col="green"))
  # Label points with the textxy function from the calibrate plot
  # add entrez id to plot
  with(subset(dataset_subset, adj.P.Val<adj.pval_cutt_off & abs(logFC)>3), textxy(logFC, -log10(P.Value), labs=Entrez, cex=.8))
}

#apply
plot_volcano(HT_HG_U133A_MCF7, duplicate_drugs[1],  0.05, 2, duplicate_drugs[1])
plot_volcano(HG_U133A_MCF7, duplicate_drugs[1],  0.05, 2, duplicate_drugs[1])

plot_volcano(HT_HG_U133A_MCF7, duplicate_drugs[2],  0.05, 2, duplicate_drugs[2])
plot_volcano(HG_U133A_MCF7, duplicate_drugs[2],  0.05, 2, duplicate_drugs[2])

plot_volcano(HT_HG_U133A_MCF7, duplicate_drugs[3],  0.05, 2, duplicate_drugs[3])
plot_volcano(HG_U133A_MCF7, duplicate_drugs[3],  0.05, 2, duplicate_drugs[3])

plot_volcano(HT_HG_U133A_MCF7, duplicate_drugs[4],  0.05, 2, duplicate_drugs[4])
plot_volcano(HG_U133A_MCF7, duplicate_drugs[4],  0.05, 2, duplicate_drugs[4])

plot_volcano(HT_HG_U133A_MCF7, duplicate_drugs[5],  0.05, 2, duplicate_drugs[5])
plot_volcano(HG_U133A_MCF7, duplicate_drugs[5],  0.05, 2, duplicate_drugs[5])

plot_volcano(HT_HG_U133A_MCF7, duplicate_drugs[6],  0.05, 2, duplicate_drugs[6])
plot_volcano(HG_U133A_MCF7, duplicate_drugs[6],  0.05, 2, duplicate_drugs[6])

plot_volcano(HT_HG_U133A_MCF7, duplicate_drugs[7],  0.05, 2, duplicate_drugs[7])
plot_volcano(HG_U133A_MCF7, duplicate_drugs[7],  0.05, 2, duplicate_drugs[7])

plot_volcano(HT_HG_U133A_MCF7, duplicate_drugs[8],  0.05, 2, duplicate_drugs[8])
plot_volcano(HG_U133A_MCF7, duplicate_drugs[8],  0.05, 2, duplicate_drugs[8])

plot_volcano(HT_HG_U133A_MCF7, duplicate_drugs[9],  0.05, 2, duplicate_drugs[9])
plot_volcano(HG_U133A_MCF7, duplicate_drugs[9],  0.05, 2, duplicate_drugs[9])



setwd(work_dir)
pdf("Trichostatin_A_across_3_chips_volcano_plot.pdf")
plot_volcano(test2, 0.05, 2, "Trichostatin A_1e-07 HT_HG_U133A_MCF7")
plot_volcano(test3, 0.05, 2, "Trichostatin A_1e-07 HG_U133A_MCF7")
dev.off()

ls()

##### CURRENT AD DRUGS #######

names(HT_HG_U133A_MCF7)[grep("galantamine", names(HT_HG_U133A_MCF7))]
grep("galantamine", names(HG_U133A_MCF7))

names(HT_HG_U133A_MCF7)[grep("memantine", names(HT_HG_U133A_MCF7))]
grep("memantine", names(HG_U133A_MCF7))

test4<-HT_HG_U133A_MCF7[[grep("galantamine", names(HT_HG_U133A_MCF7))]]
test5<-HT_HG_U133A_MCF7[[grep("memantine", names(HT_HG_U133A_MCF7))]]

plot_volcano(test4, 0.05, 2, "galantamine")
plot_volcano(test5, 0.05, 2, "memantine")

##### COUNT DEG #######

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

head(number_of_DEG_in_full_data, 10)
tail(number_of_DEG_in_full_data)

subset(number_of_DEG_in_full_data, Drug=="memantine_1.86e-05")
subset(number_of_DEG_in_full_data, Drug=="galantamine_1.08e-05")

# write

setwd(work_dir)

write.table(number_of_DEG_in_full_data, "HT_HG_U133A_MCF7_drug_DE_count.csv", row.names = F, sep=",")
