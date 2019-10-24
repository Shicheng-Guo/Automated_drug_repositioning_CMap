##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                                   RUN SPIA ON CMAP DATA                                                                #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################

## AUTHOR: HAMEL PATEL
## EMAIL: Hamel.patel@kcl.ac.uk
## DATE: 16/05/2017 

##### DESCRIPTION OF ANALYSIS ####
## RUN SPIA ON CMAP AND CREATE QUERY DATABASE
#####

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### LIBRARY #####

library(SPIA)
library(foreach)
library(doParallel)

##### SET DIRECTORIES ####

work_dir="/media/hamel/Workspace/Dropbox/Projects/CMap2/3.Method_development/5.with_pathway_knowledge/"
cmap_data_dir="/media/hamel/Workspace/Dropbox/Projects/CMap2/3.Method_development/1.single_drug_bi_directional_enrichment/"
#HT_HG_U133A_MCF7_dir<-"/media/hamel/Workspace/Dropbox/Projects/Brain_expression/1.Data/1.Re-process_with_new_pipeline/CMAP/HT_HG_U133A_MCF7_processing/DE_results"
#results_dir="/media/hamel/Workspace/Dropbox/Projects/Brain_expression/6.Query_CMap/7.6.Method_validation/CMap_query_results/"

##### LOAD CMAP DATA #####

# load cmap query database
setwd(cmap_data_dir)
load("cmap_database_for_bi_directional_enrichment_search.Rdata")
cmap_database<-subset_CMap
rm(subset_CMap)

head(cmap_database[[1]])
names(cmap_database)[1]

length(cmap_database)

#create empty list to store SPIA results
cmap_SPIA<-list()

# run spia on all compounds

##### CREATE CLUSTER #####
#ceate cluster for multi-threading
#Calculate the number of cores
no_cores <- detectCores()
registerDoParallel(no_cores)

#for windows
#cores<- detectCores(logical=T)
#cl<-makeCluster(cores)
#cl<-makeCluster(22) # each core requires 8gb ram?
#registerDoParallel(cl)

cmap_SPIA<-foreach(x=1:length(cmap_database), .verbose=T) %dopar% {
  #subset drug DE to sigDE only 
  subset_drug_DE<-subset(cmap_database[[x]], adj.P.Val<=0.05)
  #create spia input
  subset_drug_DE_SPIA_format<-subset_drug_DE$logFC
  names(subset_drug_DE_SPIA_format)<-as.vector(rownames(subset_drug_DE))
  #run SPIA
  cmap_SPIA<-spia(de=subset_drug_DE_SPIA_format, all=rownames(cmap_database[[x]]), organism="hsa", plots=F)
  return(cmap_SPIA)
}

names(cmap_SPIA)<-names(cmap_database)

gc(verbose=T)
#close cluster
stopImplicitCluster()

head(cmap_SPIA[[1]])
head(cmap_SPIA[[2]])
head(cmap_SPIA[[3]])

grep("phenoxybenzamine", names(cmap_SPIA))
head(cmap_SPIA[[792]])

cmap_SPIA[[792]]$Name

# subset pathways to sig only

#create empty list
cmap_SPIA_sig<-list()

for (x in 1:length(cmap_SPIA)) {
  # subset to sig pathways
  cmap_SPIA_sig[[x]]<-subset(cmap_SPIA[[x]], pG<=0.05)
  #add names
  names(cmap_SPIA_sig)[x]<-names(cmap_SPIA)[x]
}

head(cmap_SPIA_sig[[1]])
names(cmap_SPIA_sig[1])
length(cmap_SPIA_sig)

# subset to drugs with a pathway perturbed

cmap_SPIA_sig<-cmap_SPIA_sig[lapply(cmap_SPIA_sig,nrow)>0]
head(cmap_SPIA_sig[[1]])
names(cmap_SPIA_sig[1])
length(cmap_SPIA_sig)


##### SAVE #####

setwd(work_dir)

save(cmap_SPIA_sig, file="CMap_SPIA_subset_to_sig.Rdata")

save.image("CREATE_CMAP_PATHWAY_DATABASE.Rdata")

save(cmap_SPIA, file="cmap_pathway_database.Rdata")
