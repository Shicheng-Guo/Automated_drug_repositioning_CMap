##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                          CMap build 2.0 PROCESSING - PART 2 -                                                          #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################

# MICROARRAY PLATFORM - Affymetrix
# EXPRESSION CHIP - HG-U133A 
# CELL LINE - MCF7
# NUMBER OF SAMPLES - 171 + 47 (controls)

# NOTES
# - full dataset only - no filtering by non-expressed genes

## AUTHOR: HAMEL PATEL
## EMAIL: Hamel.patel@kcl.ac.uk
## DATE: 10/10/2016

##### DESCRIPTION OF ANALYSIS ####
## CMAP build 2 downloaded from from http://portals.broadinstitute.org/cmap/ on 4/08/2016.
## Data has been background corrected (MAS5), log2, and normalised (RSN) by expression chip and tissue.
#####

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

# name of data to be processed
data_name="HG_U133A_MCF7"

##### SET DIRECTORIES ####

#set existing directories and create new directories
data_dir="/media/hamel/Workspace/Dropbox/Projects/CMap2/1.Data/1.pre-process_data/3.Normalised_data/"

work_dir="/media/hamel/Workspace/Dropbox/Projects/CMap2/1.Data/"

pheno_info_dir<-"/media/hamel/Workspace/Dropbox/Projects/CMap2/1.Data/1.pre-process_data/2.CMap_pheno_data/"

setwd(work_dir)

# create dir
dir.create(paste(work_dir, paste(data, "_processing", sep=""), sep="/"))
processing_dir=(paste(work_dir, paste("2.", data_name, "_processing", sep=""), sep="/"))

setwd(processing_dir)

# create dir for boxplots
dir.create(paste(processing_dir,"boxplots_density_plots", sep="/"))
boxplots_density_plots_dir=paste(processing_dir,"boxplots_density_plots", sep="/")

# create dir for PCA plots
dir.create(paste(processing_dir,"PCA_plots", sep="/"))
PCA_plots_dir=paste(processing_dir,"PCA_plots", sep="/")

# create dir sample network plots
dir.create(paste(processing_dir,"sample_network", sep="/"))
sample_network_dir=paste(processing_dir,"sample_network", sep="/")

# create dir for clean data
dir.create(paste(processing_dir,"Clean_data", sep="/"))
clean_dir=paste(processing_dir,"Clean_data", sep="/")

# create directory for differential expression
dir.create(paste(processing_dir,"DE_results", sep="/"))
DE_dir=paste(processing_dir,"DE_results", sep="/")

##### LOAD LIBRARIES #####

library(WGCNA)
library(lumi)
library(sva)
library(massiR)
library(reshape)
library(ggplot2)
library(limma)
library(stringr)
library(foreach)
library(doParallel)
library(hgu133a.db)
#library(hthgu133a.db)

##### LOAD NORMALISED DATA #####

setwd(data_dir)

load(paste(data_name, "_background_corrected.Rdata", sep=""))

# assign expression object 

exprs_data<-eval(parse(text=paste(data_name, "_bc", sep="")))
exprs_data<-exprs(exprs_data)

rm(list= paste(data_name, "_bc", sep=""))

tail(exprs_data)[,1:5]
dim(exprs_data)

##### SET NEGATIVE VALUE TO ZERO #####

# set negative values to zero

exprs_data[exprs_data<0]<-0

##### READ PHENO INFO #####

#read in pheno information cmap data
setwd(pheno_info_dir)

pheno_info_drug<-read.table("cmap_instances_only_clean.txt", head=T, sep="\t", as.is=T)
pheno_info_control<-read.table("cmap_controls_only_clean.txt", head=T, sep="\t", as.is=T)
pheno_info_control_FULL<-read.table("cmap_controls_only_full.txt", head=T, sep="\t", as.is=T)

head(pheno_info_drug)
head(pheno_info_control)

# subset control to unique

colnames(pheno_info_control)

pheno_info_control_unique<-pheno_info_control[c(2, 3, 5, 7, 8, 11, 16)]
head(pheno_info_control_unique)
pheno_info_control_unique<-subset(pheno_info_control_unique, !duplicated(cell_id))
head(pheno_info_control_unique)
dim(pheno_info_control_unique)

table(pheno_info_control_unique$array3)
table(pheno_info_control_unique$cell2, pheno_info_control_unique$array3)

#keep only pheno info for samples available

subset_pheno_drug<-subset(pheno_info_drug, pheno_info_drug$cell_id %in% colnames(exprs_data))
subset_pheno_control<-subset(pheno_info_control_unique, pheno_info_control_unique$cell_id %in% colnames(exprs_data))

#create unique column with drug_concentration
subset_pheno_drug$drug_conc<-paste(subset_pheno_drug$cmap_name, subset_pheno_drug$concentration..M., sep="_")
subset_pheno_control$drug_conc<-"Control"

#create group
subset_pheno_drug$group<-"Drug"
subset_pheno_control$group<-"Control"

#sample numbers same?
dim(subset_pheno_drug)[1]+dim(subset_pheno_control)[1]==dim(exprs_data)[2]

# merge subsets together - keep only batch_id, cell2, array3, scanner, cell_id and create column for drug or control
subset_pheno_drug<-subset_pheno_drug[c(2, 3, 5, 7, 8, 11, 16, 17, 18)]
all(colnames(subset_pheno_drug)==colnames(subset_pheno_control))==T
pheno<-rbind(subset_pheno_drug, subset_pheno_control)
head(pheno)

dim(pheno)[1]==dim(exprs_data)[2]

##### CHECK FOR DUPLICATE SAMPLES IDs #####

anyDuplicated(colnames(exprs_data))

##### INITIAL PLOTS #####

setwd(boxplots_density_plots_dir)

jpeg(file=paste(data_name, "_normalised_boxplot.jpeg", sep=""))
boxplot(exprs_data)
dev.off()

jpeg(file=paste(data_name, "_normalised_density_plot.jpeg", sep=""))
plotDensity(exprs_data[,1:100], logMode=F, addLegend=F)
dev.off()

##### PCA PLOT 1 #####

#create function to pCA plot - will be done before and after each QC step

pca_data<-function(data, pca_pheno, legend_position){
  pca<-prcomp(data)
  #summary_pca<-summary(pca)
  # order of samples in expression data
  sample_order<-colnames(data)
  # subset pca_pheno to match expression data
  pca_pheno<-subset(pheno, cell_id %in% colnames(data))
  # match order
  ordered_pca_pheno<-pca_pheno[match(sample_order, pca_pheno$cell_id),]
  batch_id_pca_colour<-labels2colors(as.character(ordered_pca_pheno$batch_id))
  scanner_pca_colour<-labels2colors(as.character(ordered_pca_pheno$scanner))
  group_pca_colour<-labels2colors(as.character(ordered_pca_pheno$group))
  drug_name_pca_colour<-labels2colors(as.numeric(as.numeric(as.factor(ordered_pca_pheno$drug_conc))))
  #remove . from colour names
  drug_name_pca_colour<-gsub("\\..*","",drug_name_pca_colour)
  drug_conc_pca_colour<-labels2colors(as.character(ordered_pca_pheno$concentration..M.))
  # pca plot - Batch ID
  plot(pca$rotation[,1:2], main=" PCA plot coloured by Batch ID",col="black", pch=21,bg=batch_id_pca_colour)
  legend(legend_position, legend=unique(ordered_pca_pheno$batch_id), fill=unique(batch_id_pca_colour), title="Batch ID")
  # pca plot - Scanner
  plot(pca$rotation[,1:2], main=" PCA plot coloured by Scanner",col="black", pch=21,bg=scanner_pca_colour)
  legend(legend_position, legend=unique(ordered_pca_pheno$scanner), fill=unique(scanner_pca_colour), title="Scanner")
  # pca plot - group
  plot(pca$rotation[,1:2], main=" PCA plot coloured by Group",col="black", pch=21,bg=group_pca_colour)
  legend(legend_position, legend=unique(ordered_pca_pheno$group), fill=unique(group_pca_colour), title="Group")
  # pca plot - drug name
  plot(pca$rotation[,1:2], main=" PCA plot coloured by Drug name",col="black", pch=21,bg=drug_name_pca_colour)
  legend(legend_position, legend=unique(ordered_pca_pheno$cmap_name), fill=unique(drug_name_pca_colour), title="Drug Name")
  # pca plot - concentration
  plot(pca$rotation[,1:2], main=" PCA plot coloured by Drug Concentration",col="black", pch=21,bg=drug_conc_pca_colour)
  legend(legend_position, legend=unique(ordered_pca_pheno$concentration..M.), fill=unique(drug_conc_pca_colour), title="Drug Concentration")
}

#apply function
pca_data(exprs_data, pheno, 'bottomright')

#plot to pdf

setwd(PCA_plots_dir)

pdf("1.PCA_plot_before_QC.pdf")
pca_data(exprs_data, pheno, 'bottomright')
dev.off()

##### COMBAT  - ERROR CAUSED BY LOW VARINACE IN DATA - using limma #####

# move cell id to rownames
rownames(pheno)<-pheno$cell_id

#check same order # should be T
all(rownames(pheno)==colnames(exprs_data))==T

# copy pheno for pca + convert to factor
pheno_pca<-as.data.frame(unclass(pheno))
str(pheno_pca)

#extract batch info
batch<-pheno_pca$batch_id

#create model matrix
modcombat = model.matrix(~1, data=pheno_pca)
head(modcombat)

#apply combat
#exprs_data_combat<-ComBat(dat=exprs_data, batch=batch, mod=modcombat, par.prior = T, prior.plots = T)

exprs_data_combat<- removeBatchEffect(exprs_data, batch=batch, design=modcombat)

##### PCA PLOT 2 #####

#create function to pCA plot - will be done before and after each QC step

#apply function
pca_data(exprs_data_combat, pheno, 'bottomright')

#plot to pdf

setwd(PCA_plots_dir)

pdf("2.PCA_plot_after_combat.pdf")
pca_data(exprs_data_combat, pheno, 'bottomright')
dev.off()

##### SUBSET TO SAMPLES WHERE MORE THAN 1 REPEAT #####

dim(pheno)
table(pheno$group)

#subset pheno to where drug repeated at least once
repeated_drugs<-subset(as.data.frame(table(pheno$drug_conc)), Freq>1)
dim(repeated_drugs)
subset_pheno<-subset(pheno,  drug_conc%in%repeated_drugs$Var1)

dim(subset_pheno)
head(subset_pheno)

#subset expression data
exprs_data_combat_subset<-exprs_data_combat[,subset_pheno$cell_id]
dim(exprs_data_combat_subset)
dim(exprs_data_combat)

#check pheno + expression same dimensions
dim(exprs_data_combat_subset)[2]==dim(subset_pheno)[1]

##### PCA PLOT 3 #####

#create function to pCA plot - will be done before and after each QC step

#apply function
pca_data(exprs_data_combat_subset, subset_pheno, 'bottomright')

#plot to pdf

setwd(PCA_plots_dir)

pdf("3.PCA_plot_duplicates_only.pdf")
pca_data(exprs_data_combat_subset, subset_pheno, 'bottomright')
dev.off()

##### SVA #####

head(subset_pheno)
head(t(exprs_data_combat_subset))[,1:5]

# transpose exprs data
exprs_data_combat_subset_t<-as.data.frame(t(exprs_data_combat_subset))
dim(exprs_data_combat_subset_t)

# move cell id in pheno to row name
rownames(subset_pheno)<-subset_pheno$cell_id
head(subset_pheno)
dim(subset_pheno)

# add GROUP in exprs table
exprs_data_with_group<-merge(subset_pheno[9], exprs_data_combat_subset_t, by="row.names")
rownames(exprs_data_with_group)<-exprs_data_with_group$Row.names
exprs_data_with_group$Row.names<-NULL

head(exprs_data_with_group)[1:5]

# create sva function

check_SV_in_data<-function(dataset){
  # create sva compatable matrix - sample in columns, probes in rows - pheno info seperate - sort by group 1st to keep AD top
  sorted_by_group<-dataset[order(dataset$group),]
  # separate expresion and pheno
  dataset_pheno<-sorted_by_group[1]
  dataset_exprs<-t(sorted_by_group[2:dim(sorted_by_group)[2]])
  #full model matrix for group
  mod = model.matrix(~group, data=dataset_pheno)
  # check number of SV in data
  print(num.sv(dataset_exprs, mod, method="leek"))
}

# check sv and store

number_of_SV_check<-check_SV_in_data(exprs_data_with_group)

#make copy - as original will be over written SV loop
number_of_SV<-number_of_SV_check
# create function to sva adjust - if stored sva > 0 then adjust

# create function to adjust for SVA and adjust if needed - exclude gender as creating massive mbatch effect

if (number_of_SV>0){
  adjust_for_sva<-function(dataset){
    # create sva compatable matrix - sample in columns, probes in rows - pheno info seperate - sort by group 1st to keep AD top
    sorted_by_group<-dataset[order(dataset$group),]
    # separate expresion and pheno
    dataset_sva_pheno<-sorted_by_group[1]
    dataset_sva_exprs<-t(sorted_by_group[2:dim(sorted_by_group)[2]])
    #full model matrix for group
    mod = model.matrix(~group, data=dataset_sva_pheno)
    mod0 = model.matrix(~1, data=dataset_sva_pheno)
    # number of SV
    num.sv(dataset_sva_exprs, mod, method="leek")
    n.sv=num.sv(dataset_sva_exprs, mod, method="leek")
    # exit if n.sv=0
    if(n.sv==0){stop("No Significant Variable found, exiting....")}
    # apply sva - removed n.sv
    svobj = sva(dataset_sva_exprs, mod, mod0, n.sv=n.sv, method="two-step")
    # adjust for sva
    X = cbind(mod, svobj$sv)
    Hat = solve(t(X) %*% X) %*% t(X)
    beta = (Hat %*% t(dataset_sva_exprs))
    P = ncol(mod)
    clean_data<-dataset_sva_exprs - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),])
    # merge clean data with pheno
    clean_data_with_pheno<-merge(dataset_sva_pheno, as.data.frame(t(clean_data)), by="row.names")
    rownames(clean_data_with_pheno)<-clean_data_with_pheno$Row.names
    clean_data_with_pheno$Row.names<-NULL
    # check SVA on adjusted data
    cat("\n")
    cat("number of surrogate variables after adjustment:")
    cat("\n")
    print(num.sv(clean_data, mod, method="leek"))
    # return clean data with pheno
    return(clean_data_with_pheno)
  }
  while (number_of_SV>0){
    iteration=1
    exprs_data_sva<-adjust_for_sva(exprs_data_with_group)
    exprs_data_with_group<-exprs_data_sva
    print(paste("Iteration: ", iteration, sep=""))
    number_of_SV<-check_SV_in_data(exprs_data_with_group)
    iteration=iteration+1}} else
      {
        exprs_data_sva<-exprs_data_with_group
        check_SV_in_data(exprs_data_sva)
        }

#double check SV
check_SV_in_data(exprs_data_sva)

##### PCA PLOT 4 #####

head(exprs_data_sva)[1:5]
dim(exprs_data_sva)

pca_data(t(exprs_data_sva[2:dim(exprs_data_sva)[2]]), pheno ,'topright')

setwd(PCA_plots_dir)

pdf("4.PCA_plot_after_sva.pdf")
pca_data(t(exprs_data_sva[2:dim(exprs_data_sva)[2]]), pheno, 'topright')
dev.off()

##### SAMPLE NETWORK ANALYSIS #####

# filtering control samples only. testing on drug samples to see netwrok.

# sample plot function - taken from steve expression pipeline # adjusted - removed dendogram as error in C stack due to large number of samples.

sampleNetwork_plot <- function(dataset) {
  datExprs<-t(dataset[2:dim(dataset)[2]])
  diagnosis<-dataset[1]
  gp_col <- "group"
  cat(" setting up data for qc plots","\r","\n")
  ## expression matrix and IAC
  cat(" expression matrix and IAC","\r","\n")
  IAC <- cor(datExprs)
  IAC_d <- 1-IAC
  samle_names <- colnames(datExprs)
  IAC=cor(datExprs, method="p",use="p")
  diag(IAC)=0
  A.IAC=((1+IAC)/2)^2 ## ADJACENCY MATRIX
  cat(" fundamentalNetworkConcepts","\r","\n")
  FNC=fundamentalNetworkConcepts(A.IAC) ## WGCNA
  K2=FNC$ScaledConnectivity
  Z.K=(K2-mean(K2))/sd(K2)
  Z.C=(FNC$ClusterCoef-mean(FNC$ClusterCoef))/sd(FNC$ClusterCoef)
  rho <- signif(cor.test(Z.K,Z.C,method="s")$estimate,2)
  rho_pvalue <- signif(cor.test(Z.K,Z.C,method="s")$p.value,2)
  # set colours
  cat(" colorvec [",paste(gp_col),"]","\r","\n")
  if(gp_col=="chip") { colorvec <- labels2colors(as.character(pData(eset)$Sentrix.Barcode)) }
  if(gp_col=="group") { colorvec <- labels2colors(diagnosis[1]) }
  mean_IAC <- mean(IAC[upper.tri(IAC)])
  ## samplenetwork
  local(
    {colLab <<- function(n,treeorder) {
      if(is.leaf(n)) {
        a <- attributes(n)
        i <<- i+1
        attr(n, "nodePar") <- c(a$nodePar, list(lab.col = colorvec[treeorder][i], lab.font = i%%3))
      }
      n
    }
    i <- 0
    })
  cat(" begin SampleNetwork plots","\r","\n")
  group_colours<-unique(cbind(colorvec, diagnosis))
  ## Cluster for pics
  cluster1 <- hclust(as.dist(1-A.IAC),method="average")
  cluster1order <- cluster1$order
  cluster2 <- as.dendrogram(cluster1,hang=0.1)
  ######cluster3 <- dendrapply(cluster2,colLab,cluster1order)
  ## PLOTS
  ## cluster IAC
  par(mfrow=c(2,2))
  par(mar=c(5,6,4,2))
  # plot(cluster3,nodePar=list(lab.cex=1,pch=NA),
  #      main=paste("Mean ISA = ",signif(mean(A.IAC[upper.tri(A.IAC)]),3),sep=""),
  #      xlab="",ylab="1 - ISA",sub="",cex.main=1.8,cex.lab=1.4)
  # mtext(paste("distance: 1 - ISA ",sep=""),cex=0.8,line=0.2)
  ## Connectivity
  par(mar=c(5,5,4,2))
  plot(Z.K,main="Connectivity", ylab="Z.K",xaxt="n",xlab="Sample",type="n",cex.main=1.8,cex.lab=1.4)
  text(Z.K,labels=samle_names,cex=0.8,col=colorvec)
  abline(h=-2)
  abline(h=-3)
  par(mar=c(5,5,4,2))
  plot(Z.K,Z.C,main="Connectivity vs ClusterCoef",xlab="Z.K",ylab="Z.C",col=colorvec ,cex.main=1.8,cex.lab=1.4)
  abline(lm(Z.C~Z.K),col="black",lwd=2)
  mtext(paste("rho = ",signif(cor.test(Z.K,Z.C,method="s")$estimate,2)," p = ",signif(cor.test(Z.K,Z.C,method="s")$p.value,2),sep=""),cex=0.8,line=0.2)
  abline(v=-2,lty=2,col="grey")
  abline(h=-2,lty=2,col="grey")
  ##blank plot for legend
  par(mar=c(5,5,4,2))
  plot(1, type="n", axes=F, xlab="", ylab="")
  legend(0.6, 1.4, unique(diagnosis[,1]), fill=unique(colorvec))
} #taken from steves expression pipeline

# create functio to ID outliers

names_of_outliers<-function(dataset, threshold){
  datExprs<-t(dataset[2:dim(dataset)[2]])
  IAC = cor(datExprs, method = "p", use = "p")
  diag(IAC) = 0
  A.IAC = ((1 + IAC)/2)^2  ## ADJACENCY MATRIX
  # fundamentalNetworkConcepts
  FNC = fundamentalNetworkConcepts(A.IAC)  ## WGCNA
  K2 = FNC$ScaledConnectivity
  Z.K = round((K2 - mean(K2))/sd(K2), 3)
  Z.K_outliers <- Z.K < threshold
  Z.K_outliers <- names(Z.K_outliers[Z.K_outliers == TRUE])
  n_outliers <- length(Z.K_outliers)
  return(Z.K_outliers)
}

# create function to run network analysis on each expression dataset, plot and remove bad samples 

run_sample_network_plot<-function(dataset, threshold){
  #sample network plot
  sampleNetwork_plot(dataset)
  #identify sample below Z.K threshold
  dataset_removal_1<-names_of_outliers(dataset, threshold)
  # #create empty count list to record samples removed
  count<-dataset_removal_1
  # reiterate over datae till no samples fall below threshold
  if(length(dataset_removal_1)>0){
  while (length(dataset_removal_1)>0) {
    # remove bad samples 
    dataset_QC<-dataset[!(rownames(dataset)%in%dataset_removal_1),]
    #plot
    sampleNetwork_plot(dataset_QC)
    #identify sample below Z.K threshold
    dataset_removal_1<-names_of_outliers(dataset_QC, threshold)
    #record samples removed
    count<-c(count, dataset_removal_1)
    # rename dataset
    dataset<-dataset_QC
  }}
  else {dataset_QC<-dataset}
  # print to screen number of samples removed
  cat("\n")
  print(c("Total number of samples removed...", length(count)))
  # return clean expression set
  return(dataset_QC)
}

# run sample network on full data

exprs_data_sva_drug_QC<-run_sample_network_plot(exprs_data_sva_drug, -10)
exprs_data_sva_control_QC<-run_sample_network_plot(exprs_data_sva_control, -3)

#merge QC data together
any(colnames(exprs_data_sva_control_QC)==colnames(exprs_data_sva_drug_QC))==F

exprs_data_sva_QC<-rbind(exprs_data_sva_control_QC, exprs_data_sva_drug_QC)

dim(exprs_data_sva_QC)
head(exprs_data_sva_QC)[1:5]

##### PLOT SAMPLE NETWORK ANALYSIS TO PDF #####

setwd(sample_network_dir)

pdf("drug_sample_network_analysis_check.pdf")
exprs_data_sva_drug_QC<-run_sample_network_plot(exprs_data_sva_drug, -10)
dev.off()

pdf("control_sample_network_analysis.pdf")
exprs_data_sva_control_QC<-run_sample_network_plot(exprs_data_sva_control, -3)
dev.off()

##### PCA PLOT 5 #####

setwd(PCA_plots_dir)

pdf("5.PCA_plot_after_sample_network_removal_unfiltered.pdf")
pca_data(t(exprs_data_sva_QC[2:dim(exprs_data_sva_QC)[2]]), pheno, 'bottomright')
dev.off()

##### CONVERT TO ENTREZ GENE ID ######

# Get the probe identifiers that are mapped to an ENTREZ Gene ID using hgu133a.db
mapped_probes <- mappedkeys(hgu133aENTREZID)

# Convert to a list
hgu133a.db_mapping <- as.data.frame(hgu133aENTREZID[mapped_probes])
# arrange order of column by entrezgene probe_id
hgu133a.db_mapping<-hgu133a.db_mapping[c(2,1)]
colnames(hgu133a.db_mapping)[1]<-"entrezgene"

head(hgu133a.db_mapping)
dim(hgu133a.db_mapping)

#check any duplicated probe IDs
anyDuplicated(hgu133a.db_mapping$probe_id)

#check any duplicated entrezgene IDs
anyDuplicated(hgu133a.db_mapping$entrezgene)

# add group into list
hgu133a.db_mapping<-rbind(hgu133a.db_mapping, c("group"))

# create convert_probe_id_to_entrez_id function 

convert_probe_id_to_entrez_id <- function(expression_dataset, probe_mapping_file){
  # transform dataset # - removed this step
  # expression_dataset_t<-as.data.frame(expression_dataset)
  # keep only probes which appear in probe_mapping_file
  data_frame_in_probe_mapper<-expression_dataset[colnames(expression_dataset)%in%probe_mapping_file$probe_id]
  # match probe id in data_frame_in_probe_mapper to that in probe_mapping_file and convert to entrez id
  colnames(data_frame_in_probe_mapper)<-probe_mapping_file$entrezgene[match(colnames(data_frame_in_probe_mapper), probe_mapping_file$probe_id)]
  return(data_frame_in_probe_mapper)
}

# apply to exprs_data_sva_QC

exprs_data_sva_QC_entrez_id<-convert_probe_id_to_entrez_id(exprs_data_sva_QC, hgu133a.db_mapping)
dim(exprs_data_sva_QC)
dim(exprs_data_sva_QC_entrez_id)
length(which(duplicated(colnames(exprs_data_sva_QC_entrez_id))))

head(exprs_data_sva_QC_entrez_id)[1:5]


##### COLLAPSE MULTIPPLE ENTREZ ID BY SELECTING ONE WITH HIGHEST AVERAGE EXPRESSION ACROSS SAMPLES ######

select_duplicate_probe_by_top_expr <- function(exprs) {
  # extract gouping
  group<-exprs[1]
  # transpose data frame - keep as dataframe
  exprs_t<-as.data.frame(t(exprs[2:dim(exprs)[2]]))
  # calculate mean expression per probe across samples - create new column - probe mean column
  exprs_t$probe_mean_expression<-rowMeans(exprs_t)
  #copy rownames (probe id) to column and truncate
  exprs_t$trunc_entrez_id<-trunc(as.numeric(as.character(rownames(exprs_t))))
  # order data frame by truncated probe id and then expression level
  exprs_t<-exprs_t[order(exprs_t$trunc_entrez_id, -exprs_t$probe_mean_expression), ]
  # remove all duplicate probe id - keep one with highest mean expression
  exprs_t_unique<-exprs_t[!duplicated(exprs_t$trunc_entrez_id),]
  #unique entrez column back to row name
  rownames(exprs_t_unique)<-exprs_t_unique$trunc_entrez_id
  #remove unwanted column
  exprs_t_unique$trunc_entrez_id<-NULL
  #remove unwanted column
  exprs_t_unique$probe_mean_expression<-NULL
  #transpose dataframe back
  exprs_unique<-as.data.frame(t(exprs_t_unique))
  # bind grouping back
  exprs_unique<-merge(group, exprs_unique, by="row.names")
  #move row.names to rownames
  rownames(exprs_unique)<-exprs_unique$Row.names
  #remove unwanted row.names
  exprs_unique$Row.names<-NULL
  return(exprs_unique)
}

# apply function to exprs_data_sva_QC_entrez_id

exprs_data_sva_QC_entrez_id_unique<-select_duplicate_probe_by_top_expr(exprs_data_sva_QC_entrez_id)
dim(exprs_data_sva_QC_entrez_id_unique)
length(which(duplicated(colnames(exprs_data_sva_QC_entrez_id_unique))))

head(exprs_data_sva_QC_entrez_id_unique[1:5])

##### SAVE CLEAN DATA IN ENTREZ GENE FORMAT #####

setwd(clean_dir)

save(exprs_data_sva_QC_entrez_id_unique, file="exprs_data_sva_QC_entrez_id.Rdata")

write.csv(subset_pheno, file="subset_pheno_info.txt", row.names=T)

##### DIFFERENTIAL EXPRESSION ANALYSIS - FULL PROBES #####

#pair up same drugs and same concentrations and run through limma

#subset pheno 
pheno_in_full_data<-subset(subset_pheno, cell_id %in% rownames(exprs_data_sva_QC_entrez_id_unique))
dim(pheno_in_full_data)[1]==dim(exprs_data_sva_QC_entrez_id_unique)[1]

# create table of count
full_drug_count<-as.data.frame(table(pheno_in_full_data$drug_conc))
#change colnames
colnames(full_drug_count)[1]<-"Drug"
head(full_drug_count)
# check each is replicated at least once
any(full_drug_count$Freq==1)
#remove drug from list
full_drug_count<-subset(full_drug_count, Drug!="Control")
#check control removed
grep("Control", full_drug_count$Drug)
# get control cell_ids
control_cell_id<-(subset(pheno_in_full_data, group=="Control"))$cell_id
length(control_cell_id)

#create cluster to run limma
# Calculate the number of cores
no_cores <- detectCores() - 1
registerDoParallel(no_cores)

#create function

full_DE_results<-foreach(x=1:dim(full_drug_count)[1]) %dopar% {
  #extract drug
  drug_to_extract<-full_drug_count$Drug[x]
  #extract cell id
  cell_id_to_extract<-(subset(pheno_in_full_data, drug_conc %in% drug_to_extract))$cell_id
  # add controls to this gourp
  cell_id_to_extract<-c(cell_id_to_extract, control_cell_id)
  #subset expression table
  subset_exprs<-subset(exprs_data_sva_QC_entrez_id_unique, rownames(exprs_data_sva_QC_entrez_id_unique) %in% cell_id_to_extract)
  #split exprs by extract group
  group<-subset_exprs[,1]
  exprs<-subset_exprs[2:dim(subset_exprs)[2]]
  #replace sample name to numbers
  rownames(exprs)<-c(1:dim(exprs)[1])
  # set up design
  design <- model.matrix(~0 + group)
  # change colnames to Control + Drug
  colnames(design)<-c("Control", "Drug")
  # transpose dataset, convert to numeric 
  transposed_exprs<-t(exprs)
  #run limma
  dataset_exprs_fit <- lmFit(transposed_exprs, design, method="robust")
  exprs_contrast_matrix<- makeContrasts(Drug-Control, levels=design)
  exprs_contrast_fit <-contrasts.fit(dataset_exprs_fit, exprs_contrast_matrix)
  exprs_contrast_ebayes <- eBayes(exprs_contrast_fit, robust=T)
  exprs_top_genes <- topTable(exprs_contrast_ebayes, number=(dim(transposed_exprs)[1]), coef=1, adjust.method="fdr", confint=TRUE) 
  return(exprs_top_genes)
}

#stop cluster
stopImplicitCluster()

# add names of drugs to list
for(x in 1:dim(full_drug_count)[1]){
  names(full_DE_results)[x]<-as.character(full_drug_count$Drug[x])
}

#check output file
names(full_DE_results)
head(full_DE_results[[1]])
head(full_DE_results[[2]])

#save
setwd(DE_dir)

save(full_DE_results, file="full_DE_results.Rdata")

##### SIG DE ANALYSIS COUNT #####

grep("trichostatin A", names(full_DE_results))

head(full_DE_results[[39]])

dim(subset(full_DE_results[[39]], adj.P.Val<=0.05))

dim(full_DE_results[[39]])

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
dim(exprs_data_sva_QC_entrez_id_unique)
head(number_of_DEG_in_full_data)
tail(number_of_DEG_in_full_data)

subset(number_of_DEG_in_full_data, sig_DEG_count=="0")

##### SAVE IMAGE #####

setwd(processing_dir)
save.image(paste(data_name, "_processing.Rdata", sep=""))
#load(paste(data, "_processing.Rdata", sep=""))
