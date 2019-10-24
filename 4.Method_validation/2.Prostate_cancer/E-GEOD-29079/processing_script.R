
##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                              E-GEOD-29079 PROCESSING                                                                   #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################

# QC PIPELINE VERSION: 2.2 single dataset, single tissue
# DATE: 26/02/2018
# ARRAY EXPRESS NUMBER: E-GEOD-29079
# DISORDER: prostate cancer
# MICROARRAY PLATFORM: Affymetrix
# EXPRESSION CHIP: HG Exon 1.0 ST
# NUMBER OF SAMPLES: 95
# TISSUE: 
#
# NOTES - 
# script modified as all males
# using pre-processed data - already RMA bg and log2
# Not enough genes to verify gene detection threshold - 

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### LOAD LIBRARIES ####

library(Biobase)
library(GEOquery)
library(RCurl)
library(XML)
library(ArrayExpress)
library(affy)
library(lumi)
library(WGCNA)
library(pamr)
library(sva)
library(ggplot2)
library(reshape)
library(massiR)
library(limma)
library(gridExtra)

#library(pd.hg.u133.plus.2)
#library(hgu133plus2.db)
#library(hgu133a.db)
#library(hgu133b.db)
library(hugene10sttranscriptcluster.db)
#library(illuminaHumanv4.db)
#library(illuminaHumanv3.db)

##### SET DIRECTORIES ####

work_dir="/media/hamel/Workspace/Dropbox/Projects/CMap/1.Data/Validation_datasets/Prostate_Cancer"

setwd(work_dir)

# create directory for raw data
dir.create(paste(work_dir,"Raw_Data", sep="/"))
raw_dir=paste(work_dir,"Raw_Data", sep="/")

# create directory for all Plots
dir.create(paste(work_dir,"Preprocessing_Plots", sep="/"))
plots_dir=paste(work_dir,"Preprocessing_Plots", sep="/")

# create directory for boxplots + density plots
dir.create(paste(plots_dir,"Boxplots_Density_plots", sep="/"))
boxplots_density_plots_dir=paste(plots_dir,"Boxplots_Density_plots", sep="/")

# create directory for PCA plots
dir.create(paste(plots_dir,"PCA_Plots", sep="/"))
pca_dir=paste(plots_dir,"PCA_Plots", sep="/")

# create directory for sample netowk plots
dir.create(paste(plots_dir,"Sample_Network_Plots", sep="/"))
sample_network_dir=paste(plots_dir,"Sample_Network_Plots", sep="/")

# create directory non-expressed threshold plots
dir.create(paste(plots_dir,"Probe_Detection_Plots", sep="/"))
probe_detection_plots_dir=paste(plots_dir,"Probe_Detection_Plots", sep="/")

# create directory for QC's expression data
dir.create(paste(work_dir,"Clean_Data", sep="/"))
clean_data_dir=paste(work_dir,"Clean_Data", sep="/")

# create directory Differential Expression Results
dir.create(paste(work_dir,"DE_results", sep="/"))
DE_dir=paste(work_dir,"DE_results", sep="/")

##### TEMP FIX FOR getAE FUNCTION ######

# EBI changed to https - their getAE function broken. here is temp fix

https_getAE <- function (accession, path = getwd(), type = "full", extract = TRUE, local = FALSE, sourcedir = path) {
  if (!local) {
    baseURL = "https://www.ebi.ac.uk/arrayexpress/xml/v2/files"
    xmlURL = getURL(paste(baseURL, accession, sep = "/"))
    xml = xmlTreeParse(xmlURL, useInternalNodes = TRUE, isURL=FALSE)
    sdrfURL = xpathSApply(xml, "/files/experiment/file[kind='sdrf' and extension='txt']/url", 
                          xmlValue)
    sdrfFile = xpathSApply(xml, "/files/experiment/file[kind='sdrf' and extension='txt']/name", 
                           xmlValue)
    idfURL = xpathSApply(xml, "/files/experiment/file[kind='idf' and extension='txt']/url", 
                         xmlValue)
    idfFile = xpathSApply(xml, "/files/experiment/file[kind='idf' and extension='txt']/name", 
                          xmlValue)
    adfURL = xpathApply(xml, "/files/experiment/file[kind='adf' and extension='txt']/url", 
                        xmlValue)
    adfFiles = xpathApply(xml, "/files/experiment/file[kind='adf' and extension='txt']/name", 
                          xmlValue)
    rawArchiveURL = xpathApply(xml, "/files/experiment/file[kind='raw' and extension='zip']/url", 
                               xmlValue)
    procArchiveURL = xpathApply(xml, "/files/experiment/file[kind='processed' and extension='zip']/url", 
                                xmlValue)
  }
  else {
    allfiles = list.files(sourcedir)
    sdrfFile = allfiles[grep(paste(accession, ".sdrf.txt$", 
                                   sep = ""), allfiles)]
    if (length(sdrfFile) == 0) 
      stop("SDRF file not found in directory ", sourcedir)
    sdrfURL = paste("file:/", sourcedir, sdrfFile, sep = "/")
    idfFile = allfiles[grep(paste(accession, ".idf.txt$", 
                                  sep = ""), allfiles)]
    if (length(idfFile) == 0) 
      warning("IDF file not found in directory ", sourcedir)
    idfURL = paste("file:/", sourcedir, idfFile, sep = "/")
    ph = try(read.AnnotatedDataFrame(sdrfFile, path = sourcedir, 
                                     row.names = NULL, blank.lines.skip = TRUE, fill = TRUE, 
                                     varMetadata.char = "$"))
    if (inherits(ph, "try-error")) {
      warning("Unable to retrieve ADF reference from SDRF. Reading any ADF in directory.")
      adfFiles = allfiles[grep(".adf.txt$", allfiles)]
    }
    else {
      adr = unique(pData(ph)[, getSDRFcolumn("ArrayDesignREF", 
                                             varLabels(ph))])
      adfFiles = paste(adr, ".adf.txt", sep = "")
    }
    if (all(file.exists(file.path(sourcedir, adfFiles)))) {
      adfURL = paste("file:/", sourcedir, adfFiles, sep = "/")
      downloadADF = FALSE
    }
    else {
      filesURL = "https://www.ebi.ac.uk/arrayexpress/files"
      adfURL = paste(filesURL, adr, adfFiles, sep = "/")
      downloadADF = TRUE
    }
    rawArchiveURL = NULL
    procArchiveURL = NULL
    rawArchive = allfiles[grep(paste(accession, ".raw.[0-9]{1,}.zip", 
                                     sep = ""), allfiles)]
    if (length(rawArchive) != 0) 
      rawArchiveURL = paste("file:/", sourcedir, rawArchive, 
                            sep = "/")
    else warning("No raw files found in directory ", sourcedir)
    processedArchive = allfiles[grep(paste(accession, ".processed.[0-9]{1,}.zip", 
                                           sep = ""), allfiles)]
    if (length(processedArchive) != 0) 
      procArchiveURL = paste("file:/", sourcedir, processedArchive, 
                             sep = "/")
    else warning("No processed data files found in directory ", 
                 sourcedir)
  }
  if (length(sdrfURL) > 1) {
    warning("Found two SDRF files: \n", paste(sdrfURL, "\n"))
    hybSDRF = grep("hyb.sdrf", sdrfURL)
    if (length(hybSDRF) > 0) {
      message("Choosing ", sdrfURL[hybSDRF])
      sdrfURL = sdrfURL[hybSDRF]
      sdrfFile = sdrfFile[hybSDRF]
    }
    else {
      warning("Unable to choose SDRF file. Please report experiment to miamexpress@ebi.ac.uk")
    }
  }
  if (!local || path != sourcedir || downloadADF) {
    adfFiles <- lapply(adfURL, function(url) {
      filedest = paste(path, basename(url), sep = "/")
      dnld = try(download.file(url, filedest, mode = "wb"))
      if (inherits(dnld, "try-error") || file.info(filedest)$size == 
          0) {
        warning(paste(url, " does not exist or is empty. \n"), 
                sep = "")
        adffile = NULL
      }
      else {
        adffile = basename(filedest)
      }
      return(adffile)
    })
    if (!is.null(adfFiles)) 
      adfFiles = unlist(adfFiles)
  }
  if (!local || path != sourcedir) {
    sdrfFileDest = paste(path, sdrfFile, sep = "/")
    dnld = try(download.file(sdrfURL, sdrfFileDest, mode = "wb"))
    if (inherits(dnld, "try-error") || file.info(sdrfFileDest)$size == 
        0) {
      warning(paste(sdrfFile, " does not exist or is empty. The object will not have featureData or phenoData. \n"), 
              sep = "")
      sdrfFile = NULL
      adffile = NULL
    }
    idfFileDest = paste(path, idfFile, sep = "/")
    dnld = try(download.file(idfURL, idfFileDest, mode = "wb"))
    if (inherits(dnld, "try-error") || file.info(idfFileDest)$size == 
        0) {
      warning(paste(idfFile, " does not exist or is empty. \n"), 
              sep = "")
      idfFile = NULL
    }
    rawArchive = NULL
    processedArchive = NULL
    if (type != "mageFilesOnly" && !is.null(rawArchiveURL) && 
        (type == "full" || type == "raw")) {
      message("Copying raw data files\n")
      rawArchive <- lapply(rawArchiveURL, function(url) {
        filedest = paste(path, basename(url), sep = "/")
        dnld = try(download.file(url, filedest, mode = "wb"))
        if (inherits(dnld, "try-error") || file.info(filedest)$size == 
            0) {
          warning(paste(url, " does not exist or is empty. \n"), 
                  sep = "")
        }
        else {
          return(filedest)
        }
      })
      if (!is.null(rawArchive)) {
        rawArchive = unlist(rawArchive)
        rawArchive = basename(rawArchive)
      }
    }
    if ((type != "mageFilesOnly" && type == "full" || type == 
         "processed") && !is.null(procArchiveURL)) {
      message("Copying processed data files\n")
      processedArchive <- lapply(procArchiveURL, function(url) {
        filedest = paste(path, basename(url), sep = "/")
        dnld = try(download.file(url, filedest, mode = "wb"))
        if (inherits(dnld, "try-error") || file.info(filedest)$size == 
            0) {
          warning(paste(url, " does not exist or is empty. \n"), 
                  sep = "")
        }
        else {
          return(filedest)
        }
      })
      if (!is.null(processedArchive)) {
        processedArchive = unlist(processedArchive)
        processedArchive = basename(processedArchive)
      }
    }
  }
  rawFiles = NULL
  processedFiles = NULL
  if (extract) {
    message("Unpacking data files")
    if (!is.null(rawArchive)) 
      rawFiles <- lapply(rawArchive, function(zipfile) {
        rawfiles = extract.zip(file = paste(path, zipfile, 
                                            sep = "/"))
        return(rawfiles)
      })
    if (!is.null(processedArchive)) 
      processedFiles <- lapply(processedArchive, function(zipfile) {
        procfiles = extract.zip(file = paste(path, zipfile, 
                                             sep = "/"))
        return(procfiles)
      })
    if (!is.null(rawFiles)) 
      rawFiles = unlist(rawFiles)
    if (!is.null(processedFiles)) 
      processedFiles = unlist(processedFiles)
  }
  res = list(path = path, rawFiles = rawFiles, rawArchive = rawArchive, 
             processedFiles = processedFiles, processedArchive = processedArchive, 
             sdrf = sdrfFile, idf = idfFile, adf = adfFiles)
  return(res)
}

##### DOWNLOAD RAW DATA #####

setwd(raw_dir)

#raw data only - https
data_raw=https_getAE("E-GEOD-29079", type = "full")

#raw data only
#data_raw=getAE("E-GEOD-29079", type = "raw")

#processed data only
#data_raw=getAE("E-GEOD-29079", type = "processed")

#all data
#data_raw=getAE("E-GEOD-29079", type = "full")

##### CREATE R EXPRESSION OBJECT FROM RAW DATA #####

# METHOD 1 - convert MAGE-TAB files into expresssion set  - USING RAW DATA

expression_data = ae2bioc(mageFiles = data_raw)

expression_data

# # METHOD 2 - convert MAGE-TAB files into expresssion set - USING RAW DATA
# 
#expression_data<-ReadAffy()
# 
# METHOD 3 - convert MAGE-TAB files into expresssion set - USING PROCESSED DATA
# 
# cnames=getcolproc(data_raw)
# # 
cnames
#
expression_data=procset(data_raw, cnames[2])
#
expression_data

##### SET DATA PARAMETERS #####

##
## dataset name to save as/use
##

dataset="E-GEOD-29079"

##
## disease
##

disease="Prostate_Cancer"

##
## Affymetrix or Illumina
##

Microarray_platform="Affymetrix"

##
## Raw or pre-processed
##

Data_format="Raw"

##
## probe detection thresohld to use
##

#Probe_Detection_Threshold=0.9
#Probe_Detection_Threshold=0.8 
#Probe_Detection_Threshold=0.7 
Probe_Detection_Threshold=0.4 # not enough geens to verify cut-off threshold

##
## expression chip to use 
##

#expression_chip="hgu133plus2"
#expression_chip="hgu133a"
#expression_chip="hgu133b"
expression_chip="hugene10sttranscriptcluster"
#expression_chip="illuminaHumanv4"
#expression_chip="illuminaHumanv3"

##
## sample network threshold to use - samples less than Z.K threshold will be removed - need manual check
##

sample_network_ZK_threshold=-3

##
## massi R chip to use
##

#massi_R_chip="illumina_humanwg_6_v1" 
#massi_R_chip="illumina_humanwg_6_v2" 
#massi_R_chip="illumina_humanwg_6_v1" 
#massi_R_chip="illumina_humanht_12"   
massi_R_chip="affy_hugene_1_0_st_v1" 
#massi_R_chip="affy_hg_u133_plus_2"  

##
## Phenotype info
##

phenotype_data<-pData(expression_data)
head(phenotype_data)
names(phenotype_data)

Ethnicity<-phenotype_data[1]
colnames(Ethnicity)<-"Ethnicity"

Tissue<-phenotype_data[3]
colnames(Tissue)<-"Tissue"

Age<-phenotype_data[1]
colnames(Age)<-"Age"

Diagnosis<-phenotype_data[38]
colnames(Diagnosis)<-"Diagnosis"

Gender<-phenotype_data[1]
colnames(Gender)<-"Gender"

#standardise Ethnicity
Ethnicity$Ethnicity<-as.character(Ethnicity$Ethnicity)

table(Ethnicity)
Ethnicity[1]<-"unknown"
table(Ethnicity)

#standardise Age

table(Age)
Age[1]<-"unknown"
table(Age)

#standardise diagnosis- case control

Diagnosis$Diagnosis<-as.character(Diagnosis$Diagnosis)

table(Diagnosis)
Diagnosis[grep("benign", Diagnosis[,1]),]<-"control"
Diagnosis[grep("Tumor", Diagnosis[,1]),]<-"case"
table(Diagnosis)

#standardise tissue

Tissue$Tissue<-as.character(Tissue$Tissue)

table(Tissue)

#standardise Gender

Gender$Gender<-as.character(Gender$Gender)

table(Gender)
Gender[1]<-"male"
table(Gender)

##### PLOTS OF RAW DATA #####

setwd(boxplots_density_plots_dir)

boxplot(head(exprs(expression_data)))
pdf(file="raw_data_boxplot.pdf")
boxplot(exprs(expression_data))
dev.off()

plotDensity(head(exprs(expression_data)), logMode=F, addLegend=F)
pdf(file="raw_data_density_plot.pdf")
plotDensity(exprs(expression_data), logMode=F, addLegend=F)
dev.off()

##### PRE-PROCESS ######

library(oligo)

#normalise

expression_data_background_corrected <- rma(expression_data, target = "core")

#background correct
# 
# expression_data_background_corrected<-mas5(expression_data)
# 
# #normalise
# 
expression_data_normalised<-rsn(exprs(expression_data_background_corrected))
# #expression_data_normalised<-rsn(log2(exprs(expression_data)))

# set negative values to zero

#expression_data_normalised<-exprs(expression_data)
expression_data_normalised[expression_data_normalised<0]<-1

#convert to data.frame

expression_data_normalised_as_data_frame<-as.data.frame(expression_data_normalised)

##### PLOTS OF PRE_PROCESSED DATA #####

setwd(boxplots_density_plots_dir)

boxplot(expression_data_normalised)
pdf(file="pre-processed_data_boxplot.pdf")
boxplot(expression_data_normalised)
dev.off()

plotDensity(expression_data_normalised, logMode=F, addLegend=F)
pdf(file="pre-processed_data_density_plot.pdf")
plotDensity(expression_data_normalised, logMode=F, addLegend=F)
dev.off()

setwd(work_dir)

##### CHECK FOR DUPLICATE SAMPLES IDs #####

anyDuplicated(rownames(phenotype_data))

##### PCA PLOT 1 #####

pca_data<-function(data, legend_position){
  #run PCA
  pca<-prcomp(data)
  # order of samples in expression data
  sample_order<-colnames(data)
  # merge pheno info together
  pheno<-cbind(Diagnosis, Gender, Ethnicity, Age)
  # subset pca_pheno to match expression data
  pca_pheno<-subset(pheno, rownames(pheno) %in% colnames(data))
  # match order
  ordered_pca_pheno<-pca_pheno[match(sample_order, rownames(pca_pheno)),]
  Diagnosis_pca_colour<-labels2colors(as.character(ordered_pca_pheno$Diagnosis))
  Gender_pca_colour<-labels2colors(as.character(ordered_pca_pheno$Gender))
  Ethnicity_pca_colour<-labels2colors(as.character(ordered_pca_pheno$Ethnicity))
  Age_pca_colour<-labels2colors(as.character(ordered_pca_pheno$Age))
  # pca plot - Diagnosis
  plot(pca$rotation[,1:2], main=" PCA plot coloured by Diagnosis",col="black", pch=21,bg=Diagnosis_pca_colour)
  legend(legend_position, legend=unique(ordered_pca_pheno$Diagnosis), fill=unique(Diagnosis_pca_colour), title="Diagnosis")
  # pca plot - Gender
  plot(pca$rotation[,1:2], main=" PCA plot coloured by Clinical Gender",col="black", pch=21,bg=Gender_pca_colour)
  legend(legend_position, legend=unique(ordered_pca_pheno$Gender), fill=unique(Gender_pca_colour), title="Gender")
  # pca plot - Ethnicity
  plot(pca$rotation[,1:2], main=" PCA plot coloured by Ethnicity",col="black", pch=21,bg=Ethnicity_pca_colour)
  legend(legend_position, legend=unique(ordered_pca_pheno$Ethnicity), fill=unique(Ethnicity_pca_colour), title="Ethnicity")
  # pca plot - Age
  plot(pca$rotation[,1:2], main=" PCA plot coloured by Age",col="black", pch=21,bg=Age_pca_colour)
  legend(legend_position, legend=unique(ordered_pca_pheno$Age), fill=unique(Age_pca_colour), title="Age")
}

#apply function
pca_data(expression_data_normalised_as_data_frame, 'bottomright')

#plot to pdf
setwd(pca_dir)

pdf("1.PCA_plot_before_QC.pdf")
pca_data(expression_data_normalised_as_data_frame, 'bottomright')
dev.off()

# ##### GENDER CHECK ##### 
# 
# # get Y choromosome genes
# data(y.probes)
# names(y.probes)
# 
# y_chromo_probes <- data.frame(y.probes[massi_R_chip])
# 
# # extract Y chromosome genes from dataset
# eset.select.out <- massi_select(expression_data_normalised_as_data_frame, y_chromo_probes)
# 
# #massi_y_plot(eset.select.out)
# massi_cluster_plot(eset.select.out)
# 
# # check sex bias - should have at least 15% male samples and minimum 6 samples
# dip.result <- massi_dip(eset.select.out)
# 
# # run gender predict
# eset.results <- massi_cluster(eset.select.out)
# 
# #extract gender prediction
# predicted_gender<-(eset.results$massi.results)[c(1,5)]
# rownames(predicted_gender)<-predicted_gender$ID
# predicted_gender$ID<-NULL
# colnames(predicted_gender)<-"Predicted_Gender"
# 
# #compare to clinical Gender
# 
# # merge
# gender_comparison<-merge(Gender, predicted_gender, by="row.names")
# rownames(gender_comparison)<-gender_comparison$Row.names
# gender_comparison$Row.names<-NULL
# colnames(gender_comparison)<-c("Clinical_Gender", "Predicted_Gender")
# head(gender_comparison)
# 
# # dip test <0.08 - sex bias - change gender to unknown
# gender_comparison<-Gender
# gender_comparison$Predicted_Gender<-"Unknown"
# colnames(gender_comparison)[1]<-"Clinical_Gender"
# 
# # #separae male/female IDs - use predicted
# # female_samples<-subset(gender_comparison, Predicted_Gender=="female")
# # male_samples<-subset(gender_comparison, Predicted_Gender=="male")
# #
# # head(female_samples)
# # head(male_samples)
# #
# # # gender miss-matches
# Gender_Missmatch<-gender_comparison[gender_comparison$Clinical_Gender!=gender_comparison$Predicted_Gender,]
# 
# Gender_Missmatch

gender_comparison<-cbind(Gender, Gender)
colnames(gender_comparison)<-c("Clinical_Gender", "Predicted_Gender")
gender_comparison$Predicted_Gender<-"Unknown"
head(gender_comparison)

##### PROBE ID DETECTION #####

# separate case control - Factor.Value..disease. column 

case_ID<-rownames(subset(Diagnosis, Diagnosis=="case"))

control_ID<-rownames(subset(Diagnosis, Diagnosis=="control"))

case_ID

control_ID

case_exprs<-expression_data_normalised_as_data_frame[,colnames(expression_data_normalised_as_data_frame)%in%case_ID]
control_exprs<-expression_data_normalised_as_data_frame[,colnames(expression_data_normalised_as_data_frame)%in%control_ID]

head(case_exprs)
dim(case_exprs)
head(control_exprs)
dim(control_exprs)

# calculate 90th percentile for each sample in each group

extract_good_probe_list<-function(dataset, probe_percentile_threshold) {
  # dataset - expression dataset as dataframe
  # probe_percentile_threshold - percentile at which to use as cut-off for detected probes 
  # number of samples in which probe must be expressed in - fixed at 0.8 - i.e 80% of samples
  # calculate quantile threshold for each sample
  sample_quantiles<-apply(dataset, 2, quantile, probs=c(probe_percentile_threshold))
  # count length of quantile - will be number of samples
  number_of_samples<-length(sample_quantiles)
  # convert probes values to NA in for each sample if probe expression value below sample percentile cut-off
  for (x in 1:number_of_samples) {
    is.na(dataset[x]) <- dataset[x] >= sample_quantiles[x]
  }
  # convert to dataframe
  dataset_dataframe<-as.data.frame(dataset)
  # count number of NA
  dataset_count<-as.data.frame(rowSums(is.na(dataset_dataframe)))
  colnames(dataset_count)<-"count"
  # subset good probes
  good_probes<-rownames(subset(dataset_count, dataset_count$count >= (number_of_samples*0.8)))
  #print threshold used
  print(as.data.frame(sample_quantiles))
  boxplot(as.data.frame(sample_quantiles))
  # return good probes
  return(good_probes)
}

# apply function to case samples

case_exprs_expressed_probes_list<-extract_good_probe_list(case_exprs, Probe_Detection_Threshold)
length(case_exprs_expressed_probes_list)

control_exprs_expressed_probes_list<-extract_good_probe_list(control_exprs, Probe_Detection_Threshold)
length(control_exprs_expressed_probes_list)


# #apply by case control on
# case_exprs_expressed_probes_list<-extract_good_probe_list(case_exprs, Probe_Detection_Threshold)
# length(case_exprs_expressed_probes_list)
# 
# control_exprs_expressed_probes_list<-extract_good_probe_list(control_exprs, Probe_Detection_Threshold)
# length(control_exprs_expressed_probes_list)


# merge list of good probes from both case + control, sort and keep unique values

good_probe_list<-unique(sort(c(case_exprs_expressed_probes_list,
                               control_exprs_expressed_probes_list)))

length(good_probe_list)

# extract good probes from dataset

data_exprs_good_probes<-expression_data_normalised_as_data_frame[rownames(expression_data_normalised_as_data_frame)%in%good_probe_list,]

data_case_exprs_good_probes<-case_exprs[rownames(case_exprs)%in%good_probe_list,]

data_control_exprs_good_probes<-control_exprs[rownames(control_exprs)%in%good_probe_list,]

head(data_exprs_good_probes)[1:5]
dim(expression_data_normalised)
dim(data_exprs_good_probes)
dim(data_case_exprs_good_probes)
dim(data_control_exprs_good_probes)

##### PROBE DETECTION THRESHOLD PLOTS #####

# using dataframe before probe removal

# get gene symbol list for chip

Gene_symbols_probes <- mappedkeys(eval(parse(text = paste(expression_chip, "SYMBOL", sep=""))))

# Convert to a list
Gene_symbols <- as.data.frame(eval(parse(text = paste(expression_chip, "SYMBOL", sep="")))[Gene_symbols_probes])

head(Gene_symbols)
dim(Gene_symbols)

#Expressed in females only
XIST_probe_ID<-subset(Gene_symbols, symbol=="XIST")
XIST_probe_ID

#Expressed in Males only
PRKY_probe_ID<-subset(Gene_symbols, symbol=="PRKY")
PRKY_probe_ID

RPS4Y1_probe_ID<-subset(Gene_symbols, symbol=="RPS4Y1")
RPS4Y1_probe_ID

KDM5D_probe_ID<-subset(Gene_symbols, symbol=="KDM5D")
KDM5D_probe_ID

USP9Y_probe_ID<-subset(Gene_symbols, symbol=="USP9Y")
USP9Y_probe_ID

UTY_probe_ID<-subset(Gene_symbols, symbol=="UTY")
UTY_probe_ID

# HK genes expressed in all cells + males + females
MKRN1_probe_ID<-subset(Gene_symbols, symbol=="MKRN1")
MKRN1_probe_ID

ADIPOR1_probe_ID<-subset(Gene_symbols, symbol=="ADIPOR1")
ADIPOR1_probe_ID

BNIP3L_probe_ID<-subset(Gene_symbols, symbol=="BNIP3L")
BNIP3L_probe_ID

# RNF10_probe_ID<-subset(Gene_symbols, symbol=="RNF10")
# RNF10_probe_ID


# merge all genes 
gene_list<-rbind(XIST_probe_ID,
                 PRKY_probe_ID,
                 RPS4Y1_probe_ID,
                 KDM5D_probe_ID,
                 USP9Y_probe_ID,
                 UTY_probe_ID,
                 MKRN1_probe_ID,
                 ADIPOR1_probe_ID,
                 BNIP3L_probe_ID)
                 #RNF10_probe_ID)

gene_list

# create table of genes and state if KH or gender specific
gene_table<-read.table(text =
                         "Gene Expressed_in
                       ADIPOR1 All
                       BNIP3L All
                       KDM5D Males
                       MKRN1 All
                       PRKY Males
                       RPS4Y1 Males
                       USP9Y Males
                       UTY Males
                       XIST Females" , header=T)

gene_table

#create function to plot
plot_gender_specific_genes<-function(Expression_table, Gender, genes_to_extract, threshold, boxplot_title){
  #extract gene of interest
  Expression_table_gene_check<-as.data.frame(t(Expression_table[rownames(Expression_table)%in% genes_to_extract$probe_id,]))
  # change colnames TO GENE SYMBOL using genes to extract file
  for (x in 1:dim(Expression_table_gene_check)[2]){
    colnames(Expression_table_gene_check)[x]<-gene_list[genes_to_extract$probe_id==colnames(Expression_table_gene_check)[x],2]
  }
  # add in gender information
  Expression_table_gene_check_gender<-merge(Gender, Expression_table_gene_check, by="row.names")
  rownames(Expression_table_gene_check_gender)<-Expression_table_gene_check_gender$Row.names
  Expression_table_gene_check_gender$Row.names<-NULL
  #melt dataframe for plot
  Expression_table_gene_check_gender_melt<-melt(Expression_table_gene_check_gender, by=Gender)
  # change variable colun from factor to character
  Expression_table_gene_check_gender_melt$variable<-as.character(Expression_table_gene_check_gender_melt$variable)
  # order dataframe by variable
  Expression_table_gene_check_gender_melt<-Expression_table_gene_check_gender_melt[order(Expression_table_gene_check_gender_melt$variable),]
  # calculate user defined percentie threshold
  sample_quantiles<-apply(Expression_table, 2, quantile, probs=threshold)
  # mean of used defined threshold across samples
  mean_threshold=mean(sample_quantiles)
  #plot
  plot1<-qplot(variable, value, colour=get(colnames(Gender)), data = Expression_table_gene_check_gender_melt, geom = c("boxplot", "jitter")) + 
    geom_hline(yintercept = mean_threshold) +
    ggtitle(boxplot_title) +
    theme(text = element_text(size=20), axis.text.x = element_text(angle=45, hjust=1)) +
    labs(x="Gene",y="Expression", colour = colnames(Gender)) 
  # 2nd legend
  plot2<-tableGrob(gene_table, rows=NULL)
  # plot
  grid.arrange(plot1, plot2,
               nrow=2,
               heights=c(3,1))
}
# plot

plot_gender_specific_genes(case_exprs, gender_comparison[1], gene_list, Probe_Detection_Threshold, paste(dataset, "case_samples", sep="_"))
plot_gender_specific_genes(control_exprs, gender_comparison[1], gene_list, Probe_Detection_Threshold, paste(dataset, "control_samples", sep="_"))

setwd(probe_detection_plots_dir)

pdf("Probe_detection_threshold_based_on_gender_specific_and_house_keeping_genes.pdf", height=15, width = 12)
plot_gender_specific_genes(case_exprs, gender_comparison[1], gene_list, Probe_Detection_Threshold, paste(dataset, "case_samples", sep="_"))
plot_gender_specific_genes(control_exprs, gender_comparison[1], gene_list, Probe_Detection_Threshold, paste(dataset, "control_samples", sep="_"))
dev.off()

##### PCA PLOT 2 #####

#plot to pdf
setwd(pca_dir)

pca_data(data_exprs_good_probes, 'bottomright')

pdf("2.PCA_plot_after_probe_detection.pdf")
pca_data(data_exprs_good_probes, 'bottomright')
dev.off()

##### SVA ####

# add diagnosis in

dim(data_case_exprs_good_probes)
dim(data_control_exprs_good_probes)

data_case_exprs_good_probes<-cbind(Diagnosis ="case", as.data.frame(t(data_case_exprs_good_probes)))
data_control_exprs_good_probes<-cbind(Diagnosis = "control", as.data.frame(t(data_control_exprs_good_probes)))

head(data_case_exprs_good_probes)[1:5]
head(data_control_exprs_good_probes)[1:5]

# check order of gene in each dataframe

any(colnames(data_case_exprs_good_probes)==colnames(data_control_exprs_good_probes))==F

# merge case control in to one dataset per data region

data_exprs_good_probes_diagnosis<-rbind(data_case_exprs_good_probes, data_control_exprs_good_probes)
head(data_exprs_good_probes_diagnosis)[1:5]

# create sva function

check_SV_in_data<-function(dataset){
  # create sva compatable matrix - sample in columns, probes in rows - pheno info seperate - sort by diagnosis 1st to keep AD top
  sorted_by_diagnosis<-dataset[order(dataset$Diagnosis),]
  # separate expresion and pheno
  dataset_pheno<-sorted_by_diagnosis[1]
  dataset_exprs<-t(sorted_by_diagnosis[2:dim(sorted_by_diagnosis)[2]])
  #full model matrix for Diagnosis
  mod = model.matrix(~Diagnosis, data=dataset_pheno)
  # check number of SV in data
  print(num.sv(dataset_exprs, mod, method="leek"))
}

# check sv

number_of_SV<-check_SV_in_data(data_exprs_good_probes_diagnosis)

# create function to sva adjust

# create function to adjust for SVA and adjust if needed - exclude gender as creating massive mbatch effect

if (number_of_SV>0){
  adjust_for_sva<-function(dataset){
    # create sva compatable matrix - sample in columns, probes in rows - pheno info seperate - sort by diagnosis 1st to keep AD top
    sorted_by_diagnosis<-dataset[order(dataset$Diagnosis),]
    # separate expresion and pheno
    dataset_sva_pheno<-sorted_by_diagnosis[1]
    dataset_sva_exprs<-t(sorted_by_diagnosis[2:dim(sorted_by_diagnosis)[2]])
    #full model matrix for Diagnosis
    mod = model.matrix(~Diagnosis , data=dataset_sva_pheno)
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
  data_exprs_good_probes_diagnosis_sva<-adjust_for_sva(data_exprs_good_probes_diagnosis)
  check_SV_in_data(data_exprs_good_probes_diagnosis_sva)
}  else {
  data_exprs_good_probes_diagnosis_sva<-data_exprs_good_probes_diagnosis
  check_SV_in_data(data_exprs_good_probes_diagnosis_sva)
}

check_SV_in_data(data_exprs_good_probes_diagnosis_sva)

##### PCA PLOT 3 #####

setwd(pca_dir)

pca_data(as.data.frame(t(data_exprs_good_probes_diagnosis_sva[3:dim(data_exprs_good_probes_diagnosis_sva)[2]])), 'bottomright')

pdf("3.PCA_plot_after_SVA.pdf")
pca_data(as.data.frame(t(data_exprs_good_probes_diagnosis_sva[3:dim(data_exprs_good_probes_diagnosis_sva)[2]])), 'bottomright')
dev.off()

##### SAMPLE NETWORK PLOT #####

# separate dataframe into case and control
data_exprs_good_probes_diagnosis_sva_case<-data_exprs_good_probes_diagnosis_sva[data_exprs_good_probes_diagnosis_sva$Diagnosis=="case",]

data_exprs_good_probes_diagnosis_sva_control<-data_exprs_good_probes_diagnosis_sva[data_exprs_good_probes_diagnosis_sva$Diagnosis=="control",]

# sample plot function - taken from steve expression pipeline

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
  cluster3 <- dendrapply(cluster2,colLab,cluster1order)
  ## PLOTS
  ## cluster IAC
  par(mfrow=c(2,2))
  par(mar=c(5,6,4,2))
  plot(cluster3,nodePar=list(lab.cex=1,pch=NA),
       main=paste("Mean ISA = ",signif(mean(A.IAC[upper.tri(A.IAC)]),3),sep=""),
       xlab="",ylab="1 - ISA",sub="",cex.main=1.8,cex.lab=1.4)
  mtext(paste("distance: 1 - ISA ",sep=""),cex=0.8,line=0.2)
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
  #remove samples with ZK below threshold
  dataset_QC<-dataset[!(rownames(dataset)%in%dataset_removal_1),]
  #sample network plot
  sampleNetwork_plot(dataset_QC)
  #create empty count list to record samples removed
  count<-dataset_removal_1
  # reiterate above till no samples fall below threshold
  while (length(dataset_removal_1)>0) {
    # remove bad samples - 1st iteration removes none
    dataset_QC<-dataset_QC[!(rownames(dataset_QC)%in%dataset_removal_1),]
    #identify sample below Z.K threshold
    dataset_removal_1<-names_of_outliers(dataset_QC, threshold)
    #record samples removed
    count<-c(count, dataset_removal_1)
  }
  #final network plot
  sampleNetwork_plot(dataset_QC)
  # print to screen number of samples removed
  cat("\n")
  print(c("Total number of samples removed...", length(count)))
  # return clean expression set
  return(dataset_QC)
}

# run sample network on entorhinal Cortex - on dataframe without gender

data_case_exprs_good_probes_QC<-run_sample_network_plot(data_exprs_good_probes_diagnosis_sva_case, sample_network_ZK_threshold)
data_control_exprs_good_probes_QC<-run_sample_network_plot(data_exprs_good_probes_diagnosis_sva_control, sample_network_ZK_threshold)

##### PLOT SAMPLE NETWORK ANALYSIS TO PDF #####

setwd(sample_network_dir)

pdf("case_sample_network_analysis.pdf")
run_sample_network_plot(data_exprs_good_probes_diagnosis_sva_case, sample_network_ZK_threshold)
dev.off()

pdf("control_sample_network_analysis.pdf")
run_sample_network_plot(data_exprs_good_probes_diagnosis_sva_control, sample_network_ZK_threshold)
dev.off()

##### CREATE QC'd DATASET #####

# extract sample ID's from QC'd sample network file

# check colnames same in all dataframes - should be FALSE
any(colnames(data_case_exprs_good_probes_QC)==colnames(data_control_exprs_good_probes_QC))==F

# cbind all dataframes

expression_data_QCd<-rbind(data_case_exprs_good_probes_QC,
                           data_control_exprs_good_probes_QC)

dim(expression_data_QCd)

##### PCA PLOT 4 #####

setwd(pca_dir)

pca_data(as.data.frame(t(expression_data_QCd[2:dim(expression_data_QCd)[2]])), 'bottomright')

pdf("4.PCA_plot_after_sample_removal.pdf")
pca_data(as.data.frame(t(expression_data_QCd[2:dim(expression_data_QCd)[2]])), 'bottomright')
dev.off()

##### CONVERT PROBE ID TO ENTREZ ID #####

# Get the probe identifiers that are mapped to an ENTREZ Gene ID using hgu133a.db
mapped_probes <- mappedkeys(eval(parse(text = paste(expression_chip, "ENTREZID", sep=""))))

# Convert to a list
probe_entrez_mapping <- as.data.frame(eval(parse(text = paste(expression_chip, "ENTREZID", sep="")))[mapped_probes])
# arrange order of column by entrezgene probe_id
probe_entrez_mapping<-probe_entrez_mapping[c(2,1)]
colnames(probe_entrez_mapping)[1]<-"entrezgene"

head(probe_entrez_mapping)
dim(probe_entrez_mapping)

#check any duplicated probe IDs
anyDuplicated(probe_entrez_mapping$probe_id)

#check any dupliacted entrezgene IDs
anyDuplicated(probe_entrez_mapping$entrezgene)

# create convert_probe_id_to_entrez_id function 

convert_probe_id_to_entrez_id <- function(expression_dataset, probe_mapping_file){
  # transform dataset # - removed this step
  expression_dataset_t<-as.data.frame(expression_dataset)
  # keep only probes which appear in probe_mapping_file
  data_frame_in_probe_mapper<-expression_dataset_t[colnames(expression_dataset_t)%in%probe_mapping_file$probe_id]
  # match probe id in data_frame_in_probe_mapper to that in probe_mapping_file and convert to entrez id
  colnames(data_frame_in_probe_mapper)<-probe_mapping_file$entrezgene[match(colnames(data_frame_in_probe_mapper), probe_mapping_file$probe_id)]
  return(data_frame_in_probe_mapper)
}

# using probe_entrez_mapping file

expression_data_QCd_entrez_id<-convert_probe_id_to_entrez_id(expression_data_QCd[2:dim(expression_data_QCd)[2]], probe_entrez_mapping)
dim(expression_data_QCd)
dim(expression_data_QCd_entrez_id)
length(which(duplicated(colnames(expression_data_QCd_entrez_id))))

head(expression_data_QCd_entrez_id[100:110])

##### COLLAPSE MULTIPPLE ENTREZ ID BY SELECTING ONE WITH HIGHEST AVERAGE EXPRESSION ACROSS SAMPLES ######

## create function

select_duplicate_probe_by_top_expr <- function(x) {
  # transpose data frame - keep as dataframe
  x_t<-as.data.frame(t(x))
  # calculate mean expression per probe across samples - create new column - probe mean column
  x_t$probe_mean_expression<-rowMeans(x_t)
  #copy rownames (probe id) to column and truncate
  x_t$trunc_entrez_id<-trunc(as.numeric(as.character(rownames(x_t))))
  # order data frame by truncated probe id and then expression level
  x_t<-x_t[order(x_t$trunc_entrez_id, -x_t$probe_mean_expression), ]
  # remove all duplicate probe id - keep one with highest mean expression
  x_t_unique<-x_t[!duplicated(x_t$trunc_entrez_id),]
  #unique entrez column back to row name
  rownames(x_t_unique)<-x_t_unique$trunc_entrez_id
  #remove unwanted column
  x_t_unique$trunc_entrez_id<-NULL
  #remove unwanted column
  x_t_unique$probe_mean_expression<-NULL
  #transpose dataframe back
  x_unique<-as.data.frame(t(x_t_unique))
  return(x_unique)
}

## apply function to dataframes - check number of probes - check duplicates

expression_data_QCd_entrez_id_unique<-select_duplicate_probe_by_top_expr(expression_data_QCd_entrez_id)
dim(expression_data_QCd_entrez_id_unique)
length(which(duplicated(colnames(expression_data_QCd_entrez_id_unique))))

head(expression_data_QCd_entrez_id_unique[1:5])

##### ATTACH DIAGNOSIS AND data REGION #####

# create function to merge multiple dataframes

MyMerge       <- function(x, y){
  df            <- merge(x, y, by= "row.names", all.x= F, all.y= F)
  rownames(df)  <- df$Row.names
  df$Row.names  <- NULL
  return(df)
}

# create phenotype infor to attach - diagnosis + gender + Age + Ethnicity + Tissue
phenotype_to_attach<-Reduce(MyMerge, list(Diagnosis, gender_comparison, Age, Ethnicity, Tissue))

dim(phenotype_to_attach)
head(phenotype_to_attach)

# attach pheno to exprs table

expression_data_QCd_entrez_id_unique_pheno<-merge(phenotype_to_attach, expression_data_QCd_entrez_id_unique, by="row.names")

rownames(expression_data_QCd_entrez_id_unique_pheno)<-expression_data_QCd_entrez_id_unique_pheno$Row.names

expression_data_QCd_entrez_id_unique_pheno$Row.names<-NULL

head(expression_data_QCd_entrez_id_unique_pheno)[1:10]

# rows should be same in exprs table - should be TRUE
dim(expression_data_QCd_entrez_id_unique)[1]==dim(expression_data_QCd_entrez_id_unique_pheno)[1]

##### DE #####

#function for DE analysis
run_diff_exprs_analysis<-function(dataset, x){
  #sort dataset by Diagnosis column - want Case samples 1st
  dataset_sorted<-dataset[order(dataset$Diagnosis),]
  #split dataset into expres and diagnosis
  dataset_exprs<-dataset_sorted[7:dim(dataset_sorted)[2]]
  dataset_pheno<-dataset_sorted[,1]
  #replace sample name to numbers
  rownames(dataset_exprs)<-c(1:dim(dataset)[1])
  # setup experimental desgin
  design <- model.matrix(~0 + dataset_pheno)
  # change colnames to Case and Control
  colnames(design)<-c("Case", "Control")
  # transpose dataset, convert to numeric 
  transposed_dataset_exprs<-t(dataset_exprs)
  #run diff expression
  dataset_exprs_fit <- lmFit(transposed_dataset_exprs, design, method="robust", maxit=x)
  dataset_exprs_contrast_matrix<- makeContrasts(Case-Control, levels=design)
  dataset_exprs_contrast_fit <-contrasts.fit(dataset_exprs_fit, dataset_exprs_contrast_matrix)
  dataset_exprs_contrast_ebayes <- eBayes(dataset_exprs_contrast_fit, robust=T)
  dataset_exprs_top_genes <- topTable(dataset_exprs_contrast_ebayes, number=(dim(transposed_dataset_exprs)[1]), coef=1, adjust.method="fdr", confint=TRUE) 
  return(dataset_exprs_top_genes)
}

#apply function
DE_results<-run_diff_exprs_analysis(expression_data_QCd_entrez_id_unique_pheno, 1000)

#sig_DE_results
sig_DE_results<-subset(DE_results, adj.P.Val<=0.05)
dim(sig_DE_results)

#write out DE results
setwd(DE_dir)
write.table(DE_results, file="Full_Differential_Expression_Results.txt")

##### SUMMARY #####

print(c("dataset:", dataset), quote=F)
print(c("Disease:", disease), quote=F)
print(c("Microarray Platform:", Microarray_platform), quote=F)
print(c("Expression Chip:", expression_chip), quote=F)
print(c("Data Format:", Data_format), quote=F)
print(c("Tissue:", as.character(Tissue[1,1])), quote=F)
print(c("Case Number:", length(case_ID)), quote=F)
print(c("Control Number:", length(control_ID)), quote=F)
print(c("Probe Detection Threshold:", Probe_Detection_Threshold), quote=F)
print(c("Number of SV:", number_of_SV), quote=F)
print(c("Gender-Missmatch:", dim(Gender_Missmatch)[1]), quote=F)
print(c("Samples Removed:", dim(expression_data_normalised_as_data_frame)[2]-dim(expression_data_QCd_entrez_id_unique_pheno)[1]), quote=F)
print(c("Final Case numbers:", length(expression_data_QCd_entrez_id_unique_pheno[expression_data_QCd_entrez_id_unique_pheno$Diagnosis=="case",1])), quote=F)
print(c("Final Control Numbers:", length(expression_data_QCd_entrez_id_unique_pheno[expression_data_QCd_entrez_id_unique_pheno$Diagnosis=="control",1])), quote=F)
print(c("Initial Probe Numbers:", dim(expression_data_normalised_as_data_frame)[1]), quote=F)
print(c("Final Probe Numbers:", dim(expression_data_QCd_entrez_id_unique_pheno)[2]), quote=F)
print(c("Number of sig DEG:", dim(sig_DE_results)[1]), quote=F)

#write out summary report
setwd(work_dir)
cat("Processing Date:", strftime(Sys.Date(),"%Y-%m-%d"), "\n", file="Summary_Report.txt")
cat("Dataset:", dataset, "\n", file="Summary_Report.txt", append=T)
cat("Disease:", disease, "\n", file="Summary_Report.txt", append=T)
cat("Microarray Platform:", Microarray_platform, "\n", file="Summary_Report.txt", append=T)
cat("Expression Chip:", expression_chip, "\n", file="Summary_Report.txt", append=T)
cat("Data Format:", Data_format, "\n", file="Summary_Report.txt", append=T)
cat("Tissue:", as.character(Tissue[1,1]), "\n", file="Summary_Report.txt", append=T)
cat("Case Numbers before QC:", length(case_ID), "\n", file="Summary_Report.txt", append=T)
cat("Control Numbers before QC:", length(control_ID), "\n", file="Summary_Report.txt", append=T)
cat("Probe Detection Threshold:", Probe_Detection_Threshold, "\n", file="Summary_Report.txt", append=T)
cat("Number of SV:", number_of_SV, "\n", file="Summary_Report.txt", append=T)
cat("Gender-Missmatch:", dim(Gender_Missmatch)[1], "\n", file="Summary_Report.txt", append=T)
cat("Samples Removed:", dim(expression_data_normalised_as_data_frame)[2]-dim(expression_data_QCd_entrez_id_unique_pheno)[1], "\n", file="Summary_Report.txt", append=T)
cat("Case numbers after QC:", length(expression_data_QCd_entrez_id_unique_pheno[expression_data_QCd_entrez_id_unique_pheno$Diagnosis=="case",1]), "\n", file="Summary_Report.txt", append=T)
cat("Control Numbers after QC:", length(expression_data_QCd_entrez_id_unique_pheno[expression_data_QCd_entrez_id_unique_pheno$Diagnosis=="control",1]), "\n", file="Summary_Report.txt", append=T)
cat("Probe Numbers before QC:", dim(expression_data_normalised_as_data_frame)[1], "\n", file="Summary_Report.txt", append=T)
cat("Probe Numbers after QC:", dim(expression_data_QCd_entrez_id_unique_pheno)[2], "\n", file="Summary_Report.txt", append=T)
cat("Number of sig DEG:", dim(sig_DE_results)[1], "\n", file="Summary_Report.txt", append=T) 

##### SAVE EXPRESSION DATAFRAME #####

setwd(clean_data_dir)

write.table(expression_data_QCd_entrez_id_unique_pheno, file=paste(dataset, "processing_data.txt", sep="_"), sep="\t")

##### SAVE IMAGE #####

setwd(work_dir)

save.image(file=paste(dataset, "processing_data.Rdata", sep="_"))

