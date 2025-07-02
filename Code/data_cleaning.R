###############################################################################
###############################################################################
###############################################################################
#Creator: Kylee Rahm
#Advisor: Matthew Bailey
#Last Updated: April 4, 2025
#MCP_Microglial_Classifier for Cell-type classifier for HMC3 cells
#Published: doi::
###############################################################################
###############################################################################

#==============================================================================
##---------------LIBRARIES-----------------------------------------------------
#==============================================================================
library(tidyverse)
library(dplyr)
library(utils)
library(multiclassPairs)

#NOTE: conda install for M3 or M4 Macs doesn't work well. Best to use a linux 
#system

#===============================================================================
##------------------LOAD IN DATA------------------------------------------------
#===============================================================================

#data for primary cells / classifier 
rawDataHansen <- read_csv("Data/hansenData_1250.csv")
dfGalFinal <- read_csv("Data/galatroDataFormatted.csv")
rawDataPrimaryBrain <- read_csv("Data/primary_brain_data.csv")
rawDataGoss <- read_csv("Data/tpm_cleaned_gosselin.csv")
rawDataVerd <- read.delim("Data/verdera_featureCounts_results.txt", comment.char="#")

#data for cell lines / classifier
rawDataChandra <- read_tsv("Data/chandra_counts_TPM.tsv")
rawGupta <- read_tsv("Data/gupta_norm_counts_TPM_GRCh38.p13_NCBI.tsv")
rawWei <- read_tsv("Data/wei_norm_counts_TPM_GRCh38.p13_NCBI.tsv")
rawYu <- read_tsv("Data/yu_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")
rawKaya <- read_tsv("Data/kaya_norm_counts_TPM_GRCh38.p13_NCBI.tsv")
rawGoyal <- read_tsv("Data/goyal_norm_counts_TPM_GRCh38.p13_NCBI.tsv")
rawZimmer <- read_tsv("Data/zimmer_norm_counts_TPM_GRCh38.p13_NCBI.tsv")
rawDataLu <- read.delim("Data/counts_lu.txt", header = TRUE, sep = "\t")
rawDataSunMOLM <- read_csv("Data/sun_MOLM13_counts.tpm.csv")
rawDataSunTHP <- read_csv("Data/sun_THP1_counts.tpm.csv")
rawBarbMOLM <- read.delim("Data/barbosa_MOLM13_TPM.txt", header = TRUE, sep = "\t")
rawBarbU937 <- read.delim("Data/barbosa_U937_TPM.txt", header = TRUE, sep = "\t")
rawSubo <- read_tsv("Data/subo_norm_counts_TPM.tsv")

#validation / HMC3 datasets 
rawDataAbud <- read_tsv("Data/abud_TPMcounts.tsv")
rawDataADAB <- read_tsv("Data/GSE181153_ADAB_geneCounts.tsv")
rawDataRow <-read_csv("Data/dataRowitch.csv")
rawDataSagar <- read_csv("Data/dataSagar.csv")
rawDataChai <- read_tsv("Data/hmc3_chai.tsv")
rawDataArman <- read_tsv("Data/arman_summary_star_genes_tpm.stranded.annotated.tsv")
human_geneannotation <- read_tsv("Data/abud_human_geneannotation.tsv")

#===============================================================================
#                  DATASETS THAT WILL BE USED IN THE CLASSIFIER 
#                       SPECIFICALLY -- PRIMARY CELLS
#===============================================================================

#-------------------------------------------------------------------------------
#                         CREATE HANSEN DF
#-------------------------------------------------------------------------------

#load in data - Hansen renamed for first author Srinivasan
rawDfHansen <- rawDataHansen
dimsHansen <- dim(rawDfHansen)
cellTypeHansen <- rawDfHansen%>%  #takes cell type names from df
  slice(dimsHansen[1]) %>%
  unlist(., use.names=FALSE)
cellTypeHansen = cellTypeHansen[2:dimsHansen[2]] #removes "celltype" name
dfHansen <- rawDfHansen[-(dimsHansen[1]),] #removes cell type row from raw df

#replace Gene names containing - or . to underscores
dfHansen$Gene <- gsub('-', '_', dfHansen$Gene)
dfHansen$Gene <- gsub('\\.','_', dfHansen$Gene)

#we label myeloid as microglia cells
len = length(cellTypeHansen)
for (i in 1:len) {
  if (str_detect(cellTypeHansen[i], regex("myeloid", ignore_case=TRUE))) {
    cellTypeHansen[i] = "Microglia"
  }
  else if (str_detect(cellTypeHansen[i], regex("astrocyte", ignore_case=TRUE))) {
    cellTypeHansen[i] = "Astrocyte"
  }
  else if (str_detect(cellTypeHansen[i], regex("endothelial", ignore_case=TRUE))) {
    cellTypeHansen[i] = "Endothelial"
  }
  else if (str_detect(cellTypeHansen[i], regex("neuron", ignore_case=TRUE))) {
    cellTypeHansen[i] = "Neuron"
  }
}
#create class labels for HANSEN data
classLabelsHansen = cellTypeHansen

#create platform labels for HANSEN data
lenH <- length(cellTypeHansen)
platformHansen <- rep("Srinivasan Data", lenH) #first author name

##labels: classLabelsHansen and platformHansen 
## final df is dfHansen


#-------------------------------------------------------------------------------
#                     CREATE GALATRO DF
#-------------------------------------------------------------------------------

##create class and platform labels 
desiredLength <- length(dfGalFinal) - 1 #because we don't want to include the Gene column
classLabelsGal <- rep("Microglia", desiredLength)
platformGal <- rep("Galatro Data", desiredLength)

#replace Gene names containing - or . to underscores
dfGalFinal$Gene <- gsub('-', '_', dfGalFinal$Gene)
dfGalFinal$Gene <- gsub('\\.','_', dfGalFinal$Gene)

#Galatro data was already previously formatted and had ENSEMBL genes changed to HUGO in previous code

## labels: classLabelsGal and platformGal
#final dataframe is dfGalFinal


#-------------------------------------------------------------------------------
#               CREATE PRIMARY BRAIN DF / ZHANG 
#-------------------------------------------------------------------------------

#load in data
rawDfPrimaryBrain <- rawDataPrimaryBrain
dfPrimaryBrain <- rawDfPrimaryBrain

#replace Gene names containing - or . to underscores
dfPrimaryBrain$Gene <- gsub('-', '_', dfPrimaryBrain$Gene)
dfPrimaryBrain$Gene <- gsub('\\.','_', dfPrimaryBrain$Gene)

#take out diseased / unwanted samples 
unwantedSamples <- c("GBM", "sclerotic", "Oligodendrocyte", "cortex")
keepCols <- grep(paste(unwantedSamples, collapse = "|"), colnames(dfPrimaryBrain), invert=TRUE)

#subset the dataframe
healthyDfPrimaryBrain <- dfPrimaryBrain[,keepCols]

dfPrimaryBrain <- healthyDfPrimaryBrain

#pull colnames as a vector, change them to use as class labels 
columnNamesPrimaryBrain = colnames(dfPrimaryBrain)

len = length(columnNamesPrimaryBrain)
for (i in 1:len) {
  if (str_detect(columnNamesPrimaryBrain[i], regex("astro", ignore_case=TRUE))) {
    columnNamesPrimaryBrain[i] = "Astrocyte"
  }
  else if (str_detect(columnNamesPrimaryBrain[i], regex("oligo", ignore_case=TRUE))) {
    columnNamesPrimaryBrain[i] = "Oligodendrocyte"
  }
  else if (str_detect(columnNamesPrimaryBrain[i], regex("microg", ignore_case=TRUE))) {
    columnNamesPrimaryBrain[i] = "Microglia"
  }
  else if (str_detect(columnNamesPrimaryBrain[i], regex("cortex", ignore_case=TRUE))) {
    columnNamesPrimaryBrain[i] = "Whole.cortex"
  }
  else if (str_detect(columnNamesPrimaryBrain[i], regex("neuron", ignore_case=TRUE))) {
    columnNamesPrimaryBrain[i] = "Neuron"
  }
}

#create class labels for PrimaryBrain data
columnNamesPrimaryBrain = columnNamesPrimaryBrain[2:len] #excludes gene column
classLabelsPrimaryBrain = columnNamesPrimaryBrain

#create platform labels for PrimaryBrain data
lenB <- length(classLabelsPrimaryBrain) #get times of repeat for platform
platformBrain <- rep("Zhang Data", lenB)

##labels: classLabelsPrimaryBrain and platformBrain
## df is dfPrimaryBrain


#-------------------------------------------------------------------------------
#                 CREATE GOSSELIN DF
#-------------------------------------------------------------------------------
#load in data
rawDfGoss <- rawDataGoss
dfGoss <- rawDfGoss[,-grep("Lysate", colnames(rawDfGoss))] #remove lysate samples, verified it works 
#below code verifies that lysate samples are gone
#print("DF GOSS AFTER REMOVING LYSATE COLUMNS:")
#head(dfGoss)

#replace Gene names containing - or . to underscores
dfGoss$Gene <- gsub('-', '_', dfGoss$Gene)
dfGoss$Gene <- gsub('\\.','_', dfGoss$Gene)
#to verify the gene names changed correctly
#print("DF GOSS GENE COLUMN AFTER REPLACING DASHES AND PERIODS")
#print(dfGoss$Gene)

#grab and change column names for class labels 
lastN <- ncol(dfGoss) #number of columns
precolumnNamesGoss = colnames(dfGoss)
precolumnNamesGoss = precolumnNamesGoss[2:lastN] #takes cell col names, excludes "Gene" col

columnNamesGoss = precolumnNamesGoss

#we label microglia and monocytes 
len = length(columnNamesGoss)
for (i in 1:len) {
  if (str_detect(columnNamesGoss[i], regex("Microglia", ignore_case=TRUE))) {
    columnNamesGoss[i] = "Microglia"
  }
  else if (str_detect(columnNamesGoss[i], regex("monocyte", ignore_case=TRUE))) {
    columnNamesGoss[i] = "Monocyte"
  }
}

#create class labels for Goss data
classLabelsGoss = columnNamesGoss

#create platform labels for ADAB data
lenC <- length(classLabelsGoss)
platformGoss <- rep("Gosselin Data", lenC)

##labels: classLabelsGoss and platformGoss 
##df is dfGoss


#-------------------------------------------------------------------------------
#                         CREATE VERDERA DF
#-------------------------------------------------------------------------------

#load in data
rawDfVerd <- rawDataVerd

#replace Gene names containing - or . to underscores
rawDfVerd$Geneid <- gsub('-', '_', rawDfVerd$Geneid)
rawDfVerd$Geneid <- gsub('\\.','_', rawDfVerd$Geneid)

#rename Geneid column to uniform Gene name
rawDfVerd <- rename(rawDfVerd, Gene = Geneid)

#get only untreated samples
dfVerd <- select(rawDfVerd, matches("Gene|UT"))

#prepare to create class labels
colsVerd <- colnames(dfVerd)[-1]
cellTypeVerd <- c()

#create labels
lenV = length(colsVerd)
for (i in 1:lenV) {
  if (str_detect(colsVerd[i], "Astro")) {
    cellTypeVerd[i] = "Astrocyte"
  }
  else if (str_detect(colsVerd[i], "Neu")) {
    cellTypeVerd[i] = "Neuron"
  }
}
#create class labels for HANSEN data
classLabelsVerd = cellTypeVerd

#create platform labels for HANSEN data
lenV <- length(cellTypeVerd)
platformVerd <- rep("Verdera Data", lenV) #first author name

##labels: classLabelsVerd and platformVerd
## final df is dfVerd


#=========================================================================================
#          CREATING THE PRIMARY CELL CLASSIFIER DATAFRAME 
#=========================================================================================

#ORDER

# ALL VARIABLES: order top to bottom 

##labels: classLabelsHansen and platformHansen 
## final df is dfHansen

## labels: classLabelsGal and platformGal
#final dataframe is dfGalFinal

##labels: classLabelsPrimaryBrain and platformBrain
## df is dfPrimaryBrain

##labels: classLabelsGoss and platformGoss 
##df is dfGoss

##labels: classLabelsVerd and platformVerd
## final df is dfVerd


allClassLabels <- c(classLabelsHansen, classLabelsGal, classLabelsPrimaryBrain, classLabelsGoss,
                    classLabelsVerd)

allPlatformLabels <- c(platformHansen, platformGal, platformBrain, platformGoss,
                       platformVerd)

#merge all dataframes into one huge one
distinctHansen <- distinct(dfHansen, Gene, .keep_all=TRUE)
distinctGal <- distinct(dfGalFinal, Gene, .keep_all=TRUE)
distinctPrimaryBrain <- distinct(dfPrimaryBrain, Gene, .keep_all=TRUE)
distinctGoss <- distinct(dfGoss, Gene, .keep_all=TRUE)
distinctVerd <- distinct(dfVerd, Gene, .keep_all=TRUE)

distinct_df_list <- list(distinctHansen, distinctGal, distinctPrimaryBrain, distinctGoss,
                         distinctVerd)
full_distinctDf <- reduce(distinct_df_list, full_join, by="Gene")
dimsdistinct <- dim(full_distinctDf)


#AGAIN SEE BELOW >> I'm going to make somee changes to the filal dataset for df because I don't think it was done in the safestt way. See below for the updated changes to the combined dataset. 


print(paste("The full dataframe has", dimsdistinct[2], "columns. The labels should have",
            dimsdistinct[2]-1, "length. The lengths of the class labels and platforms labels are",
            length(allClassLabels),"and", length(allPlatformLabels)))

#writing all data to files for easy use 

write.table(allClassLabels, "Processed_data/classlabels_primary.csv",
            row.names=FALSE,col.names=FALSE,quote=FALSE, sep=",")

write.table(allPlatformLabels, "Processed_data/platformlabels_primary.csv",
            row.names=FALSE,col.names=FALSE,quote=FALSE, sep=",")

#dataframes
write_csv(full_distinctDf, "Processed_data/cleandata_combined_primary.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfHansen, "Processed_data/cleandata_Hansen.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfGalFinal, "Processed_data/cleandata_Gal.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfPrimaryBrain, "Processed_data/cleandata_Zhang.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfGoss, "Processed_data/cleandata_Gosselin.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfVerd, "Processed_data/cleandata_Verdera.csv", append=FALSE,
          col_names=TRUE,quote="none")


#=========================================================================================
#                  DATASETS THAT WILL BE USED IN THE CLASSIFIER 
#                       SPECIFICALLY -- CELL LINES
#=========================================================================================

#-------------------------------------------------------------------------------
#                       CREATE CHANDRA DF - U87 samples
#-------------------------------------------------------------------------------

#load in data
rawChandra <- rawDataChandra

#only keep wanted samples for training (control HMC3 samples)
control_samples_chandra <- c("GSM8252622", "GSM8252623")
controlonlyDf_chandra <- dplyr::select(rawChandra, c("GeneID", any_of(control_samples_chandra)))

#get HUGO gene names
genes_onlyDf <- select(human_geneannotation, GeneID:Symbol)
gene_names <- pull(genes_onlyDf, Symbol)

extragenesDf_chandra <- merge(genes_onlyDf, controlonlyDf_chandra, by="GeneID", sort=FALSE)

#check that gene names are correct 
sum(gene_names == extragenesDf_chandra$Symbol) == dim(extragenesDf_chandra)[1]

#take out extra column
dfChandra <- select(extragenesDf_chandra, -GeneID) %>%
  rename(Gene=Symbol)

#replace Gene names containing - or . to underscores
dfChandra$Gene <- gsub('-', '_', dfChandra$Gene)
dfChandra$Gene <- gsub('\\.','_', dfChandra$Gene)

#create class labels for data
classLabelsChandra = c("U87", "U87")

#create platform labels for data
platformChandra <- rep("Bhaskara Data", 2)

##labels: classLabelsChandra and platformChandra 
#final is dfChandra


#-------------------------------------------------------------------------------
#                       CREATE GUPTA DF - U87 samples
#-------------------------------------------------------------------------------

#DATA = rawGupta

#only keep wanted samples for training (control, unaltered U87 samples)
control_samples_gupta <- c("GSM5640594", "GSM5640595", "GSM5640596", "GSM5640597")
controlonlyDf_gupta <- dplyr::select(rawGupta, c("GeneID", any_of(control_samples_gupta)))

#get HUGO gene names
genes_onlyDf <- select(human_geneannotation, GeneID:Symbol)
gene_names <- pull(genes_onlyDf, Symbol)

extragenesDf_gupta <- merge(genes_onlyDf, controlonlyDf_gupta, by="GeneID", sort=FALSE)

#check that gene names are correct 
sum(gene_names == extragenesDf_gupta$Symbol) == dim(extragenesDf_gupta)[1]

#take out extra column
dfGupta <- select(extragenesDf_gupta, -GeneID) %>%
  rename(Gene=Symbol)

#replace Gene names containing - or . to underscores
dfGupta$Gene <- gsub('-', '_', dfGupta$Gene)
dfGupta$Gene <- gsub('\\.','_', dfGupta$Gene)

#create class labels for data
classLabelsGupta = rep("U87", 4)

#create platform labels for data
platformGupta <- rep("Gupta Data", 4)

##labels: classLabelsGupta and platformGupta
#final is dfGupta



#-------------------------------------------------------------------------------
#                       CREATE WEIWEI DF - U87 samples
#-------------------------------------------------------------------------------

#DATA = rawWei

#only keep wanted samples for training (control U87 samples)
control_samples_wei <- c("GSM7510194", "GSM7510195", "GSM7510196", "GSM7510197",
                         "GSM7510201", "GSM7510202")
controlonlyDf_wei <- dplyr::select(rawWei, c("GeneID", any_of(control_samples_wei)))

#get HUGO gene names
genes_onlyDf <- select(human_geneannotation, GeneID:Symbol)
gene_names <- pull(genes_onlyDf, Symbol)

extragenesDf_wei <- merge(genes_onlyDf, controlonlyDf_wei, by="GeneID", sort=FALSE)

#check that gene names are correct 
sum(gene_names == extragenesDf_wei$Symbol) == dim(extragenesDf_wei)[1]

#take out extra column
dfWei <- select(extragenesDf_wei, -GeneID) %>%
  rename(Gene=Symbol)

#replace Gene names containing - or . to underscores
dfWei$Gene <- gsub('-', '_', dfWei$Gene)
dfWei$Gene <- gsub('\\.','_', dfWei$Gene)

#create class labels for data
classLabelsWei = rep("U87", 6)

#create platform labels for data
platformWei <- rep("Weiwei Data", 6)

##labels: classLabelsWei and platformWei 
#final is dfWei


#-------------------------------------------------------------------------------
#                       CREATE YU DF - U87 samples
#-------------------------------------------------------------------------------

#DATA = rawYu

#only keep wanted samples for training (control U87 samples)
control_samples_yu <- c("GSM7846269", "GSM7846270", "GSM7846271", "GSM7846275",
                         "GSM7846276", "GSM7846277") #69 sample is not there
controlonlyDf_yu <- dplyr::select(rawYu, c("GeneID", any_of(control_samples_yu)))

#get HUGO gene names
genes_onlyDf <- select(human_geneannotation, GeneID:Symbol)
gene_names <- pull(genes_onlyDf, Symbol)

extragenesDf_yu <- merge(genes_onlyDf, controlonlyDf_yu, by="GeneID", sort=FALSE)

#check that gene names are correct 
sum(gene_names == extragenesDf_yu$Symbol) == dim(extragenesDf_yu)[1]

#take out extra column
dfYu <- select(extragenesDf_yu, -GeneID) %>%
  rename(Gene=Symbol)

#replace Gene names containing - or . to underscores
dfYu$Gene <- gsub('-', '_', dfYu$Gene)
dfYu$Gene <- gsub('\\.','_', dfYu$Gene)

#create class labels for data
classLabelsYu = rep("U87", 5)

#create platform labels for data
platformYu <- rep("Yu Data", 5)

##labels: classLabelsYu and platformYu 
#final is dfYu



#-------------------------------------------------------------------------------
#                       CREATE KAYA DF - A172 samples
#-------------------------------------------------------------------------------
#DATA = rawKaya

#only keep wanted samples for training (control A172 samples)
control_samples_kaya <- c("GSM3937531", "GSM3937532")
controlonlyDf_kaya <- dplyr::select(rawKaya, c("GeneID", any_of(control_samples_kaya)))

#get HUGO gene names
genes_onlyDf <- select(human_geneannotation, GeneID:Symbol)
gene_names <- pull(genes_onlyDf, Symbol)

extragenesDf_kaya <- merge(genes_onlyDf, controlonlyDf_kaya, by="GeneID", sort=FALSE)

#check that gene names are correct 
sum(gene_names == extragenesDf_kaya$Symbol) == dim(extragenesDf_kaya)[1]

#take out extra column
dfKaya <- select(extragenesDf_kaya, -GeneID) %>%
  rename(Gene=Symbol)

#replace Gene names containing - or . to underscores
dfKaya$Gene <- gsub('-', '_', dfKaya$Gene)
dfKaya$Gene <- gsub('\\.','_', dfKaya$Gene)

#create class labels for data
classLabelsKaya = c("A172", "A172")

#create platform labels for data
platformKaya <- rep("Kaya Data", 2)

##labels: classLabelsKaya and platformKaya
#final is dfKaya


#-------------------------------------------------------------------------------
#                       CREATE GOYAL DF - A172 samples
#-------------------------------------------------------------------------------
#DATA = rawGoyal

#only keep wanted samples for training (control A172 samples)
control_samples_goyal <- c("GSM6395107", "GSM6395109", "GSM6395111")
controlonlyDf_goyal <- dplyr::select(rawGoyal, c("GeneID", any_of(control_samples_goyal)))

#get HUGO gene names
genes_onlyDf <- select(human_geneannotation, GeneID:Symbol)
gene_names <- pull(genes_onlyDf, Symbol)

extragenesDf_goyal <- merge(genes_onlyDf, controlonlyDf_goyal, by="GeneID", sort=FALSE)

#check that gene names are correct 
sum(gene_names == extragenesDf_goyal$Symbol) == dim(extragenesDf_goyal)[1]

#take out extra column
dfGoyal <- select(extragenesDf_goyal, -GeneID) %>%
  rename(Gene=Symbol)

#replace Gene names containing - or . to underscores
dfGoyal$Gene <- gsub('-', '_', dfGoyal$Gene)
dfGoyal$Gene <- gsub('\\.','_', dfGoyal$Gene)

#create class labels for data
classLabelsGoyal = c("A172", "A172", "A172")

#create platform labels for data
platformGoyal <- rep("Goyal Data", 3)

##labels: classLabelsGoyal and platformGoyal
#final is dfGoyal


#-------------------------------------------------------------------------------
#                       CREATE ZIMMER DF - A172 samples
#-------------------------------------------------------------------------------
#DATA = rawZimmer

#only keep wanted samples for training (control A172 samples)
# GSM6012714 is supposed to be a control, but it's not in the count matrix data 
control_samples_zimmer <- c("GSM6012712", "GSM6012713", "GSM6012715")
controlonlyDf_zimmer <- dplyr::select(rawZimmer, c("GeneID", any_of(control_samples_zimmer)))

#get HUGO gene names
genes_onlyDf <- select(human_geneannotation, GeneID:Symbol)
gene_names <- pull(genes_onlyDf, Symbol)

extragenesDf_zimmer <- merge(genes_onlyDf, controlonlyDf_zimmer, by="GeneID", sort=FALSE)

#check that gene names are correct 
sum(gene_names == extragenesDf_zimmer$Symbol) == dim(extragenesDf_zimmer)[1]

#take out extra column
dfZimmer <- select(extragenesDf_zimmer, -GeneID) %>%
  rename(Gene=Symbol)

#replace Gene names containing - or . to underscores
dfZimmer$Gene <- gsub('-', '_', dfZimmer$Gene)
dfZimmer$Gene <- gsub('\\.','_', dfZimmer$Gene)

#create class labels for data
classLabelsZimmer = c("A172", "A172", "A172")

#create platform labels for data
platformZimmer <- rep("Zimmermannova Data", 3)

##labels: classLabelsZimmer and platformZimmer
#final is dfZimmer


#-------------------------------------------------------------------------------
#                       CREATE LU DF - MOLM13 samples
#-------------------------------------------------------------------------------
#load in data
rawLu <- rawDataLu

#only keep wanted samples for training (control HMC3 samples)
control_samples_lu <- str_subset(colnames(rawLu), "^C\\d")
controlonlyDf_lu <- dplyr::select(rawLu, c("Geneid", any_of(control_samples_lu)))

#correct column name 
dfLu <- rename(controlonlyDf_lu, Gene=Geneid)

#replace Gene names containing - or . to underscores
dfLu$Gene <- gsub('-', '_', dfLu$Gene)
dfLu$Gene <- gsub('\\.','_', dfLu$Gene)

#create class labels for data
classLabelsLu = c("MOLM13", "MOLM13", "MOLM13")

#create platform labels for data
platformLu <- rep("Chen Data", 3)

##labels: classLabelsLu and platformLu 
#final is dfLu


#-------------------------------------------------------------------------------
#                   CREATE SUN DF - MOLM13 & THP1 samples
#-------------------------------------------------------------------------------
#load in data
rawSunMOLM <- rawDataSunMOLM

#replace Gene names containing - or . to underscores
rawSunMOLM$Gene <- gsub('-', '_', rawSunMOLM$Gene)
rawSunMOLM$Gene <- gsub('\\.','_', rawSunMOLM$Gene)

#create class labels for data
classLabelsSunMOLM = c("MOLM13", "MOLM13")

#create platform labels for data
platformSunMOLM <- rep("Sun Data", 2)

dfSunMOLM <- rawSunMOLM

##labels: classLabelsSunMOLM and platformSunMOLM
#final is dfSunMOLM

#load in data
rawSunTHP <- rawDataSunTHP

#replace Gene names containing - or . to underscores
rawSunTHP$Gene <- gsub('-', '_', rawSunTHP$Gene)
rawSunTHP$Gene <- gsub('\\.','_', rawSunTHP$Gene)

#create class labels for data
classLabelsSunTHP = c("THP1", "THP1")

#create platform labels for data
platformSunTHP <- rep("Sun Data", 2)

dfSunTHP <- rawSunTHP

##labels: classLabelsSunTHP and platformSunTHP
#final is dfSunTHP

#join the Sun data
dfSun <- full_join(dfSunMOLM, dfSunTHP, by = "Gene")
classLabelsSun <- c(classLabelsSunMOLM, classLabelsSunTHP)
platformSun <- c(platformSunMOLM, platformSunTHP)

##labels: classLabelsSun and platformSun
#final is dfSun



#-------------------------------------------------------------------------------
#                   CREATE BARBOSA DF - MOLM13 & U937 samples
#-------------------------------------------------------------------------------
#data is rawBarbMOLM and rawBarbU937

#MOLM13 data
#-------------------------------------------------------------------------------
#only keep wanted samples for training (control MOLM13 samples)
control_samples_barbmolm <- str_subset(colnames(rawBarbMOLM), "NTC")
dfBarbMOLM <- dplyr::select(rawBarbMOLM, c("gene_name", any_of(control_samples_barbmolm)))

#rename Gene column
dfBarbMOLM <- rename(dfBarbMOLM, Gene = gene_name)

#replace Gene names containing - or . to underscores
dfBarbMOLM$Gene <- gsub('-', '_', dfBarbMOLM$Gene)
dfBarbMOLM$Gene <- gsub('\\.','_', dfBarbMOLM$Gene)

#create class labels for data
classLabelsBarbMOLM = c("MOLM13", "MOLM13", "MOLM13")

#create platform labels for data
platformBarbMOLM <- rep("Barbosa Data", 3)

##labels: classLabelsBarbMOLM and platformBarbMOLM 
#final is dfBarbMOLM

#U937 data
#-------------------------------------------------------------------------------
#only keep wanted samples for training (control U937 samples)
control_samples_barbu937 <- str_subset(colnames(rawBarbU937), "NTC")
dfBarbU937 <- dplyr::select(rawBarbU937, c("gene_name", any_of(control_samples_barbu937)))

#rename Gene column
dfBarbU937 <- rename(dfBarbU937, Gene = gene_name)

#replace Gene names containing - or . to underscores
dfBarbU937$Gene <- gsub('-', '_', dfBarbU937$Gene)
dfBarbU937$Gene <- gsub('\\.','_', dfBarbU937$Gene)

#create class labels for data
classLabelsBarbU937 = c("U937", "U937", "U937")

#create platform labels for data
platformBarbU937 <- rep("Barbosa Data", 3)

##labels: classLabelsBarbU937 and platformBarbU937 
#final is dfBarbU937


#-------------------------------------------------------------------------------
#                       CREATE SUBO DF - MOLM13 samples
#-------------------------------------------------------------------------------
#rawSubo

#only keep wanted samples for training (control MOLM13 samples)
# sample GSM7869893 not found in their count matrix 
control_samples_subo <- c("GSM7869894")
controlonlyDf_subo <- dplyr::select(rawSubo, c("GeneID", any_of(control_samples_subo)))

#get HUGO gene names
genes_onlyDf <- select(human_geneannotation, GeneID:Symbol)
gene_names <- pull(genes_onlyDf, Symbol)

extragenesDf_subo <- merge(genes_onlyDf, controlonlyDf_subo, by="GeneID", sort=FALSE)

#check that gene names are correct 
sum(gene_names == extragenesDf_subo$Symbol) == dim(extragenesDf_subo)[1]

#take out extra column
dfSubo <- select(extragenesDf_subo, -GeneID) %>%
  rename(Gene=Symbol)

#replace Gene names containing - or . to underscores
dfSubo$Gene <- gsub('-', '_', dfSubo$Gene)
dfSubo$Gene <- gsub('\\.','_', dfSubo$Gene)

#create class labels for data
classLabelsSubo = c("MOLM13")

#create platform labels for data
platformSubo <- rep("S. Zhang Data", 1)

##labels: classLabelsSubo and platformSubo 
#final is dfSubo


#=========================================================================================
#          CREATING THE CELL LINE CLASSIFIER DATAFRAME 
#=========================================================================================

#ORDER

# ALL VARIABLES: order top to bottom 

##labels: classLabelsChandra and platformChandra 
#final is dfChandra

##labels: classLabelsGupta and platformGupta
#final is dfGupta

##labels: classLabelsWei and platformWei 
#final is dfWei

##labels: classLabelsYu and platformYu 
#final is dfYu

##labels: classLabelsKaya and platformKaya
#final is dfKaya

##labels: classLabelsGoyal and platformGoyal
#final is dfGoyal

##labels: classLabelsZimmer and platformZimmer
#final is dfZimmer

##labels: classLabelsADAB and platformADAB 
##df is dfADAB

##labels: classLabelsLu and platformLu 
#final is dfLu

##labels: classLabelsSun and platformSun
#final is dfSun

##labels: classLabelsBarbMOLM and platformBarbMOLM 
#final is dfBarbMOLM

##labels: classLabelsBarbU937 and platformBarbU937 
#final is dfBarbU937

##labels: classLabelsSubo and platformSubo 
#final is dfSubo


allClassLabels_line <- c(classLabelsChandra, classLabelsGupta, classLabelsWei,
                         classLabelsYu, classLabelsLu, classLabelsSunMOLM,
                         classLabelsBarbMOLM, classLabelsSubo)

allPlatformLabels_line <- c(platformChandra, platformGupta, platformWei,
                            platformYu, platformLu, platformSunMOLM, platformBarbMOLM,
                            platformSubo)

#merge all dataframes into one huge one
distinctChandra <- distinct(dfChandra, Gene, .keep_all=TRUE)
distinctGupta <- distinct(dfGupta, Gene, .keep_all=TRUE)
distinctWei <- distinct(dfWei, Gene, .keep_all=TRUE)
distinctYu <- distinct(dfYu, Gene, .keep_all=TRUE)
distinctKaya <- distinct(dfKaya, Gene, .keep_all=TRUE)
distinctGoyal <- distinct(dfGoyal, Gene, .keep_all=TRUE)
distinctZimmer <- distinct(dfZimmer, Gene, .keep_all=TRUE)
distinctLu <- distinct(dfLu, Gene, .keep_all=TRUE)
distinctSunMOLM <- distinct(dfSunMOLM, Gene, .keep_all=TRUE)
distinctBarb <- distinct(dfBarbMOLM, Gene, .keep_all=TRUE)
distinctSubo <- distinct(dfSubo, Gene, .keep_all=TRUE)

distinct_df_list_line <- list(distinctChandra, distinctGupta, distinctWei, distinctYu,
                              distinctLu, distinctSunMOLM, distinctBarb,
                              distinctSubo)

#MHB: I have a slightly different way that I would like to go about this cleaning business
#Please see the bottom of the code and I'm going to run through a couple different ways to
#clean these data. Mostly, I'm concerned about the NA values in the data, number of genes that don't have any variation, and creating a clean dataset to begin with before jumping right in. Take a look at the bottom. I'm going to use the same list distinct, but take a it a be farther. 



full_distinctDf_line <- reduce(distinct_df_list_line, full_join, by="Gene")
dimsdistinct_line <- dim(full_distinctDf_line)

print(paste("The full dataframe has", dimsdistinct_line[2], "columns. The labels should have",
            dimsdistinct_line[2]-1, "length. The lengths of the class labels and platforms labels are",
            length(allClassLabels_line),"and", length(allPlatformLabels_line)))

#writing all data to files for easy use


write.table(allClassLabels_line, "Processed_data/classlabels_line.csv",
            row.names=FALSE,col.names=FALSE,quote=FALSE, sep=",")

write.table(allPlatformLabels_line, "Processed_data/platformlabels_line.csv",
            row.names=FALSE,col.names=FALSE,quote=FALSE, sep=",")

#dataframes
write_csv(full_distinctDf_line, "Processed_data/cleandata_combined_line.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfChandra, "Processed_data/cleandata_Chandra.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfGupta, "Processed_data/cleandata_Gupta.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfWei, "Processed_data/cleandata_Wei.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfYu, "Processed_data/cleandata_Yu.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfKaya, "Processed_data/cleandata_Kaya.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfGoyal, "Processed_data/cleandata_Goyal.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfZimmer, "Processed_data/cleandata_Zimmer.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfLu, "Processed_data/cleandata_Lu.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfSun, "Processed_data/cleandata_Sun_complete.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfSunMOLM, "Processed_data/cleandata_SunMOLM.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfSunTHP, "Processed_data/cleandata_SunTHP.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfBarbMOLM, "Processed_data/cleandata_Barb_MOLM.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfBarbU937, "Processed_data/cleandata_Barb_U937.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfSubo, "Processed_data/cleandata_Subo.csv", append=FALSE,
          col_names=TRUE,quote="none")



#===============================================================================
#                           HMC3 & VALIDATION DATASETS 
#                       HMC3, Primary microglia and astrocyte
#===============================================================================

#-------------------------------------------------------------------------------
#                    CREATE ABUD DF - microglia and monocytes
#-------------------------------------------------------------------------------
rawDfAbud <- rawDataAbud

#rename samples so that cell type and condition are known
dfAbud_renamed <- rawDfAbud %>%
  rename(imgl.c1=GSM2360252, imgl.c2=GSM2360253, imgl.c3=GSM2360254, imgl.c4=GSM2360255,
         imgl.c5=GSM2360256, imgl.c6=GSM2360257, ihpc.c1=GSM2360258, ihpc.c2=GSM2360259,
         ihpc.c3=GSM2360260, ips.c1=GSM2360261, ips.c2=GSM2360262,ips.c3=GSM2360263,
         ips.c4=GSM2360264, monocyteclassic.c1=GSM2360265, monocyteclassic.c2=GSM2360266,
         monocyteclassic.c3=GSM2360267, monocyteclassic.c4=GSM2360268, monocyteclassic.c5=GSM2360269,
         monocyte.nonclassic.c1=GSM2360270, monocyte.nonclassic.c2=GSM2360271,
         monocyte.nonclassic.c3=GSM2360272, monocyte.nonclassic.c4=GSM2360273, fetalmic.c1=GSM2360274,
         fetalmic.c2=GSM2360275, fetalmic.c3=GSM2360276, adultmic.c1=GSM2360277, adultmic.c2=GSM2360278,
         adultmic.c3=GSM2360279, myeloid.dendritic.c1=GSM2360280, myeloid.dendritic.c2=GSM2360281,
         myeloid.dendritic.c3=GSM2360282, rat.imgl.c1=GSM2360283, rat.imgl.c2=GSM2360284,
         rat.imgl.c3=GSM2360285, rat.imgl.c4=GSM2360286, rat.imgl.c5=GSM2360287, rat.imgl.c6=GSM2360288,
         imgl.c7=GSM2445478, imgl.c8=GSM2445479, imgl.c9=GSM2445480, imgl.nTGFB1=GSM2551704,
         imgl.nTGFB2=GSM2551705, imgl.nTGFB3=GSM2551706)

#get HUGO gene names
genes_onlyDf <- select(human_geneannotation, GeneID:Symbol)
gene_names_abud <- pull(genes_onlyDf, Symbol)

extragenesDf_abud <- merge(genes_onlyDf, dfAbud_renamed, by="GeneID", sort=FALSE)

#check that gene names are correct 
sum(gene_names_abud == extragenesDf_abud$Symbol) == dim(extragenesDf_abud)[1]

#take out extra column
correctgenesDf_abud <- select(extragenesDf_abud, -GeneID) %>%
  rename(Gene=Symbol)

#only keep wanted samples for training (microglia and monocytes CD14+/CD16-)
training_samples <- c()
sample_names_abud <- colnames(correctgenesDf_abud)
for (name in sample_names_abud)
{
  if (str_detect(name,"mic") | str_detect(name, "monocyteclassic")) {
    training_samples = c(training_samples, name)
  }
}

trainingDfAbud <- dplyr::select(correctgenesDf_abud, c("Gene",training_samples))


#get IMGL samples for testing
imgl_samples_abud <- c()
sample_names_abud <- colnames(correctgenesDf_abud)
imgl_samples_abud <- sample_names_abud[str_detect(sample_names_abud,"imgl.c") & !str_starts(sample_names_abud, "rat")]
imglDfAbud <- dplyr::select(correctgenesDf_abud,c("Gene",imgl_samples_abud))

#final 
dfAbud <- trainingDfAbud

## final df is dfAbud (microglia and monocyte samples) and imglDfAbud 


#-------------------------------------------------------------------------------
#              CREATE ADAB DF - U87 & THP1 - HCM3 & iMGL samples
#-------------------------------------------------------------------------------
#load in data
rawDfADAB <- rawDataADAB
rawDfADAB <- rawDfADAB[,-1] #removes ensembl gene ID

#change column from hgnc to Gene title 
colNamesRawADAB <- colnames(rawDfADAB)
sampleNamesRawADAB <- colNamesRawADAB[2:67] #keeps all sample names
colnames(rawDfADAB) <- c("Gene", sampleNamesRawADAB) #replaces col with genes + original sample names

#replace Gene names containing - or . to underscores
rawDfADAB$Gene <- gsub('-', '_', rawDfADAB$Gene)
rawDfADAB$Gene <- gsub('\\.','_', rawDfADAB$Gene)

#remove iPSC and iHPC
donotkeep_samples <- c()
sample_names <- colnames(rawDfADAB)
for (name in sample_names)
{
  if (str_detect(name,"iPSC") | str_detect(name, "iHPC") | str_detect(name, "PBMC")) {
    donotkeep_samples = c(donotkeep_samples, name)
  }
}
newDfADAB <- dplyr::select(rawDfADAB, -donotkeep_samples)

#Separate the HMC3 cell line as its own df for later tests 
HMC3_df <- dplyr::select(newDfADAB,ADAB_RNA_HMC3_WT_NOAB_NONE_NONE_NOsrm_24hr_S_1.1.1:ADAB_RNA_HMC3_WT_OLGO_01uM_NONE_NOsrm_24hr_S_1.2.1)
dfADAB <- dplyr::select(newDfADAB, -colnames(HMC3_df)) #new df ADAB data without HMC3 line
#^I could clean this up with an str_detect function 

#verify HMC3 is gone - should be a difference of 4 columns
paste("Original ADAB df:", dim(newDfADAB), "New ADAB df dims:", dim(dfADAB), "HMC3 df dims:", dim(HMC3_df), sep=" ")

##take out iMGL cells
imgl_names = c()
firstcolnames <- colnames(dfADAB)
len_names <- length(firstcolnames)
for (i in 1:len_names) {
  if (str_detect(firstcolnames[i], regex("iMGL", ignore_case=TRUE))) {
    imgl_names = c(imgl_names, firstcolnames[i])
  }
}

dfImgl <- dplyr::select(dfADAB, imgl_names)
dfADAB_noimgl <- dplyr::select(dfADAB, -imgl_names)

#verify iMGL is gone - should be a difference of 8 columns
paste("Original ADAB df:", dim(dfADAB), "New ADAB df dims:", dim(dfADAB_noimgl), "iMGL df dims:", dim(dfImgl), sep=" ")

dfADAB <- dfADAB_noimgl

#grab and change column names for class labels 
lastN <- ncol(dfADAB) #number of columns
columnNamesADAB = str_split_fixed(colnames(dfADAB), "_", n= Inf)[,3] #takes cell type from column
columnNamesADAB = columnNamesADAB[2:lastN] #takes cell col names, excludes "Gene" col
#unique(columnNamesADAB)

#label cell type
len = length(columnNamesADAB)
for (i in 1:len) {
  if (str_detect(columnNamesADAB[i], regex("u\\d", ignore_case=TRUE))) {
    columnNamesADAB[i] = "U87"
  }
  if (str_detect(columnNamesADAB[i], regex("THP", ignore_case=TRUE))) {
    columnNamesADAB[i] = "THP1"
  }
}

#create class labels for ADAB data for cell line classifier
classLabelsADAB = columnNamesADAB

#create platform labels for ADAB data
lenA <- length(classLabelsADAB)
platformADAB <- rep("Quiroga Data", lenA)

##df's = HMC3_df, dfImgl
##labels: classLabelsADAB and platformADAB 
##df is dfADAB


#-------------------------------------------------------------------------------
#                       CREATE CHAI DF - HMC3 samples
#-------------------------------------------------------------------------------
#load in data
rawChai <- rawDataChai

#only keep wanted samples for training (control HMC3 samples)
control_hmc3_samples <- c("GSM4703735", "GSM4703736", "GSM4703737")
controlonly_chai <- dplyr::select(rawChai, c("GeneID", any_of(control_hmc3_samples)))


#get HUGO gene names
genes_onlyDf <- select(human_geneannotation, GeneID:Symbol)
gene_names <- pull(genes_onlyDf, Symbol)

extragenesDf_chai <- merge(genes_onlyDf, controlonly_chai, by="GeneID", sort=FALSE)

#check that gene names are correct 
sum(gene_names == extragenesDf_chai$Symbol) == dim(extragenesDf_chai)[1]

#take out extra column
dfChai <- select(extragenesDf_chai, -GeneID) %>%
  rename(Gene=Symbol)

#replace Gene names containing - or . to underscores
dfChai$Gene <- gsub('-', '_', dfChai$Gene)
dfChai$Gene <- gsub('\\.','_', dfChai$Gene)

##final: dfChai
# Contains 3 HMC3 samples 



#-------------------------------------------------------------------------------
#                       CREATE ARMAN DF - HMC3 samples
#-------------------------------------------------------------------------------
#data = rawDataArman

#only keep wanted samples for training (control HMC3 samples)
controlArman <- rawDataArman %>%
  select(matches("Name|VEH")) %>%
  rename(Gene = Name) %>%
  relocate(Gene, .before=1)

#replace Gene names containing - or . to underscores
controlArman$Gene <- gsub('-', '_', controlArman$Gene)
controlArman$Gene <- gsub('\\.','_', controlArman$Gene)

dfArman <- controlArman

##final: dfArman
# Contains 3 HMC3 samples 



#-------------------------------------------------------------------------------
#                 CREATE ROWITCH DF - astrocyte validation
#-------------------------------------------------------------------------------

#laod in data
rawDfRow <- rawDataRow
#replace Gene names containing - or . to underscores
rawDfRow$Gene <- gsub('-', '_', rawDfRow$Gene)
rawDfRow$Gene <- gsub('\\.','_', rawDfRow$Gene)

#take out diseased types 
dimsRow <- dim(rawDfRow)
cellTypeRow <- colnames(rawDfRow) #takes column names which have cell type
cellTypeRow <- cellTypeRow[2:dimsRow[2]] #removes "Gene" name

lenRow = length(cellTypeRow)
diseased = c()
for (i in 1:lenRow) {
  if (str_detect(cellTypeRow[i], regex("MS", ignore_case=TRUE))) {
    diseased = c(diseased, cellTypeRow[i])
  }
}

#dataframe with only control samples, no diseased 
controlRow <- dplyr::select(rawDfRow, -diseased)

#grab only astrocyte samples 
dimsCRow <- dim(controlRow)
cellTypeCRow <- colnames(controlRow) #takes column names which have cell type

lenCRow = length(cellTypeCRow)
row_astro_samples = c()
for (i in 1:lenCRow) {
  if (str_detect(cellTypeCRow[i], regex("astro", ignore_case=TRUE))) {
    row_astro_samples = c(row_astro_samples, cellTypeCRow[i])
  }
}

astro_setRow <- dplyr::select(controlRow, c(Gene,row_astro_samples))

##final: astro_setRow or controlRow
# Contains all astrocyte primary cell data, controlRow has many cell types



#-------------------------------------------------------------------------------
#                   CREATE SAGAR DF - microglia validation
#-------------------------------------------------------------------------------

#replace Gene names containing - or . to underscores
rawDataSagar$Gene <- gsub('-', '_', rawDataSagar$Gene)
rawDataSagar$Gene <- gsub('\\.','_', rawDataSagar$Gene)

dfSagar <- rawDataSagar

##final: dfSagar
# Contains all microglia samples 



#=========================================================================================
#          WRITING ALL TEST DATA TO FILES
#=========================================================================================

#ORDER

# ALL VARIABLES:

##final: dfChai
# Contains 3 HMC3 samples 

##final: astro_setRow or controlRow
# Contains all astrocyte primary cell data, controlRow has many cell types

##final: dfSagar
# Contains all microglia samples 

##final from ADAB data: HMC3_df & dfImgl
# HMC3_df contains all HMC3 samples. dfImgl contains all iMGL samples 


#validation/test sets
write_csv(dfAbud, "Processed_data/cleandata_micmono_abud.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(imglDfAbud, "Processed_data/cleandata_imgl_abud.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfChai, "Processed_data/cleandata_Hmc3_chai.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfArman, "Processed_data/cleandata_Hmc3_arman.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfADAB, "Processed_data/cleandata_Adab.csv", append=FALSE,
          col_names=TRUE,quote="none")

write.table(classLabelsADAB, "Processed_data/classlabels_adab.csv",
            row.names=FALSE,col.names=FALSE,quote=FALSE, sep=",")

write_csv(HMC3_df, "Processed_data/cleandata_Hmc3_adab.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfImgl, "Processed_data/cleandata_imgl_adab.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(controlRow, "Processed_data/dfRow_controlsamples.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(astro_setRow, "Processed_data/dfRow_astrocyteonly.csv", append=FALSE,
          col_names=TRUE,quote="none")

write_csv(dfSagar, "Processed_data/dfSagar_microglia.csv", append=FALSE,
          col_names=TRUE,quote="none")





#MHB: Creating a new final dataset with a couple more stipulations. 
distinct_df_list #This is all of the tissue-based datasets 
distinct_df_list_line


#MHB: Tissue classifier: No na.genes and 80% sample coverage* (*nonzero expression)
distinct_df_list_MB = lapply(distinct_df_list, function(x){
     numsamp = floor((dim(x)[2]-1)*0.8) #80% of samples with non-zero RNA for a gene (-1 gets rid of Gene column)
     narows = -which(is.na(x$Gene))
     if(length(narows) != 0){
         x = x[narows,]
     }
     keeprows = which(rowSums(x[-which(colnames(x) %in% "Gene")] != 0) > numsamp) #not >= to get above floor of 80%. 
     return(x[keeprows,])
})

#MHB: Now I'm want to extract the itersection of genes in all of these datasets. 
distinct_df_list_MBg = unlist(sapply(distinct_df_list_MB, function(x){
    x$Gene
}))

distinct_df_list_MBg = unlist(sapply(distinct_df_list_MB, function(x){
    x$Gene
}))
gcounts = sort(table((distinct_df_list_MBg)))
tissuesgenes = names(gcounts)[which(gcounts == length(distinct_df_list_MB))]


#MHB: Cell-line classifier: No no.genes and 80 sample coverage* (*nonzero expression)
distinct_df_list_line_MB = lapply(distinct_df_list_line, function(x){
     numsamp = floor((dim(x)[2]-1)*0.9) #80% of samples with non-zero RNA for a gene (-1 gets rid of Gene column)
     narows = -which(is.na(x$Gene))
     if(length(narows) != 0){
         x = x[narows,]
     }
     keeprows = which(rowSums(x[-which(colnames(x) %in% "Gene")] != 0) > numsamp) #not >= to get above floor of 80%. 
     return(x[keeprows,])
})
#MHB: Now I'm going to extract the genes that are in all of these datasets that have at least 80% coverage. 
distinct_df_list_line_MBg = unlist(sapply(distinct_df_list_line_MB, function(x){
    x$Gene
}))
gcounts_line = sort(table((distinct_df_list_line_MBg)))
cellgenes = names(gcounts_line)[which(gcounts_line == length(distinct_df_list_line_MB))]


#MHB: Genes that are well-covered in both sets
same_gene = intersect(cellgenes,tissuesgenes)


#MHB: Now I want to make a geneset with fewer genes that have a good amount of overlap in these studies to be safe. 
clean_gene_tissue_dfs = lapply(distinct_df_list_MB, function(x){
    return(x %>% filter(Gene %in% same_gene))
})

clean_gene_cell_dfs = lapply(distinct_df_list_line_MB, function(x){
    return(x %>% filter(Gene %in% same_gene))
})



#MHB: So here is the note for the paper. 8723 genes were included. We will rely on MCP to test the gene pairs in the testing data. 
tissuenames <- c("Hansen", "Gal", "PrimaryBrain", "Goss", "Verd")
if(length(tissuenames) != length(distinct_df_list )){
    stop("It looks like you made a change to which tissue-line datasets that were include")
}

cellnames <- c("Chandra", "Gupta", "Wei", "Yu", "Lu", "SunMOLM", "Barb","Subo")
if(length(cellnames) != length(distinct_df_list_line)){
    stop("It looks like you made a change to which cell-line datasets that were include")
}


full_cleanMHB_tissue <- reduce(clean_gene_tissue_dfs, full_join, by="Gene")
full_cleanMHB_line <- reduce(clean_gene_cell_dfs, full_join, by="Gene")

write.table(full_cleanMHB_tissue, "Processed_data/cleandata_combined_training_tissue.MHB.csv", append=FALSE, quote=F, row.names=F,sep=",")
write.table(full_cleanMHB_line, "Processed_data/cleandata_combined_training_line.MHB.csv", append=FALSE, quote=F, row.names=F,sep=",")

