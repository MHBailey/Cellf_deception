###############################################################################
###############################################################################
#Creator: Matthew Bailey
#Last Updated June 16, 2025
#MCP_Microglial_Classifier for Cell-type classifier for HMC3 cells
#Published: doi:
###############################################################################
###############################################################################

#This is too is meant to classify tissue of origin of cell lines from depmap. 
#I'm going to start by trying to classify BRAIN from Myeloid lineages 
#And I'll probably try to do a full comparision across all tissue of origin to 
#See how it works in the DepMap dataset. 

#NOTE: I think this would also be interesting to see if the BatchCorrected or 
#Raw data perform better, but I likely don't have time for that. 


#==============================================================================
##---------------LIBRARIES-----------------------------------------------------
#==============================================================================
library(tidyverse)
library(data.table)
library(dplyr)
library(utils)
library(multiclassPairs)
library(reshape2)
library(RColorBrewer)
library(Biobase)
library(switchBox)


#-------------------------------------------------------------------------------
##----------------------CREATING DATA MATRIX FOR MCP----------------------------
#-------------------------------------------------------------------------------

#Depmap files
meta = "Data/Model.csv"

depmap_file_path = "Data/OmicsExpressionProteinCodingGenesTPMLogp1.csv"

depmap <- tryCatch({
  read.csv(depmap_file_path)
}, error = function(e) {
  message("ERROR: File '", depmap_file_path, "' not found. You'll need to download this from DEPMAP. Likely at https://depmap.org/portal/data_page/?tab=allData")
  NULL  # Return NULL or handle accordingly
})

#Genelist Restriction: 
comboDf_fname_MHB = "Processed_data/cleandata_combined_training_tissue.MHB.csv"

#tests and validations 
hmc3_df_fname = "Processed_data/cleandata_Hmc3_adab.csv"
adabDf_fname = "Processed_data/cleandata_Adab.csv"
imglDf_fname = "Processed_data/cleandata_imgl_adab.csv"

#Read many of the other datasets 
hmc3_chai = read_csv("Processed_data/cleandata_Hmc3_chai.csv")
astro_setRow = read_csv("Processed_data/dfRow_astrocyteonly.csv")
dfSagar = read_csv("Processed_data/dfSagar_microglia.csv")
dfAbud = read_csv("Processed_data/cleandata_micmono_abud.csv")
imgl_abud = read_csv("Processed_data/cleandata_imgl_abud.csv")
dfGupta = read_csv("Processed_data/cleandata_Gupta.csv")
dfChandra = read_csv("Processed_data/cleandata_Chandra.csv")

#ONE MORE TEST MHBAILEY that we did for cell lines but not primary
dfArman = read_csv("Processed_data/cleandata_Hmc3_arman.csv")


#-------------------------------------------------------------------------------
##------------d----------LOAD IN THE TRAINING DATA ------------------------------
#-------------------------------------------------------------------------------

metadf <- fread(meta)
metadf$Sex <- ifelse(metadf$Sex == "", "Unknown", metadf$Sex) #A little cleaning 

dmdf <- fread(depmap)
fullDf <- read_csv(comboDf_fname_MHB)

#Get the gene list that I used for training the first model
genelist_subset = fullDf$Gene

#Subset the depmap expression panel to these genes 
#Change the names first 
names(dmdf) = gsub(pattern = " \\(.*\\)", replacement = "", x = names(dmdf))
dmdf_subgene <- dmdf %>% select(any_of(c("V1",genelist_subset)))


#I need to transpose the data to have genes for rows and samples for columns. 
dmdf_trans <- dcast(reshape2::melt(dmdf_subgene), variable~V1)
dmdf_trans <- rename(dmdf_trans,"Gene"=variable)

#Submeta 
metadf_sub <- metadf %>% filter(ModelID %in% names(dmdf_trans))


#Create class and platform labels
L1 <- pull(metadf_sub,OncotreeLineage)
P1 <- pull(metadf_sub,Sex) 

#SUBSET TO CNS and Myeloid  #It's too big otherwise. 
NL1 = L1[sort(c(which(L1== "CNS/Brain"),which(L1 == "Myeloid")))]
wL1 = sort(c(which(L1== "CNS/Brain"),which(L1 == "Myeloid")))
WP1 = P1[wL1]


cols2get <- c(1,wL1+1)#This is offset the labels based on the first column 

dmdf_trans_cns_mye <- dmdf_trans[,cols2get]
#Turn this into a matrix for the MCP algorithm

fullMatrix <- as.matrix(dmdf_trans_cns_mye[,-1])
rownames(fullMatrix) = dmdf_trans_cns_mye$Gene
L1 = NL1 #Overwrite so I don't have to change code below
P1 = WP1 #Overwrite these variable so I don't need to change below 

#-------------------------------------------------------------------------------
##----------------------LOAD IN HCM3 TESTDATA ----------------------------------
#-------------------------------------------------------------------------------

#read in data
hmc3Df <- read_csv(hmc3_df_fname)
dfADAB <- read_csv(adabDf_fname)
imglDf <- read_csv(imglDf_fname)

#create the matrices 
#get gene names from main dataset
gene_names <- pull(dfADAB, "Gene")

#create HMC3 matrix
HMC3_matrix <- as.matrix(hmc3Df)
rownames(HMC3_matrix) <- gene_names

#create iMGL matrix
imgl_matrix <- as.matrix(imglDf)
rownames(imgl_matrix) <- gene_names

#-------------------------------------------------------------------------------
##---------------------- MULTICLASS PAIRS TRAINING -----------------------------
#-------------------------------------------------------------------------------


## create the data object
object <- ReadData(Data = fullMatrix,
                   Labels = L1,
                   Platform = P1,
                   verbose = FALSE)
object

#------------------------------------------------------
#------------------START TRAINING DATA-----------------
#------------------------------------------------------

## split the data
# 60% as training data and 40% as testing data

#get dims of matrix
n <- ncol(fullMatrix)
nr <- nrow(fullMatrix)

set.seed(1234)
training_samples <- sample(1:n, size = n*.6) #randomly pulls 60% of indices for training data

#split matrix into train and test matrix
train <- fullMatrix[1:nr,training_samples]
test  <- fullMatrix[1:nr, -training_samples]
#assign labels to split data
L1.train <- L1[training_samples]
P1.train <- P1[training_samples]
L1.test <- L1[-training_samples]
P1.test <- P1[-training_samples]
#dim(train)
#dim(test)

# to be sure there are no shared samples between the training and testing data
result = sum(colnames(test) %in% colnames(train)) == 0
sprintf("There is no training and testing data overlap. Our data: %s", result)

## create the data objects
objectTrain <- ReadData(Data = train,
                        Labels = L1.train,
                        Platform = P1.train,
                        verbose = FALSE)
objectTrain


# -------------------------------------------------------------------------
#TEST OBJECT

objectTest <- ReadData(Data = test,
                       Labels = L1.test,
                       Platform = P1.test,
                       verbose = FALSE)
objectTest

#------------------------------------------------------
#----------------------RANDOM FOREST SCHEME------------
#------------------------------------------------------

# sorting genes
genes_RF <- sort_genes_RF(data_object = objectTrain,
                          rank_data = TRUE,
                          platform_wise = TRUE,
                          num.trees = 2000, # more features, more tress are recommended. 
                          #14431 for all genes 
                          seed=123456, # for reproducibility
                          write.forest = TRUE,
                          verbose = TRUE)
genes_RF # sorted genes object

#-------------------RULES------------------------
# to get an idea of how many genes we will use
# and how many rules will be generated
#SO this is just a helpful function for me, the user to know how many genes
# I must have to get a specific # of rules 
#this is just to optimize our paramters 
summary_genes <- summary_genes_RF(sorted_genes_RF = genes_RF,
                                  genes_altogether = c(50, 75, 85,100,150,200),
                                  genes_one_vs_rest = c(50, 75, 85,100,150,200))
knitr::kable(summary_genes)

# to give enough number of  rules and unique genes for our classes
# Now let's run sort_rules_RF to create the rules and sort them
rules_RF <- sort_rules_RF(data_object = objectTrain,
                          sorted_genes_RF = genes_RF,
                          genes_altogether = 85,
                          genes_one_vs_rest = 100,
                          num.trees = 1000,
                          seed=123456,
                          verbose = TRUE)



rules_RF # sorted rules object



#------------------MODEL TRAINING----------------------

# prepare the simple data.frame for the parameters I want to test
# names of arguments as column names
# this df has three sets (3 rows) of parameters
parameters <- data.frame(
  gene_repetition=c(1,1,1),
  rules_one_vs_rest=c(85, 100, 100),
  rules_altogether=c(100, 85,100),
  run_boruta=c(TRUE, TRUE,TRUE),
  plot_boruta = FALSE,
  num.trees=c(1000,1000,1000),
  stringsAsFactors = FALSE)

# parameters
# for overall and byclass possible options, check the help files
para_opt <- optimize_RF(data_object = objectTrain,
                        sorted_rules_RF = rules_RF,
                        parameters = parameters,
                        test_object = NULL,
                        overall = c("Accuracy","Kappa"), # wanted overall measurements 
                        byclass = c("F1"), # wanted measurements per class
                        verbose = TRUE)

para_opt # results object
para_opt$summary # the df of with summarized information

# train the final model
# it is preferred to increase the number of trees and rules in case you have
# large number of samples and features
# for quick example, we have small number of trees and rules here
# based on the optimize_RF results we will select the parameters
RF_classifier <- train_RF(data_object = objectTrain,
                          sorted_rules_RF = rules_RF,
                          gene_repetition = 1,
                          rules_altogether = 100,
                          rules_one_vs_rest = 100,
                          run_boruta = TRUE,
                          plot_boruta = FALSE,
                          probability = TRUE,
                          num.trees = 1000,
                          boruta_args = list(),
                          verbose = TRUE)

#-------------------model TRAINING---------------------------
# plot proximity matrix of the out-of-bag samples
# Note: this takes a lot of time if the data is big

pdf("Figures/PredAccur_depmap_training.pdf",heigh=10, width=10,useDingbats=F)
proximity_matrix_RF(object = objectTrain,
                    classifier = RF_classifier,
                    plot = T,
                    return_matrix = FALSE, # if we need to get the matrix itself
                    title = "Prediction Accuracy on Training Data",
                    cluster_cols = TRUE)

dev.off()

#-------------------------------------------------------------
#----------------TRAINING ACCURACY----------------------------
#-------------------------------------------------------------

# get the prediction labels from the trained model
# if the classifier trained using probability = FALSE
training_pred <- RF_classifier$RF_scheme$RF_classifier$predictions
# if (is.factor(training_pred)) {
#   x <- as.character(training_pred)
# }

# if the classifier trained using probability = TRUE
if (is.matrix(training_pred)) {
  x <- colnames(training_pred)[max.col(training_pred)]
}

# training accuracy
caret::confusionMatrix(data =factor(x),
                       reference = factor(L1.train),
                       mode = "everything")


#-------------------------------------------------------------
#---------------PREDICTIONS FOR TEST
#------------------------------------
#-------------------------------------------------------------

# apply on test data
results <- predict_RF(classifier = RF_classifier,
                      Data = objectTest,
                      impute = TRUE) # can handle missed genes by imputation


# get the prediction labels
# if the classifier trained using probability   = FALSE
test_pred <- results$predictions
# if (is.factor(test_pred)) {
#   x2 <- as.character(test_pred)
# }

# if the classifier trained using probability   = TRUE
if (is.matrix(test_pred)) {
  x2 <- colnames(test_pred)[max.col(test_pred)]
}

# training accuracy
caret::confusionMatrix(data = factor(x2),
                       reference = factor(L1.test),
                       mode = "everything")

pdf("Figures/PredAccur_depmap_testing.pdf",height=10,width=10,useDingbats=F)
proximity_matrix_RF(object = objectTest,
                    classifier = RF_classifier,
                    plot = TRUE,
                    return_matrix = FALSE, # if we need to get the matrix itself
                    title = "Prediction Accuracy on Testing Data",
                    cluster_cols = FALSE)
dev.off()

#------------------VISUALIZATION--------------
#visualize the binary rules in training dataset
pdf("Figures/PredAccur_depmap_training_rules.pdf",height=10,width=10,useDingbats=F)
plot_binary_RF(Data = objectTrain,
               classifier = RF_classifier,
               prediction = NULL,
               as_training = TRUE, # to extract the scores from the model
               show_scores = TRUE,
               top_anno = "ref",
               show_predictions = TRUE,
               cluster_cols = TRUE, #otherwise there's an issue with clustering too small
               #margin = c(0,5,0,8),
               title = "Training data")
dev.off()

# visualize the binary rules in testing dataset
pdf("Figures/PredAccur_depmap_testing_rule.pdf",height=10,width=10,useDingbats=F)
plot_binary_RF(Data = objectTest,
               ref = L1.test, # Get ref labels from the test data
               classifier = RF_classifier,
               prediction = results,
               as_training = FALSE,
               show_scores = TRUE,
               top_anno = "ref",
               show_predictions = TRUE,
               title = "Testing data",
               cluster_cols = FALSE)
dev.off()




#=========================================================================================
#                           TEST & VALIDATION DATASETS 
#                       2 HMC3, astrocyte, microglia
#=========================================================================================

#-------------------------------------------------------------------------------
#                     HMC3 samples from Chai data
#-------------------------------------------------------------------------------
#read in data
#hmc3_chai
#create the matrix 
#get gene names from main dataset
hmc3_chai = read_csv("Processed_data/cleandata_Hmc3_chai.csv")
gene_names_chai <- pull(hmc3_chai, "Gene")

#create HMC3 chai matrix
prematrix_chai <- select(hmc3_chai, -"Gene")
HMC3_matrix_chai <- as.matrix(prematrix_chai)
rownames(HMC3_matrix_chai) <- gene_names_chai

#    HMC3 CHAI DATA
#-------------------------------------------------------------------------------
objectHMC3chai <- ReadData(Data = HMC3_matrix_chai,
                           Labels = c("HMC3", "HMC3", "HMC3"),
                           Platform = NULL,
                           verbose = FALSE)

resultsHMC3chai <- predict_RF(classifier = RF_classifier,
                              Data = objectHMC3chai,
                              impute = TRUE) # can handle missed genes by imputation

hmc3_pred_chai <- resultsHMC3chai$predictions

# # if the classifier trained using probability   = TRUE
if (is.matrix(hmc3_pred_chai)) {
  x3chai <- colnames(hmc3_pred_chai)[max.col(hmc3_pred_chai)]
}


hmc3_predlonger_chai <- reshape2::melt(hmc3_pred_chai)%>%
  mutate(fixedValues = round(value, 3))

sample_names_chai = colnames(hmc3_chai)[2:4]

#use barplot to visualize predictions by numbers

chaiplot <- ggplot(hmc3_predlonger_chai, aes(x=Var2, y=value, fill=Var1), color="black") +
  labs(title = "HMC3 Depmap Cell Prediction (Chai Data)", x="Cell Type", y="Cell Prediction Score") +
  geom_col(position="dodge") + theme_bw() +
  geom_text(aes(label = fixedValues), vjust=-0.5, position = position_dodge(.9), size=1.79) +
  scale_x_discrete(labels= c("astrocyte"="Astrocyte",
                             "microglia" = "Microglia",
                             "Whole.cortex" = "Whole cortex",
                             "oligodendrocyte" = "Oligodendrocyte")) +
  ylim(0,1) +
  theme(plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = c(.72,.7), legend.text=element_text(size=6.5)) +
  scale_fill_viridis_d(labels = sample_names_chai)


pdf("Figures/Depmap_chi_plot.pdf",height=6,width=6,useDingbats=F)
print(chaiplot)
dev.off()

#-------------------------------------------------------------------------------
#                     HMC3 & iMGL samples from Quiroga data
#-------------------------------------------------------------------------------
#read in data
hmc3Df <- read_csv(hmc3_df_fname)
dfADAB <- read_csv(adabDf_fname)
imglDf <- read_csv(imglDf_fname)

#create the matrices 
#get gene names from main dataset
gene_names <- pull(dfADAB, "Gene")

#create HMC3 matrix
HMC3_matrix <- as.matrix(hmc3Df)
rownames(HMC3_matrix) <- gene_names

#create iMGL matrix
imgl_matrix <- as.matrix(imglDf)
rownames(imgl_matrix) <- gene_names


#    HMC3 DATA
#-------------------------------------------------------------------------------
objectHMC3 <- ReadData(Data = HMC3_matrix,
                       Labels = c("HMC3", "HMC3", "HMC3", "HMC3"),
                       Platform = NULL,
                       verbose = FALSE)

resultsHMC3 <- predict_RF(classifier = RF_classifier,
                          Data = objectHMC3,
                          impute = TRUE) # can handle missed genes by imputation

hmc3_pred <- resultsHMC3$predictions

# # if the classifier trained using probability   = TRUE
if (is.matrix(hmc3_pred)) {
  x3 <- colnames(hmc3_pred)[max.col(hmc3_pred)]
}


# # show results
# visualize the binary rules in testing dataset
# plot_binary_RF(Data = objectHMC3,
#                ref = c("HMC3", "HMC3","HMC3","HMC3"), # Get ref labels from the test data
#                classifier = RF_classifier,
#                prediction = resultsHMC3, 
#                as_training = FALSE, 
#                show_scores = TRUE,
#                top_anno = "ref",
#                show_predictions = TRUE,
#                title = "HMC3 data",
#                cluster_cols = FALSE)

hmc3_predlonger <- reshape2::melt(hmc3_pred)%>%
  mutate(fixedValues = round(value, 3))



#use barplot to visualize predictions by numbers
hmc3plot <- ggplot(hmc3_predlonger, aes(x=Var2, y=value, fill=Var1), color="black") +
  labs(title = "HMC3 Depmap Cell Prediction (Quiroga Data)", x="Cell Type", y="Cell Prediction Score") +
  geom_col(position="dodge") + theme_bw() +
  geom_text(aes(label = fixedValues), vjust=-0.5, position = position_dodge(.9), size=1.79) +
  scale_x_discrete(labels= c("astrocyte"="Astrocyte",
                             "microglia" = "Microglia",
                             "Whole.cortex" = "Whole cortex",
                             "oligodendrocyte" = "Oligodendrocyte")) +
  ylim(0,1) +
  theme(plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = c(.72,.7), legend.text=element_text(size=6.5)) +
  scale_fill_viridis_d(labels =c("HMC3 ctl 1", "HMC3 1uM AB 1", "HMC3 ctl 2", "HMC3 1uM AB 2"))


pdf("Figures/Depmap_hmc3_plot.pdf",height=6,width=6,useDingbats=F)
print(hmc3plot)
dev.off()


#    IMGL DATA
#-------------------------------------------------------------------------------
objectIMGL <- ReadData(Data = imgl_matrix,
                       Labels = c("iMGL", "iMGL", "iMGL", "iMGL", "iMGL", "iMGL", "iMGL", "iMGL"),
                       Platform = NULL,
                       verbose = FALSE)

#results
resultsIMGL <- predict_RF(classifier = RF_classifier,
                          Data = objectIMGL,
                          impute = TRUE) # can handle missed genes by imputation

imgl_pred <- resultsIMGL$predictions

# # if the classifier trained using probability   = TRUE
if (is.matrix(imgl_pred)) {
  x3 <- colnames(imgl_pred)[max.col(imgl_pred)]
}

imgl_predlonger <- reshape2::melt(imgl_pred)%>%
  mutate(fixedValues = round(value, 3))


#use barplot to visualize predictions by numbers
imglplot <- ggplot(imgl_predlonger, aes(x=Var2, y=value, fill=Var1), color="black") +
  labs(title = "iMGL Depmap Cell Prediction (Quiroga Data)", x="Cell Type", y="Cell Prediction Score") +
  geom_col(position="dodge") + theme_bw() +
  geom_text(aes(label = fixedValues), vjust=-0.5, position = position_dodge(.9), size=1.79) +
  scale_x_discrete(labels= c("astrocyte"="Astrocyte",
                             "Microglia" = "Microglia",
                             "Whole.cortex" = "Whole cortex",
                             "oligodendrocyte" = "Oligodendrocyte")) +
  ylim(0,1) +
  theme(plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = "right", legend.text=element_text(size=6.5)) +
  scale_fill_viridis_d(labels =c("Sample 1","Sample 2","Sample 3","Sample 4","Sample 5",
                                 "Sample 6","Sample 7","Sample 8"))

pdf("Figures/Depmap_imgl_plot.pdf",height=6,width=6,useDingbats=F)
print(imglplot)
dev.off()


#-------------------------------------------------------------------------------
#                     ASTROCYTE samples from Rowitch data
#-------------------------------------------------------------------------------
#create the matrix
gene_names_row <- pull(astro_setRow, "Gene")

preRow <- dplyr::select(astro_setRow, -"Gene")
astroMatrixRow <- as.matrix(preRow)
rownames(astroMatrixRow) <- gene_names_row

#create MCP object 
lengthRow = ncol(preRow)
labels_astroRow = rep("Astrocyte", lengthRow)

object_astroRow <- ReadData(Data = astroMatrixRow,
                            Labels = labels_astroRow,
                            Platform = NULL,
                            verbose = FALSE)

#ROW RF DATA ---------------
resultsRow_astro <- predict_RF(classifier = RF_classifier,
                               Data = object_astroRow,
                               impute = TRUE,
                               impute_reject = .94)
#can handle missed genes by imputation

astroRow_pred <- resultsRow_astro$predictions

# # if the classifier trained using probability   = TRUE
if (is.matrix(astroRow_pred)) {
  x3 <- colnames(astroRow_pred)[max.col(astroRow_pred)]
}

astroRow_predlonger <- reshape2::melt(astroRow_pred)%>%
  mutate(fixedValuesRowAstro = round(value, 3))

#use barplot to visualize predictions by numbers
astroplot <- ggplot(astroRow_predlonger, aes(x=Var2, y=value, fill=Var1), color="black") +
  labs(title = "Astrocyte Depmap Cell Prediction (Schirmer Data)", x="Cell Type", y="Cell Prediction Score") +
  geom_col(position="dodge") + theme_bw() +
  geom_text(aes(label = fixedValuesRowAstro), vjust=-0.5, position = position_dodge(.9), size=1.5) +
  scale_x_discrete(labels= c("astrocyte"="Astrocyte",
                             "Microglia" = "Microglia",
                             "Whole.cortex" = "Whole cortex",
                             "oligodendrocyte" = "Oligodendrocyte")) +
  theme(plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position.inside = c(.72,.7), legend.text=element_text(size=6.5)) +
  scale_fill_viridis_d(labels =c("C7Astroctye", "C9Astroctye", "C2Astroctye", "C8Astroctye",
                                 "C4Astroctye", "C5Astroctye", "C1Astroctye", "C6Astroctye", "C3Astroctye"))

pdf("Figures/Depmap_astro_plot.pdf",height=6,width=6,useDingbats=F)
print(astroplot)
dev.off()

#-------------------------------------------------------------------------------
#                     MICROGLIA samples from Sagar data
#-------------------------------------------------------------------------------
#create matrix
gene_names_sagar <- pull(dfSagar, "Gene")

noNameSagar <- dfSagar %>%
  dplyr::select(-"Gene")

sagarMatrix <- as.matrix(noNameSagar)
rownames(sagarMatrix) <- gene_names_sagar

lengthSag = ncol(noNameSagar)
labels = rep("Microglia", lengthSag)

objectSagar <- ReadData(Data = sagarMatrix,
                        Labels = labels, #cell line we know it's not similar to
                        Platform = NULL,
                        verbose = FALSE)

#SAGAR DATA ---------------
resultsSagar <- predict_RF(classifier = RF_classifier,
                           Data = objectSagar,
                           impute = TRUE,
                           impute_reject = 0.94) # can handle missed genes by imputation

sagar_pred <- resultsSagar$predictions

# # if the classifier trained using probability   = TRUE
if (is.matrix(sagar_pred)) {
  x3 <- colnames(sagar_pred)[max.col(sagar_pred)]
}

sagar_predlonger <- reshape2::melt(sagar_pred)%>%
  mutate(fixedValuesSag = round(value, 3))

#use barplot to visualize predictions by numbers
sagarplot <- ggplot(sagar_predlonger, aes(x=Var2, y=value, fill=Var1), color="black") +
  labs(title = "Depmap Cell Prediction (Masuda Data)", x="Cell Type", y="Cell Prediction Score") +
  geom_col(position="dodge") + theme_bw() +
  geom_text(aes(label = fixedValuesSag), vjust=-0.5, position = position_dodge(0.9), size=1) +
  scale_x_discrete(labels= c("astrocyte"="Astrocyte",
                             "microglia" = "Microglia / Macrophage",
                             "whole.cortex" = "Whole cortex",
                             "oligodendrocyte" = "Oligodendrocyte")) +
  theme(plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = c(.72,.7), legend.text=element_text(size=6.5)) +
  scale_fill_viridis_d(labels =c("Pat3.Myeloid", "Pat1.Myeloid", "Pat2.Myeloid", "Pat5.Myeloid",
                                 "MS2.Myeloid", "MS4.Myeloid", "MS5.Myeloid"))


pdf("Figures/Depmap_sagar_plot.pdf",height=6,width=6,useDingbats=F)
print(sagarplot)
dev.off()

#-------------------------------------------------------------------------------
#                     Micro & iMGL samples from ABUD data
#-------------------------------------------------------------------------------
#read in data
#dfAbud and imgl_abud

#create the matrix 
#get gene names from main dataset
gene_names_abud <- pull(dfAbud, "Gene")

#create matrix for the microglia/monocyte data 
prematrix_abud <- select(dfAbud, -"Gene")
matrix_abud <- as.matrix(prematrix_abud)
rownames(matrix_abud) <- gene_names_abud
sample_names_abudmm <- colnames(dfAbud)[2:12]

#create matrix for the iMGL data 
prematrix_abud_imgl <- select(imgl_abud, -"Gene")
matrix_abud_imgl <- as.matrix(prematrix_abud_imgl)
rownames(matrix_abud_imgl) <- gene_names_abud
sample_names_abudimgl <- colnames(imgl_abud)[2:10]

#   MICROGLIA/MONOCYTE ABUD DATA
#-------------------------------------------------------------------------------
objectAbud <- ReadData(Data = matrix_abud,
                           Labels = c(sample_names_abudmm),
                           Platform = NULL,
                           verbose = FALSE)

resultsAbud <- predict_RF(classifier = RF_classifier,
                              Data = objectAbud,
                              impute = TRUE) # can handle missed genes by imputation

pred_abud <- resultsAbud$predictions

# # if the classifier trained using probability   = TRUE
if (is.matrix(pred_abud)) {
  x3abud <- colnames(pred_abud)[max.col(pred_abud)]
}


predlonger_abud <- reshape2::melt(pred_abud)%>%
  mutate(fixedValues = round(value, 3))

#use barplot to visualize predictions by numbers
abud_plot <- ggplot(predlonger_abud, aes(x=Var2, y=value, fill=Var1), color="black") +
  labs(title = "Depmap Cell Prediction (Abud Data)", x="Cell Type", y="Cell Prediction Score") +
  geom_col(position="dodge") + theme_bw() +
  geom_text(aes(label = fixedValues), vjust=-0.5, position = position_dodge(.9), size=1.79) +
  scale_x_discrete(labels= c("astrocyte"="Astrocyte",
                             "microglia" = "Microglia",
                             "Whole.cortex" = "Whole cortex",
                             "oligodendrocyte" = "Oligodendrocyte")) +
  ylim(0,1) +
  theme(plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = "right", legend.text=element_text(size=6.5)) +
  scale_fill_viridis_d(labels = sample_names_abudmm)


pdf("Figures/Depmap_abud_plot.pdf",height=6,width=6,useDingbats=F)
print(abud_plot)
dev.off()

#   iMGL ABUD DATA
#-------------------------------------------------------------------------------
objectAbudimgl <- ReadData(Data = matrix_abud_imgl,
                       Labels = c(sample_names_abudimgl),
                       Platform = NULL,
                       verbose = FALSE)

resultsAbudimgl <- predict_RF(classifier = RF_classifier,
                          Data = objectAbudimgl,
                          impute = TRUE) # can handle missed genes by imputation

pred_abud_imgl <- resultsAbudimgl$predictions

# # if the classifier trained using probability   = TRUE
if (is.matrix(pred_abud_imgl)) {
  x3abud_imgl <- colnames(pred_abud_imgl)[max.col(pred_abud_imgl)]
}


predlonger_abudimgl <- reshape2::melt(pred_abud_imgl)%>%
  mutate(fixedValues = round(value, 3))

#use barplot to visualize predictions by numbers
abudiplot <- ggplot(predlonger_abudimgl, aes(x=Var2, y=value, fill=Var1), color="black") +
  labs(title = "iMGL Depmap Cell Prediction (Abud Data)", x="Cell Type", y="Cell Prediction Score") +
  geom_col(position="dodge") + theme_bw() +
  geom_text(aes(label = fixedValues), vjust=-0.5, position = position_dodge(.9), size=1.79) +
  scale_x_discrete(labels= c("astrocyte"="Astrocyte",
                             "microglia" = "Microglia",
                             "Whole.cortex" = "Whole cortex",
                             "oligodendrocyte" = "Oligodendrocyte")) +
  ylim(0,1) +
  theme(plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = "right", legend.text=element_text(size=6.5)) +
  scale_fill_viridis_d(labels = sample_names_abudimgl)

pdf("Figures/Depmap_budi_plot.pdf",height=6,width=6,useDingbats=F)
print(abudiplot)
dev.off()

#-------------------------------------------------------------------------------
#                     U87 samples from Gupta data
#-------------------------------------------------------------------------------
#DATA = dfGupta

#get gene names from main dataset
gene_names_gupta <- pull(dfGupta, "Gene")

#create matrix
prematrix_gupta <- select(dfGupta, -"Gene")
matrix_gupta <- as.matrix(prematrix_gupta)
rownames(matrix_gupta) <- gene_names_gupta

#   U87 GUPTA DATA
#-------------------------------------------------------------------------------
objectGupta <- ReadData(Data = matrix_gupta,
                        Labels = c("U87", "U87", "U87", "U87"),
                        Platform = NULL,
                        verbose = FALSE)

resultsGupta <- predict_RF(classifier = RF_classifier,
                           Data = objectGupta,
                           impute = TRUE) # can handle missed genes by imputation

pred_gupta <- resultsGupta$predictions

# # if the classifier trained using probability   = TRUE
if (is.matrix(pred_gupta)) {
  x3gupta <- colnames(pred_gupta)[max.col(pred_gupta)]
}


predlonger_gupta <- reshape2::melt(pred_gupta)%>%
  mutate(fixedValues = round(value, 3))

sample_names_gupta = colnames(dfGupta)[2:5]

#use barplot to visualize predictions by numbers
guptaplot <- ggplot(predlonger_gupta, aes(x=Var2, y=value, fill=Var1), color="black") +
  labs(title = "U87 Depmap Cell Prediction - Gupta Data", x="Cell Type", y="Cell Prediction Score") +
  geom_col(position="dodge") + theme_bw() +
  geom_text(aes(label = fixedValues), vjust=-0.5, position = position_dodge(.9), size=1.79) +
  scale_x_discrete(labels= c("astrocyte"="Astrocyte",
                             "microglia" = "Microglia / Macrophage",
                             "whole.cortex" = "Whole cortex",
                             "oligodendrocyte" = "Oligodendrocyte")) +
  ylim(0,1) +
  theme(plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = "right", legend.text=element_text(size=6.5)) +
  scale_fill_viridis_d(labels = sample_names_gupta)

pdf("Figures/Depmap_gupta_plot.pdf",height=6,width=6,useDingbats=F)
print(guptaplot)
dev.off()

#-------------------------------------------------------------------------------
#                     U87 samples from Chandra data
#-------------------------------------------------------------------------------
#DATA = dfChandra

#get gene names from main dataset
gene_names_chandra <- pull(dfChandra, "Gene")

#create matrix
prematrix_chandra <- select(dfChandra, -"Gene")
matrix_chandra <- as.matrix(prematrix_chandra)
rownames(matrix_chandra) <- gene_names_chandra

#   U87 DATA
#-------------------------------------------------------------------------------
objectChandra <- ReadData(Data = matrix_chandra,
                        Labels = c("U87", "U87"),
                        Platform = NULL,
                        verbose = FALSE)

resultsChandra <- predict_RF(classifier = RF_classifier,
                           Data = objectChandra,
                           impute = TRUE) # can handle missed genes by imputation

pred_chandra <- resultsChandra$predictions

predlonger_chandra <- reshape2::melt(pred_chandra)%>%
  mutate(fixedValues = round(value, 3))

sample_names_chandra = colnames(dfChandra)[2:5]

#use barplot to visualize predictions by numbers
chandraplot <- ggplot(predlonger_chandra, aes(x=Var2, y=value, fill=Var1), color="black") +
  labs(title = "U87 Depmap Cell Prediction (Chandra Data)", x="Cell Type", y="Cell Prediction Score") +
  geom_col(position="dodge") + theme_bw() +
  geom_text(aes(label = fixedValues), vjust=-0.5, position = position_dodge(.9), size=1.79) +
  scale_x_discrete(labels= c("astrocyte"="Astrocyte",
                             "microglia" = "Microglia / Macrophage",
                             "whole.cortex" = "Whole cortex",
                             "oligodendrocyte" = "Oligodendrocyte")) +
  ylim(0,1) +
  theme(plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = "right", legend.text=element_text(size=6.5)) +
  scale_fill_viridis_d(labels = sample_names_gupta)


pdf("Figures/Depmap_chandra_plot.pdf",height=6,width=6,useDingbats=F)
print(chandraplot)
dev.off()



#-------------------------------------------------------------------------------
#                     HMC3 samples from Armanville data
#-------------------------------------------------------------------------------
#read in data
#dfArman
#create the matrix 
#get gene names from main dataset
gene_names_arman <- pull(dfArman, "Gene")

#create HMC3 matrix
prematrix_arman <- select(dfArman, -"Gene")
matrix_arman <- as.matrix(prematrix_arman)
rownames(matrix_arman) <- gene_names_arman

#    HMC3 ARMAN DATA
#-------------------------------------------------------------------------------
objectHMC3arman <- ReadData(Data = matrix_arman,
                           Labels = c("HMC3", "HMC3", "HMC3"),
                           Platform = NULL,
                           verbose = FALSE)

resultsHMC3arman <- predict_RF(classifier = RF_classifier,
                              Data = objectHMC3arman,
                              impute = TRUE) # can handle missed genes by imputation

hmc3_pred_arman <- resultsHMC3arman$predictions

hmc3_predlonger_arman <- reshape2::melt(hmc3_pred_arman)%>%
  mutate(fixedValues = round(value, 3))


#use barplot to visualize predictions by numbers
armanville_pdf <- ggplot(hmc3_predlonger_arman, aes(x=Var2, y=value, fill=Var1), color="black") +
  labs(title = "HMC3 Cell Line Prediction (Armanville Data)", x="Cell Type", y="Cell Prediction Score") +
  geom_col(position="dodge") + theme_bw() +
  geom_text(aes(label = fixedValues), vjust=-0.5, position = position_dodge(.9), size=1.79) +
  scale_x_discrete(labels= c("astrocyte"="Astrocyte",
                             "microglia" = "Microglia / Macrophage",
                             "whole.cortex" = "Whole cortex",
                             "oligodendrocyte" = "Oligodendrocyte")) +
  ylim(0,1) +
  theme(plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = "right", legend.text=element_text(size=6.5)) +
  scale_fill_viridis_d(labels =c("Sample 1", "Sample 2", "Sample 3"))

pdf("Figures/Depmap_armanville.pdf",height=6,width=6)
print(armanville_pdf)
dev.off()


