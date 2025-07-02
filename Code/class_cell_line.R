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
library(reshape2)
library(RColorBrewer)
library(Biobase)
library(switchBox)

#NOTE: conda install for M3 or M4 Macs doesn't work well. Best to use a linux 
#system



#===============================================================================
##------------------READ IN DATA------------------------------------------------
#===============================================================================

#Grab the filenames
classLabels = "Processed_data/classlabels_line.csv"
platformLabels = "Processed_data/platformlabels_line.csv"
comboDf_fname = "Processed_data/cleandata_combined_line.csv"
comboDf_fname_MHB = "Processed_data/cleandata_combined_training_line.MHB.csv"

#tests and validations 
hmc3_df_fname = "Processed_data/cleandata_Hmc3_adab.csv"
adabDf_fname = "Processed_data/cleandata_Adab.csv"
labels_adab_fname = "Data/classlabels_adab.v6.csv"
imglDf_fname = "Processed_data/cleandata_imgl_adab.csv"

hmc3_chai = read_csv("Processed_data/cleandata_Hmc3_chai.csv")
astro_setRow = read_csv("Processed_data/dfRow_astrocyteonly.csv")
dfSagar = read_csv("Processed_data/dfSagar_microglia.csv")
dfGupta = read_csv("Processed_data/cleandata_Gupta.csv")
dfGoyal = read_csv("Processed_data/cleandata_Goyal.csv")
dfArman = read_csv("Processed_data/cleandata_Hmc3_arman.csv")



################################################################################
##------------------LOAD IN DATA------------------------------------------------
################################################################################
allClassLabels <- read_csv(classLabels, col_names = FALSE)
allPlatformLabels <- read_csv(platformLabels, col_names = FALSE)
#fullDf <- read_csv(comboDf_fname)
fullDf <- read_csv(comboDf_fname_MHB)

#Create class and platform labels
L1 <- pull(allClassLabels, X1)
P1 <- pull(allPlatformLabels, X1)

#verify labels and df dimensions are all correct
dimsFull <- dim(fullDf)
print(paste("The full dataframe has", dimsFull[2], "columns. The labels should have",
            dimsFull[2]-1, "length. The lengths of the class labels and platforms labels are",
            length(L1),"and", length(P1)))

#make a matrix of the data for the multiclasspairs 
genes = pull(fullDf, Gene)
matrixPrepDf <- fullDf[,-1] #remove non-numeric gene column
matrixPrepDf <- sapply(matrixPrepDf, as.numeric)
#head(matrixPrepDf)

fullMatrix <- as.matrix(matrixPrepDf)
rownames(fullMatrix) <- genes #add gene names as rows of matrix

#fullMatrix #now the final matrix with all data combined 



#-------------------------------------------------------------------------------
##---------------------- MULTICLASS PAIRS TRAINING -----------------------------
#-------------------------------------------------------------------------------

## create the data object
object <- ReadData(Data = fullMatrix,
                   Labels = L1,
                   Platform = P1,
                   verbose = FALSE)

#NOTE: Look at what this is and how to store it or save it for writing purposes. 
object


#------------------------------------------------------
#------------------START TRAINING DATA-----------------
#------------------------------------------------------

## split the data
# 70% as training data and 30% as testing data

#get dims of matrix
n <- ncol(fullMatrix)
nr <- nrow(fullMatrix)

set.seed(123)
training_samples <- sample(1:n, size = n*.7) #randomly pulls 60% of indices for training data

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
#NOTE: Look at what this is and how to store it or save it for writing purposes.
objectTrain

# -------------------------------------------------------------------------
#TEST OBJECT

objectTest <- ReadData(Data = test,
                       Labels = L1.test,
                       Platform = P1.test,
                       verbose = FALSE)
#NOTE: Look at what this is and how to store it or save it for writing purposes.
objectTest


#------------------------------------------------------
#----------------------RANDOM FOREST SCHEME------------
#------------------------------------------------------

# sorting genes
genes_RF <- sort_genes_RF(data_object = objectTrain,
                          featureNo_altogether = 50,
                          featureNo_one_vs_rest = 50,
                          rank_data = TRUE,
                          platform_wise = FALSE,
                          num.trees = 500,
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
                                  genes_altogether = c(25, 50, 75, 85,100,150),
                                  genes_one_vs_rest = c(25, 50, 75, 85,100,150))
knitr::kable(summary_genes)

# to give enough number of  rules and unique genes for our classes
# Now let's run sort_rules_RF to create the rules and sort them
rules_RF <- sort_rules_RF(data_object = objectTrain,
                          sorted_genes_RF = genes_RF,
                          genes_altogether = 50,
                          genes_one_vs_rest = 50,
                          num.trees = 1000,
                          seed=123456,
                          verbose = TRUE)



rules_RF # sorted rules object




#------------------MODEL TRAINING----------------------

# prepare the simple data.frame for the parameters I want to test
# names of arguments as column names
# this df has three sets (3 rows) of parameters
parameters <- data.frame(
  gene_repetition=c(1,1,1,1,1),
  rules_one_vs_rest=c(25, 50, 50, 85, 50),
  rules_altogether=c(50, 25, 50, 85, 50),
  run_boruta=c(TRUE, TRUE,TRUE, TRUE, TRUE),
  plot_boruta = FALSE,
  num.trees=c(500,500, 1000, 1000,500),
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
                          rules_altogether = 50,
                          rules_one_vs_rest = 50,
                          run_boruta = TRUE,
                          plot_boruta = FALSE,
                          probability = TRUE,
                          num.trees = 1000,
                          boruta_args = list(),
                          verbose = TRUE)

#-------------------model TRAINING---------------------------
# plot proximity matrix of the out-of-bag samples
# Note: this takes a lot of time if the data is big
pdf("Figures/PredAccur_cell_training.pdf",heigh=10, width=10,useDingbats=F)
proximity_matrix_RF(object = objectTrain,
                    classifier = RF_classifier,
                    plot = TRUE,
                    return_matrix = FALSE, # if we need to get the matrix itself
                    title = "Prediction Accuracy on Training Data - Cell Line Data",
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

pdf("Figures/PredAccur_cell_testing.pdf",height=10,width=10,useDingbats=F)
proximity_matrix_RF(object = objectTest,
                    classifier = RF_classifier,
                    plot = TRUE,
                    return_matrix = FALSE, # if we need to get the matrix itself
                    title = "Prediction Accuracy on Testing Data - Cell Line Data",
                    cluster_cols = FALSE)
dev.off()

#------------------VISUALIZATION--------------
#visualize the binary rules in training dataset
pdf("Figures/PredAccur_cell_training_rules.pdf",height=10,width=10,useDingbats=F)
plot_binary_RF(Data = objectTrain,
               classifier = RF_classifier,
               prediction = NULL,
               as_training = TRUE, # to extract the scores from the model
               show_scores = TRUE,
               top_anno = "ref",
               show_predictions = TRUE,
               cluster_cols = TRUE, #otherwise there's an issue with clustering too small
               #margin = c(0,5,0,8),
               title = "Training data - Cell Line Data")
dev.off()

# visualiAze the binary rules in testing dataset
pdf("Figures/PredAccur_cell_training_rules.pdf",height=10,width=10,useDingbats=F)
plot_binary_RF(Data = objectTest,
               ref = L1.test, # Get ref labels from the test data
               classifier = RF_classifier,
               prediction = results,
               as_training = FALSE,
               show_scores = TRUE,
               top_anno = "ref",
               show_predictions = TRUE,
               title = "Testing data - Cell Line Data",
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


hmc3_predlonger_chai <- melt(hmc3_pred_chai)%>%
  mutate(fixedValues = round(value, 3))


#use barplot to visualize predictions by numbers
chai_pdf<- ggplot(hmc3_predlonger_chai, aes(x=Var2, y=value, fill=Var1), color="black") +
  labs(title = "HMC3 Cell Line Prediction (Chai Data)", x="Cell Type", y="Cell Prediction Score") +
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

pdf("Figures/CL_chai.pdf",height=6,width=6)
print(chai_pdf)
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

hmc3_predlonger_arman <- melt(hmc3_pred_arman)%>%
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

pdf("Figures/CL_armanville.pdf",height=6,width=6)
print(armanville_pdf)
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

hmc3_predlonger <- melt(hmc3_pred)%>%
  mutate(fixedValues = round(value, 3))


#use barplot to visualize predictions by numbers
quiroga_pdf <- ggplot(hmc3_predlonger, aes(x=Var2, y=value, fill=Var1), color="black") +
  labs(title = "HMC3 Cell Line Prediction - Quiroga Data", x="Cell Type", y="Cell Prediction Score") +
  geom_col(position="dodge") + theme_bw() +
  geom_text(aes(label = fixedValues), vjust=-0.5, position = position_dodge(.9), size=1.79) +
  scale_x_discrete(labels= c("astrocyte"="Astrocyte",
                             "microglia" = "Microglia / Macrophage",
                             "whole.cortex" = "Whole cortex",
                             "oligodendrocyte" = "Oligodendrocyte")) +
  ylim(0,1) +
  theme(plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = c(.72,.7), legend.text=element_text(size=6.5)) +
  scale_fill_viridis_d(labels =c("HMC3 ctl 1", "HMC3 1uM AB 1", "HMC3 ctl 2", "HMC3 1uM AB 2"))


pdf("Figures/CL_quiroga.pdf",height=6,width=6)
print(quiroga_pdf)
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

imgl_predlonger <- melt(imgl_pred)%>%
  mutate(fixedValues = round(value, 3))


#use barplot to visualize predictions by numbers
imgl_quiroga_pdf <- ggplot(imgl_predlonger, aes(x=Var2, y=value, fill=Var1), color="black") +
  labs(title = "iMGL Cell Prediction - Quiroga", x="Cell Type", y="Cell Prediction Score") +
  geom_col(position="dodge") + theme_bw() +
  geom_text(aes(label = fixedValues), vjust=-0.5, position = position_dodge(.9), size=1.79) +
  scale_x_discrete(labels= c("astrocyte"="Astrocyte",
                             "microglia" = "Microglia / Macrophage",
                             "whole.cortex" = "Whole cortex",
                             "oligodendrocyte" = "Oligodendrocyte")) +
  ylim(0,1) +
  theme(plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = "right", legend.text=element_text(size=6.5)) +
  scale_fill_viridis_d(labels =c("Sample 1","Sample 2","Sample 3","Sample 4","Sample 5",
                                 "Sample 6","Sample 7","Sample 8"))

pdf("Figures/CL_quiroga_imgl.pdf",height=6,width=6)
print(imgl_quiroga_pdf)
dev.off()



#    ADAB U87 AND THP1 DATA
#-------------------------------------------------------------------------------
#create the matrix for dfADAB

adab_matrix <- as.matrix(dfADAB[,-1]) #doesn't include the Gene column 
rownames(adab_matrix) <- dfADAB$Gene #have gene names be the row names 

labels_adab <- read_csv(labels_adab_fname, col_names=FALSE)
adab_labels <- pull(labels_adab, X1) #pull the labels into a vector instead 

#create the adab object 
objectADAB <- ReadData(Data = adab_matrix,
                       Labels = adab_labels,
                       Platform = NULL,
                       verbose = FALSE)

#results
resultsADAB <- predict_RF(classifier = RF_classifier,
                          Data = objectADAB,
                          impute = TRUE) # can handle missed genes by imputation

adab_pred <- resultsADAB$predictions

# # if the classifier trained using probability   = TRUE
if (is.matrix(adab_pred)) {
  x3 <- colnames(adab_pred)[max.col(adab_pred)]
}

adab_predlonger <- melt(adab_pred)%>%
  mutate(fixedValues = round(value, 3))

#use barplot to visualize predictions by numbers
adab_pdf <- ggplot(adab_predlonger, aes(x=Var2, y=value, fill=Var1), color="black") +
  labs(title = "Cell Prediction - Quiroga Cell Line Data", x="Cell Type", y="Cell Prediction Score") +
  geom_col(position="dodge") + theme_bw() +
  geom_text(aes(label = fixedValues), vjust=-0.5, position = position_dodge(.9), size=1.79) +
  scale_x_discrete(labels= c("astrocyte"="Astrocyte",
                             "microglia" = "Microglia / Macrophage",
                             "whole.cortex" = "Whole cortex",
                             "oligodendrocyte" = "Oligodendrocyte")) +
  ylim(0,1) +
  theme(plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = "right",
        legend.key.size = unit(0.5, "cm"), legend.text=element_text(size=6.5)) +
  scale_fill_viridis_d(labels = adab_labels)


pdf("Figures/CT_adab.pdf",height=6,width=6)
print(adab_pdf)
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
                               impute_reject= 0.99)
#can handle missed genes by imputation

astroRow_pred <- resultsRow_astro$predictions

# # if the classifier trained using probability   = TRUE
if (is.matrix(astroRow_pred)) {
  x3 <- colnames(astroRow_pred)[max.col(astroRow_pred)]
}

astroRow_predlonger <- melt(astroRow_pred)%>%
  mutate(fixedValuesRowAstro = round(value, 3))

#use barplot to visualize predictions by numbers
astro_pdf <- ggplot(astroRow_predlonger, aes(x=Var2, y=value, fill=Var1), color="black") +
  labs(title = "Schirmer Data Sample Prediction", x="Cell Type", y="Cell Prediction Score") +
  geom_col(position="dodge") + theme_bw() +
  geom_text(aes(label = fixedValuesRowAstro), vjust=-0.5, position = position_dodge(.9), size=.5) +
  scale_x_discrete(labels= c("astrocyte"="Astrocyte")) +
  theme(plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position.inside = c(.72,.7), legend.text=element_text(size=6.5)) +
  scale_fill_viridis_d(labels =c("C7Astroctye", "C9Astroctye", "C2Astroctye", "C8Astroctye",
                                 "C4Astroctye", "C5Astroctye", "C1Astroctye", "C6Astroctye", "C3Astroctye"))


pdf("Figures/CT_astro.pdf",height=6,width=6)
print(astro_pdf)
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
                           impute_reject = 0.99) # can handle missed genes by imputation

sagar_pred <- resultsSagar$predictions

# # if the classifier trained using probability   = TRUE
if (is.matrix(sagar_pred)) {
  x3 <- colnames(sagar_pred)[max.col(sagar_pred)]
}

sagar_predlonger <- melt(sagar_pred)%>%
  mutate(fixedValuesSag = round(value, 3))

#use barplot to visualize predictions by numbers
masuda_pdf <- ggplot(sagar_predlonger, aes(x=Var2, y=value, fill=Var1), color="black") +
  labs(title = "Masuda Data Sample Prediction", x="Cell Type", y="Cell Prediction Score") +
  geom_col(position="dodge") + theme_bw() +
  geom_text(aes(label = fixedValuesSag), vjust=-0.5, position = position_dodge(0.9), size=.5) +
  scale_x_discrete(labels= c("astrocyte"="Astrocyte",
                             "microglia" = "Microglia / Macrophage",
                             "whole.cortex" = "Whole cortex",
                             "oligodendrocyte" = "Oligodendrocyte")) +
  theme(plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = c(.72,.7), legend.text=element_text(size=6.5)) +
  scale_fill_viridis_d(labels =c("Pat3.Myeloid", "Pat1.Myeloid", "Pat2.Myeloid", "Pat5.Myeloid",
                                 "MS2.Myeloid", "MS4.Myeloid", "MS5.Myeloid"))


pdf("Figures/CT_masuda.pdf",height=6,width=6)
print(masuda_pdf)
dev.off()

#-------------------------------------------------------------------------------
#                     A172 samples from Goyal data
#-------------------------------------------------------------------------------
#DATA = dfGoyal

#get gene names from main dataset
gene_names_goyal <- pull(dfGoyal, "Gene")

#create matrix
prematrix_goyal <- select(dfGoyal, -"Gene")
matrix_goyal <- as.matrix(prematrix_goyal)
rownames(matrix_goyal) <- gene_names_goyal

#   A172 DATA
#-------------------------------------------------------------------------------
objectGoyal <- ReadData(Data = matrix_goyal,
                        Labels = c("U87", "U87", "U87"),
                        Platform = NULL,
                        verbose = FALSE)

resultsGoyal <- predict_RF(classifier = RF_classifier,
                           Data = objectGoyal,
                           impute = TRUE) # can handle missed genes by imputation

pred_goyal <- resultsGoyal$predictions

predlonger_goyal <- melt(pred_goyal)%>%
  mutate(fixedValues = round(value, 3))

sample_names_goyal = colnames(dfGoyal)[2:4]

#use barplot to visualize predictions by numbers
goyal_pdf = ggplot(predlonger_goyal, aes(x=Var2, y=value, fill=Var1), color="black") +
  labs(title = "A172 Cell Line Prediction (Goyal Data)", x="Cell Type", y="Cell Prediction Score") +
  geom_col(position="dodge") + theme_bw() +
  geom_text(aes(label = fixedValues), vjust=-0.5, position = position_dodge(.9), size=1.79) +
  scale_x_discrete(labels= c("astrocyte"="Astrocyte",
                             "microglia" = "Microglia / Macrophage",
                             "whole.cortex" = "Whole cortex",
                             "oligodendrocyte" = "Oligodendrocyte")) +
  ylim(0,1) +
  theme(plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = "right", legend.text=element_text(size=6.5)) +
  scale_fill_viridis_d(labels = sample_names_goyal)

pdf("Figures/CL_a172.pdf",height=6,width=6)
print(goyal_pdf)
dev.off()


#-------------------------------------------------------------------------------
#                     U87 samples from Gupta data - not used in this run rn
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
  x3chai <- colnames(pred_gupta)[max.col(pred_gupta)]
}


predlonger_gupta <- melt(pred_gupta)%>%
  mutate(fixedValues = round(value, 3))

sample_names_gupta = colnames(dfGupta)[2:5]

#use barplot to visualize predictions by numbers
gupta_pdf <- ggplot(predlonger_gupta, aes(x=Var2, y=value, fill=Var1), color="black") +
  labs(title = "U87 Cell Line Prediction - Gupta Data", x="Cell Type", y="Cell Prediction Score") +
  geom_col(position="dodge") + theme_bw() +
  geom_text(aes(label = fixedValues), vjust=-0.5, position = position_dodge(.9), size=1.79) +
  scale_x_discrete(labels= c("astrocyte"="Astrocyte",
                             "microglia" = "Microglia / Macrophage",
                             "whole.cortex" = "Whole cortex",
                             "oligodendrocyte" = "Oligodendrocyte")) +
  ylim(0,1) +
  theme(plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = c(.72,.7), legend.text=element_text(size=6.5)) +
  scale_fill_viridis_d(labels = sample_names_gupta)


pdf("Figures/CL_u87.pdf",height=6,width=6)
print(gupta_pdf)
dev.off()

