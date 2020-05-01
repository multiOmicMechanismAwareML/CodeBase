library(reshape2)
library(plotly)
library(nsga2R)
library(tidyverse)
library(kernlab)
library(ranger)
library(ggplot2)
library(parallel)
library(lattice)
library(randomForest)
library(mlbench)
library(foreach)
library(GenAlgo)
library(iterators)
library(gridExtra)
library(caret)
library(stringr)
library(devtools)
library(Metrics)
library(iRF)
library(cluster)
library(doSNOW)
library(doMC)
library(purrr)
library(dplyr)
library(SGL)
registerDoMC()


set.seed(100)

#Function used to set up paralell computing 
resetCluster <- function(){
  stopCluster(cl)
  cores = detectCores(all.tests = FALSE, logical = TRUE)
  cl<-makeCluster(cores-1, outfile="Log.txt")
  registerDoSNOW(cl)
}

cores = detectCores(all.tests = FALSE, logical = TRUE)
cl<-makeCluster(cores-1, outfile="Log.txt")
registerDoSNOW(cl)


#Reading data sources 
data = read.csv('completeDataset.csv')
expressionData = read.csv('expressionOnly.csv')
metabolicExpression = read.csv('metabolic_gene_data.csv'
                               )
biomassGrowthFlux = 'r_4041';
biomassIndex = 717;
growth = data$log2relT  #Extract the growth rates
strands = data$Row #Extract the strand names 
data$Row <- NULL
data$log2relT <- NULL #
expressionData[,1] <- NULL #Remove the growth rate data
data <- data[-nearZeroVar(data)] #Here we remove any variables where the variance is close to zero

fluxData <- data[, -which(colnames(data) %in% colnames(expressionData))] #flux data is all data without gene data


#########################################################################################
#Data partition #########################################################################
#########################################################################################

#Full data partitions
trainingIndices<- createDataPartition(y = growth, p = 0.8, list = FALSE)
testingTarget <- growth[-trainingIndices]
trainingTarget <- growth[trainingIndices]
trainingData <- data[trainingIndices,]
testingData <-  data[-trainingIndices,]
normParam <- preProcess(trainingData)
normalisedTestData <- predict(normParam, testingData)
normalisedTrainingData <- predict(normParam, trainingData)

#Expression data partition 
expressionTrainingData <- expressionData[trainingIndices,]
expressionTestingData <- expressionData[-trainingIndices,]
normParamExp <- preProcess(expressionTrainingData)
normalisedExpressionTrainingData <- predict(normParamExp, expressionTrainingData)
normalisedExpressionTestingData <- predict(normParamExp, expressionTestingData)


#Flux data partition 
fluxTrainingData <- fluxData[trainingIndices,]
fluxTestingData <- fluxData[-trainingIndices,]
normParamFlux <- preProcess(fluxTrainingData)
normalisedFluxTrain <- predict(normParamFlux, fluxTrainingData)
normalisedFluxTest <- predict(normParamFlux, fluxTestingData)

#Metabolic expression data partition 

metabolicExpressionTrainingData <- metabolicExpression[trainingIndices,]
metabolicExpressionTestingData <- metabolicExpression[-trainingIndices,]
normParamMetabolicExpression <- preProcess(metabolicExpressionTrainingData)
normalisedMetabolicExpressionTrain <- predict(normParamMetabolicExpression, metabolicExpressionTrainingData)
normalisedMetabolicExpressionTest <- predict(normParamMetabolicExpression, metabolicExpressionTestingData)

#########################################################################################
#Feature Selection ######################################################################
#########################################################################################

#scale and transpose the data (so that we are clustering by feature and not by strand)
dataClust <- t(scale(trainingDataNum))

#pearsons distance matrix, used in the silhouette 
distanceMatrix <- get_dist(dataClust, method = "pearson")

#kmeans clustering applied using the distance matrix silhouette 
avg_sil <- function(k){
  km.result <- kmeans(dataClust, centers = k, nstart = 25)
  ss <- silhouette(km.result$cluster, distanceMatrix)
  mean(ss[,3])
}

#Values explored for k
k.values <- seq(250, 400, 25)

avg.sil.values <- map_dbl(k.values, avg_sil)

plot(k.values, avg.sil.values ,
           type = "b", pch = 19, frame = FALSE, 
           xlab = "Number of clusters K",
           ylab = "Average Silhouettes")



#################################################################################################################
#### Genetic Algorithms #########################################################################################



#multi objective #########################################################################################

#One of the objectives is to minimise the k nearest neighbour error score
min_rmse <- function(x){
  modelRes <- train(normalisedTrainingData[,which(x >= 0.99)],   #The x >= 0.99 only takes features with >0.99 weighting 
                    trainingTarget,  method = "knn" )
  print(modelRes$results$RMSE)
  mean(modelRes$results$RMSE)
}


ObjFun <- function (x){
  f1 <- min_rmse(x)
  f2 <- length(which(x >= 0.99))  #Second objective is to minimise the number of features used 
  c(f1=f1, f2=f2)
}

no.features <- ncol(normalisedTrainingData)
Solution <- nsga2(ObjFun, no.features, 2 , lower.bounds = rep(0, no.features), upper.bounds = rep(1, no.features),popsize = 80, generations = 250)


################################################################################################################
### Sparse group Lasso ########################################################################################

cvFit = cvSGL(trainingData , kmeans(dataClust, centers = which.min(avg.sil.values), nstart = 25)$cluster, type = "linear")

#Extract features 
sglFeatureIndicies <- which(colnames(trainingData) %in% namesOfFeaturesExtractedFromSparseGroup)

################################################################################################################
### Iterative random forest ####################################################################################

iRFfullData <- iRF(x= trainingDataNum, 
          y= trainingTarget, 
          xtest= testingDataNum, 
          ytest= testingTarget, 
          n.iter=5, 
          verbose = T,
          n.core = 7,
          interactions.return=5,
          n.bootstrap=10
)

f <- function (x)
 (substr(as.character(x), 1, 7))


#The results are returnes in concat string format so needed some playing around with to extract particular feature sets
iRFResults <- attributes(iRFfullData$interaction[[5]])
iRFResults <- str_replace(iRFResults, "_A", "*A") #We do this so that when we split by '_' we dont loose the _A 
iRFResults <- str_replace(iRFResults, "_B", "*B")
iRFResults <- unlist(strsplit(iRFResults, "_")) 
iRFResults <- unique(str_replace(iRFResults, "\\*", "_"))  #turn back to correct format

iRFColumnIndicies <- which(colnames(trainingData) %in% iRFResults)

