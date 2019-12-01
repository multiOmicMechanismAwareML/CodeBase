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
#### Kernel Matrix ##############################################################################################
source("bemkl.R")

exp_train <- as.matrix(normalisedMetabolicExpressionTrain)
reac_train <- as.matrix(normalisedFluxTrain)
colnames(exp_train) <- NULL
colnames(reac_train) <- NULL
rbf_exp_1 <- kernelMatrix(rbfdot(sigma = 0.1), exp_train)
rbf_exp_2 <- kernelMatrix(rbfdot(sigma = 0.01), exp_train)
rbf_exp_3 <- kernelMatrix(rbfdot(sigma = 0.001), exp_train)
rbf_exp_4 <- kernelMatrix(rbfdot(sigma = 0.0001), exp_train)
rbf_exp_5 <- kernelMatrix(rbfdot(sigma = 0.00001), exp_train)
rbf_exp_6 <- kernelMatrix(rbfdot(sigma = 0.000001), exp_train)
rbf_reac_1 <- kernelMatrix(rbfdot(sigma = 0.1), reac_train)
rbf_reac_2 <- kernelMatrix(rbfdot(sigma = 0.01), reac_train)
rbf_reac_3 <- kernelMatrix(rbfdot(sigma = 0.001), reac_train)
rbf_reac_4 <- kernelMatrix(rbfdot(sigma = 0.0001), reac_train)
rbf_reac_5 <- kernelMatrix(rbfdot(sigma = 0.00001), reac_train)
rbf_reac_6 <- kernelMatrix(rbfdot(sigma = 0.000001), reac_train)


exp_test <- as.matrix(normalisedMetabolicExpressionTest)
exp_reac <- as.matrix(normalisedFluxTest)
rbf_exp_1t <- kernelMatrix(rbfdot(sigma = 0.1), exp_train, y = exp_test)
rbf_exp_2t <- kernelMatrix(rbfdot(sigma = 0.01), exp_train, y = exp_test)
rbf_exp_3t <- kernelMatrix(rbfdot(sigma = 0.001), exp_train, y = exp_test)
rbf_exp_4t <- kernelMatrix(rbfdot(sigma = 0.0001), exp_train, y = exp_test)
rbf_exp_5t <- kernelMatrix(rbfdot(sigma = 0.00001), exp_train, y = exp_test)
rbf_exp_6t <- kernelMatrix(rbfdot(sigma = 0.000001), exp_train, y = exp_test)
rbf_reac_1t <- kernelMatrix(rbfdot(sigma = 0.1), reac_train, y = exp_reac)
rbf_reac_2t <- kernelMatrix(rbfdot(sigma = 0.01), reac_train, y = exp_reac)
rbf_reac_3t <- kernelMatrix(rbfdot(sigma = 0.001), reac_train, y = exp_reac)
rbf_reac_4t <- kernelMatrix(rbfdot(sigma = 0.0001), reac_train, y = exp_reac)
rbf_reac_5t <- kernelMatrix(rbfdot(sigma = 0.00001), reac_train, y = exp_reac)
rbf_reac_6t <- kernelMatrix(rbfdot(sigma = 0.000001), reac_train, y = exp_reac)


#initalize the parameters of the algorithm
parameters <- list()

#set the hyperparameters of gamma prior used for sample weights
parameters$alpha_lambda <- 1e-10
parameters$beta_lambda <- 1e+10

#set the hyperparameters of gamma prior used for intermediate noise
parameters$alpha_upsilon <- 1e-10
parameters$beta_upsilon <- 1e+10

#set the hyperparameters of gamma prior used for bias
parameters$alpha_gamma <- 1e-10
parameters$beta_gamma <- 1e+10

#set the hyperparameters of gamma prior used for kernel weights
parameters$alpha_omega <- 1e-10
parameters$beta_omega <- 1e+10

#set the hyperparameters of gamma prior used for output noise
parameters$alpha_epsilon <- 1e-10
parameters$beta_epsilon <- 1e+10

### IMPORTANT ###
#For gamma priors, you can experiment with three different (alpha, beta) values
#(1, 1) => default priors
#(1e-10, 1e+10) => good for obtaining sparsity
#(1e-10, 1e-10) => good for small sample size problems

#set the number of iterations
parameters$iteration <- 100000

#determine whether you want to store the lower bound values
parameters$progress <- 0

#set the seed for random number generator used to initalize random variables
parameters$seed <- 1606
all_kernels <- c(rbf_exp_1, rbf_exp_2, rbf_exp_3, rbf_exp_4, rbf_exp_5, rbf_exp_6, rbf_reac_1, rbf_reac_2, rbf_reac_3, rbf_reac_4, rbf_reac_5, rbf_reac_6)  
all_test_kernels <- c(rbf_exp_1t, rbf_exp_2t, rbf_exp_3t, rbf_exp_4t, rbf_exp_5t, rbf_exp_6t, rbf_reac_1t, rbf_reac_2t, rbf_reac_3t, rbf_reac_4t, rbf_reac_5t, rbf_reac_6t)  
Ktrain <- array(all_kernels, c(dim(rbf_exp_1),12))
Ktest <- array(all_test_kernels, c(dim(rbf_exp_1t), 12))
ytrain <- trainingTarget
  
#perform training
state <- bemkl_supervised_regression_variational_train(Ktrain, ytrain, parameters)

#display the kernel weights
print(state$be$mu[-1])

#perform prediction
prediction <- bemkl_supervised_regression_variational_test(Ktest, state)
RMSE <- postResample(prediction$y$mu, testingTarget)["RMSE"]



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


#ploting results 

sols <- as.data.frame(cbind(Solution$value, Solution$pareto.optimal))
colnames(sols) <- c("RMSE", "No.Features", "Pareto front")
ggplot(data = sols, aes(x = sols$RMSE, y = sols$No.Features, color = as.factor(sols$`Pareto front`) )) +
  geom_point(size = 2, alpha = 0.6) +
  labs(x = "RMSE", y = "#Features", color = "Pareto Front") + 
  theme_minimal() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                        text = element_text(size = 14),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="none")
justPareto <- sols[which(sols$`Pareto front` == 1),]
ggplot(data = justPareto , aes(x = justPareto$RMSE, y = justPareto$No.Features)) + labs(x = "RMSE", y = "#Features") + geom_line()




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


#########################################################################################
#Machine Learning #######################################################################
#########################################################################################

#Creating control setting 10 fold cross validation repeated 3 times 
ctr <- trainControl(method = "repeatedcv", number = 5, repeats =  3, verboseIter = FALSE)


########################################################################################
#SVM Gaussian ##########################################################################

#Param exploration 
svmTune <- expand.grid(C = c(10,50,100,250,500,1000, 2000) , sigma = c(0.00001, 0.0001, 0.001, 0.1, 1) )

#Flux 

#Sparse group lasso features dataset, structure as above
Fluxsvm <- train(normalisedFluxTrain, 
                trainingTarget,
                method = 'svmRadial',
                trControl=ctr,
                tuneGrid = svmTune)


FluxSvm <- svm(trainingTarget ~ ., data = normalisedFluxTrain, cost = Fluxsvm$bestTune$C, gamma =Fluxsvm$bestTune$sigma, kernel =
                "radial")
FluxSVMPred <-  predict(FluxSvm, normalisedFluxTest)
FluxSVMRMSE <- postResample(pred = FluxSVMPred, obs = testingTarget)

#Metabolic expression
#Sparse group lasso features dataset, structure as above
metExpsvm <- train(normalisedMetabolicExpressionTrain, 
                 trainingTarget,
                 method = 'svmRadial',
                 trControl=ctr,
                 tuneGrid = svmTune)


MetExpsvm <- svm(trainingTarget ~ ., data = normalisedMetabolicExpressionTrain, cost = metExpsvm$bestTune$C, gamma =metExpsvm$bestTune$sigma, kernel =
                 "radial")
MetExpSVMPred <-  predict(MetExpsvm, normalisedMetabolicExpressionTest)
MetExpSVMRMSE <- postResample(pred = MetExpSVMPred, obs = testingTarget)


#Train to learn parameters
svm_gaussianFull <- train(normalisedTrainingData, 
                    trainingTarget,
                    method = 'svmRadial',
                    trControl=ctr,
                    tuneGrid = svmTune)
#Train using the params from above 
FullDataSVM <- svm(trainingTarget ~ ., data = normalisedTrainingData, cost = svm_gaussianFull$bestTune$C, gamma = svm_gaussianFull$bestTune$sigma, kernel =
                "radial")

#Make predictions 
FullDataSVMPredictions <-  predict(FullDataSVM, normalisedTestData)
FullDataSVMRMSE <- postResample(pred = FullDataSVMPredictions, obs = testingTarget)
save(FullDataSVM, file = "fullDataSVMModel.RData")

#Train the expression only to learn params
svm_gaussianExpression <- train(normalisedExpressionTrainingData, 
                          trainingTarget,
                          method = 'svmRadial',
                          trControl=ctr,
                          tuneGrid = svmTune)

#Train using the params extracted from above 
ExpressionSVM <- svm(trainingTarget ~ ., data = normalisedExpressionTrainingData, cost = svm_gaussianExpression$bestTune$C, gamma = svm_gaussianExpression$bestTune$sigma, kernel =
                     "radial")
ExpressionSVMPred <-  predict(ExpressionSVM, normalisedExpressionTrainingData)
ExpressioniSVMRMSE <- postResample(pred = ExpressionSVMPred, obs = trainingTarget)
save(ExpressionSVM , file = "ExpressionSVMVModel.RData")

#Sparse group lasso features dataset, structure as above
SGLsvm <- train(normalisedTrainingData[,sglFeatureIndicies], 
                                    trainingTarget,
                                    method = 'svmRadial',
                                    trControl=ctr,
                                    tuneGrid = svmTune)


sglSVM <- svm(trainingTarget ~ ., data = normalisedTrainingData[,sglFeatureIndicies], cost = SGLsvm$bestTune$C, gamma =SGLsvm$bestTune$sigma, kernel =
                     "radial")
sglSVMPred <-  predict(sglSVM, normalisedTestData[,sglFeatureIndicies])
sglSVMRMSE <- postResample(pred = sglSVMPred, obs = testingTarget)


#Iterative Random Forest features-> structure as above
svm_iRFGaussian <- train(normalisedTrainingData[, iRFColumnIndicies], 
                            trainingTarget,
                            method = 'svmRadial',
                            trControl=ctr,
                            tuneGrid = svmTune   )

iRFSVM <- svm(trainingTarget ~ ., data = normalisedTrainingData[, iRFColumnIndicies], cost = svm_iRFGaussian$bestTune$C, gamma = svm_iRFGaussian$bestTune$sigma, kernel =
                "radial")
iRFSVMPredictions <-  predict(iRFSVM, normalisedTestData[, iRFColumnIndicies])
iRFRMSE <- postResample(pred = iRFSVMPredictions, obs = testingTarget)

#Genetic Algo same as above but this time we run through all of the models the gen algo highlighted
svm_gaussianGenAlgo <- list()
for (x in c(1 : nrow(justPareto))){
  print(x)
  indexes <- which(Solution$par[x,] > 0.99)
  gaussNext <- train(normalisedTrainingData[, indexes], trainingTarget,
                  method = 'svmRadial',
                  trControl = ctr,
                  verbose = F,
                  tuneGrid = svmTune
  )
  genSVM <- svm(trainingTarget ~ ., data = normalisedTrainingData[, indexes], cost = gaussNext$bestTune$C, gamma =gaussNext$bestTune$sigma, kernel =
                  "radial")
  gaussianPredict <- predict(genSVM, normalisedTestData[,indexes])
  gaussRMSE <- postResample(pred = gaussianPredict, obs = testingTarget)
  print(gaussRMSE)
  svm_gaussianGenAlgo[[x]] <- list(gaussRMSE, indexes, gaussianPredict, genSVM)
}



for (x in c(1 : nrow(justPareto))){
       print(x)
       indexes <- which(geneticAlgoSol$par[x,] > 0.99)
       t <- data[,indexes]
       write.table(t, file = paste("genfs", x, ".csv", sep=""))
   }
#########################################################################################
#Random Forest ##########################################################################

#Parameters explored in random forest
rfgrid <- expand.grid(.mtry=seq(2,500,20))

#Train against fluxes only 
rfFluxData <- train(normalisedFluxTrain, trainingTarget,
                          method = "rf",
                          trControl = ctr,
                          verbose = F,
                          tuneGrid = rfgrid
)
fluxPredictions <- predict(rfFluxData$finalModel, newdata = normalisedFluxTest)

fluxRMSE <- sqrt(mse(fluxPredictions, testingTarget ))


#train against metabolic genes only 
rfMetExpData <- train(normalisedMetabolicExpressionTrain, trainingTarget,
                    method = "rf",
                    trControl = ctr,
                    verbose = F,
                    tuneGrid = rfgrid
)
#Remember you need to rename from rfExpData Chris 
rfPredictions <- predict(rfMetExpData$finalModel, newdata = normalisedMetabolicExpressionTest)


metaExpressionOnlyRMSE <- sqrt(mse(rfPredictions, testingTarget ))

#late integration flux with metabolic expression
metExpWithFluxPredictions <- data.frame(flux = rfPredictions, rfMetPredictions)
integratedPrediction <- rowMeans(metExpWithFluxPredictions)
lateIntegratedFluxWithMetabolicExpressionRMSE <- sqrt(mse(integratedPrediction, testingTarget ))  

#late integration flux with full expression
metExpWithFluxPredictions <- data.frame(flux = fluxPredictions, expressionRfPredictions)
integratedPrediction <- rowMeans(metExpWithFluxPredictions)
lateIntegratedFluxWithMetabolicExpressionRMSE <- sqrt(mse(integratedPrediction, testingTarget ))  

#Full dataset 
rfFullData <- train(normalisedTrainingData, trainingTarget,
                     method = "rf",
                     trControl = ctr,
                     verbose = F,
                     tuneGrid = rfgrid
)
fulRMSE <- sqrt(mse(fullRfPredictions, testingTarget ))

#Train against expression only 
rfExpressionData <- train(normalisedExpressionTrainingData, trainingTarget,
                    method = "rf",
                    trControl = ctr,
                    verbose = F,
                    tuneGrid = rfgrid
)
expressionRfPredictions = predict(rfExpressionData$finalModel, newdata = normalisedExpressionTestingData)
expressionRMSE <- sqrt(mse(expressionRfPredictions, testingTarget ))


#Train using genetic algo Features 
rfGenAlgo <- list()
for (x in c(1 : nrow(justPareto))){
  print(x)
  indexes <- which(Solution$par[x,] > 0.99)
  rfNext <- train(normalisedTrainingData[, indexes], trainingTarget,
                  method = "rf",
                  trControl = ctr,
                  verbose = F,
                  tuneGrid = rfgrid
  )
  rfPredictions <- predict(rfNext$finalModel, newdata = normalisedTestData[,indexes])
  rfRMSE <- sqrt(mse(rfPredictions, testingTarget ))
  print(rfRMSE)
  rfGenAlgo[[x]] <- list(rfRMSE, indexes, rfPredictions, rfNext)
}


rffullplot <- ggplot(rfFullData) + ggtitle("Full Data Random Forest")
fullRfPredictions <- predict(rfFullData$finalModel, newdata = normalisedTestData)
expressionRfPredictions <- predict(rfExpressionData$finalModel, newdata =  normalisedTestData)

#Train using the iterative random forest dataset
rfiRFSelected <- train(normalisedTrainingData[, iRFColumnIndicies], trainingTarget,
                method = "rf",
                trControl = ctr,
                verbose = F,
                tuneGrid = rfgrid
)

rfiRFPredictions <- predict(rfiRFSelected$finalModel, newdata = normalisedTestData[,iRFColumnIndicies])
rfiRFRMSE <- sqrt(mse(rfiRFPredictions, testingTarget ))

#Train using the sparse group lasso dataset
rfsglSelected <- train(normalisedTrainingData[, sglFeatureIndicies], trainingTarget,
                       method = "rf",
                       trControl = ctr,
                       verbose = F,
                       tuneGrid = rfgrid
)

rfsglPredictions <- predict(rfsglSelected$finalModel, newdata = normalisedTestData[,sglFeatureIndicies])
rfsglRMSE <- sqrt(mse(rfsglPredictions, testingTarget ))



################################GRAPHING CODE##############################################

#Graph the biomass flux rates vs the experimental growth rates 
d <- data.frame(doubling_time = growth, biomass_flux = data[,biomassIndex] )
cor <- d %>% ggplot(aes( x = biomass_flux, y = doubling_time)) +
  geom_point(        color="blue4",
                     fill="white",
                     shape=21,
                     alpha=0.6,
                     size= 2,
                     stroke = 1) + 
  geom_smooth(method = lm, colour = 'blue4', linetype = "dashed", alpha = 0.5 , se=TRUE, size = 1) + 
  theme_minimal() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.ticks = element_line(),
        text = element_text(size=14),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="none") +
        labs(x = "Biomass Reaction Flux", y = "Experimental Doubling Time \n log(strain/w.t)")
