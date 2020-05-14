library(nsga2R)
library(doParallel)
library(MLmetrics)
library(ggplot2)
library(parallel)
library(randomForest)
library(GenAlgo)
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
library(e1071)
library(MLmetrics)
library("e1071", lib.loc="~/R/win-library/3.1")
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
#cl<-makeCluster(cores-1, outfile="Log.txt")
#registerDoSNOW(cl)


#Reading data sources 
data = readRDS('data/completeDataset.RDS')
data <- data[-nearZeroVar(data)] #Here we remove any variables where the variance is close to zero
expressionData = readRDS('data/expressionOnly.RDS')
metabolicExpression = readRDS('data/metabolic_gene_data.RDS')


biomassGrowthFlux = 'r_4041';
biomassIndex = 717;
growth = data$log2relT  #Extract the growth rates
strains = data$Row #Extract the strains names 
data$Row <- NULL
data$log2relT <- NULL #
expressionData[,1] <- NULL #Remove the growth rate data

iRFIndicies = read.csv('data/Features_Extracted_Using_iRF.csv', header = F)
SGLIndicies = read.csv('data/Features_Extracted_Using_SGL.csv', header = F) 
GeneticAlgoIndicies = read.csv('data/genetic_feature_selection_features.csv', header = FALSE,  stringsAsFactors=FALSE) 

iRFData = data[, which(colnames(data) %in% unlist(iRFIndicies[,1])) ]
SGLData = data[, which(colnames(data) %in% unlist(SGLIndicies[,1]))]
fluxData <- data[, -which(colnames(data) %in% colnames(expressionData))] #flux data is all data without gene data


#########################################################################################
#Data partition #########################################################################
#########################################################################################


#Full data partitions
testingIndices<- unlist(read.csv('testing_index.csv', header = F))
testingIndices <- testingIndices + 1 # to account for index changes between R and Python
testingTarget <- growth[testingIndices]
trainingTarget <- growth[-testingIndices]
trainingData <- data[-testingIndices,]
testingData <-  data[testingIndices,]
normParam <- preProcess(trainingData)
normalisedTestData <- predict(normParam, testingData)
normalisedTrainingData <- predict(normParam, trainingData)


#Expression data partition 
expressionTrainingData <- expressionData[-testingIndices,]
expressionTestingData <- expressionData[testingIndices,]
normParamExp <- preProcess(expressionTrainingData)
normalisedExpressionTrainingData <- predict(normParamExp, expressionTrainingData)
normalisedExpressionTestingData <- predict(normParamExp, expressionTestingData)


#Flux data partition 
fluxTrainingData <- fluxData[-testingIndices,]
fluxTestingData <- fluxData[testingIndices,]
normParamFlux <- preProcess(fluxTrainingData)
normalisedFluxTrain <- predict(normParamFlux, fluxTrainingData)
normalisedFluxTest <- predict(normParamFlux, fluxTestingData)

#Metabolic expression data partition 
metabolicExpressionTrainingData <- metabolicExpression[-testingIndices,]
metabolicExpressionTestingData <- metabolicExpression[testingIndices,]
normParamMetabolicExpression <- preProcess(metabolicExpressionTrainingData)
normalisedMetabolicExpressionTrain <- predict(normParamMetabolicExpression, metabolicExpressionTrainingData)
normalisedMetabolicExpressionTest <- predict(normParamMetabolicExpression, metabolicExpressionTestingData)

#SGL expression data partition 
SGLTrainingData <- SGLData[-testingIndices,]
SGLTestingData <- SGLData[testingIndices,]
normParamSGL <- preProcess(SGLTrainingData)
normalisedSGLTrain <- predict(normParamSGL, SGLTrainingData)
normalisedSGLTest <- predict(normParamSGL, SGLTestingData)

#iRF expression data partition 
iRFTrainingData <- iRFData[-testingIndices,]
iRFTestingData <- iRFData[testingIndices,]
normParamiRF <- preProcess(iRFTrainingData)
normalisediRFTrain <- predict(normParamiRF, iRFTrainingData)
normalisediRFTest <- predict(normParamiRF, iRFTestingData)

# Genetic algo features, need to check the normalisation is correct - come back to this
gen_features_index <- read.csv('data/genetic_feature_selection_features.csv', header = F)
normGenTrain <- vector(mode = "list", length = dim(gen_features_index)[1] - 1)
normGenTest <- vector(mode = "list", length = dim(gen_features_index)[1] - 1)

for (j in c(1:9)){
  features <- word(unlist(gen_features_index[,j]))
  print(features)
  tr = data[-testingIndices, features]
  te = data[testingIndices, features]
  normParam <- preProcess(tr)
  normGenTrain[[j]] = predict(normParam, tr)
  normGenTest[[j]] = predict(normParam, te)
}


#transfer learned 
transfer_learned_train <- read.csv('data/transfer_learned_train.csv', header = F)
transfer_learned_test <- read.csv('data/transfer_learned_test.csv', header = FALSE)
normTF <- preProcess(transfer_learned_train)
transfer_learned_train <- predict(normTF, transfer_learned_train)
transfer_learned_test <- predict(normTF, transfer_learned_test)

#####################################################################################################
##### BEMKL #########################################################################################
#####################################################################################################

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

run_bemkl_experiments <- function(){
  for (i in c(0:100)){
    #initalize the parameters of the algorithm
    parameters <- list()
    
    #set the hyperparameters of gamma prior used for sample weights
    parameters$alpha_lambda <- 1e-10
    parameters$beta_lambda <- 1e+10
    parameters$seed <- i**2 
    
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
    #set the number of iterations
    parameters$iteration <- 200
    #determine whether you want to store the lower bound values
    parameters$progress <- 0
    #set the seed for random number generator used to initalize random variables
    # parameters$seed <- 1606
    all_kernels <- c(rbf_exp_1, rbf_exp_2, rbf_exp_3,  rbf_exp_4, rbf_exp_5, rbf_exp_6, rbf_reac_1, rbf_reac_2, rbf_reac_3, rbf_reac_4, rbf_reac_5, rbf_reac_6)  
    all_test_kernels <- c(rbf_exp_1t, rbf_exp_2t, rbf_exp_3t, rbf_exp_4t, rbf_exp_5t, rbf_exp_6t, rbf_reac_1t, rbf_reac_2t, rbf_reac_3t,  rbf_reac_4t, rbf_reac_5t, rbf_reac_6t)  
    Ktrain <- array(all_kernels, c(dim(rbf_exp_1),12))
    Ktest <- array(all_test_kernels, c(dim(rbf_exp_1t), 12))
    ytrain <- trainingTarget
    print(i)
    state <- bemkl_supervised_regression_variational_train(Ktrain, ytrain, parameters)
    prediction <-  bemkl_supervised_regression_variational_test(Ktest, state)$y$mu
    print(prediction)
    write.csv(prediction, paste('predictions/BEMKL/metabolic_expression_bemkl_Predictions_', i, '.csv',sep = ''))
    RMSE <- postResample(prediction, testingTarget)["RMSE"]
    print(RMSE)
  }
  
}

run_bemkl_experiments()

########################################################################################
#SVM Gaussian ##########################################################################



#Param exploration 
svmTune <- expand.grid(C = c(10,50,100,250,500,1000, 2000) , sigma = c(0.00001, 0.0001, 0.001, 0.1, 1) )

#Training an svm 100 times 
train_100_svm <- function(train_data, train_target, test_data, test_target, data_name){
  for (i in c(0:100)){
    print(i)
    #create a list of seed, here change the seed for each resampling
    set.seed(i ** 2)
    #length is = (n_repeats*nresampling)+1
    seeds <- vector(mode = "list", length = 16)
    #(35 is the number of tuning parameters)
    for(j in 1:15) seeds[[j]]<- sample.int(n=1000, 35)
    #for the last model
    seeds[[16]]<-sample.int(1000, 1)
    ctr <- trainControl(train_target, method = "repeatedcv", number = 5, repeats =  3, verboseIter = FALSE, seeds=seeds) # we change the seed to get variation (I hope)
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    param_explore_res <- train(train_data, 
                          train_target,
                          method = 'svmRadial',
                          trControl=ctr,
                          tuneGrid = svmTune)
    stopCluster(cl)
    predictions <- predict(param_explore_res$finalModel, test_data)
    rmse <- postResample(pred = predictions, obs = test_target)
    print(cor(predictions, test_target) )
    print(MedianAE(predictions, test_target))
    print(rmse)
    write.csv(predictions, paste('predictions/SVM/', data_name, '_Predictions_', i, '.csv',sep = ''))
    }
}

# Fluxes SVM
train_100_svm(normalisedFluxTrain, trainingTarget, normalisedFluxTest, testingTarget, 'fluxes_svm')
# Metabolic Expression SVM
train_100_svm(normalisedMetabolicExpressionTrain, trainingTarget, normalisedMetabolicExpressionTest, testingTarget, 'metabolic_expression_svm')
# Expression only SVM
train_100_svm(normalisedExpressionTrainingData, trainingTarget, normalisedExpressionTestingData, testingTarget, 'expression_svm')

# Iterative RF SVM
train_100_svm(normalisediRFTrain, trainingTarget, normalisediRFTest, testingTarget, 'iRF_svm')
# SGL SVM
train_100_svm(normalisedSGLTrain, trainingTarget, normalisedSGLTest, testingTarget, 'sgl_svm')
# Concatenated data SVM 
train_100_svm(normalisedTrainingData, trainingTarget, normalisedTestData, testingTarget, 'concat_expression_fluxes_svm')

# Transfer Learned 
train_100_svm(transfer_learned_train, trainingTarget, transfer_learned_test, testingTarget, 'mm_transfer_learned_svm')



#########################################################################################
#Random Forest ##########################################################################

#Parameters explored in random forest
rfgrid <- expand.grid(.mtry=seq(2,500,20))


#Random forest training, 100 times 
train_100_rf <- function(train_data, train_target, test_data, test_target, data_name){
  for (i in c(0:100)){
    print(i)
    #create a list of seed, here change the seed for each resampling
    set.seed(i)
    #length is = (n_repeats*nresampling)+1
    seeds <- vector(mode = "list", length = 16)
    #(25 is the number of tuning parameters)
    for(j in 1:15) seeds[[j]]<- sample.int(n=1000, 25)
    #for the last model
    seeds[[16]]<-sample.int(1000, 1)
    ctr <- trainControl(index=createFolds(train_target), method = "repeatedcv", number = 5, repeats =  3, verboseIter = FALSE, seeds=seeds) # we change the seed to get variation (I hope)
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    param_explore_res <- train(train_data, 
                               train_target,
                               method = "rf",
                               trControl=ctr,
                               tuneGrid = rfgrid)
    stopCluster(cl)
    predictions <- predict(param_explore_res$finalModel, test_data)
    rmse <- postResample(pred = predictions, obs = test_target)
    print(rmse)
    write.csv(predictions, paste('predictions/RF/', data_name, '_Predictions_', i, '.csv',sep = ''))
  }
  stopCluster(cl)
}

# Fluxes RF
train_100_rf(normalisedFluxTrain, trainingTarget, normalisedFluxTest, testingTarget, 'fluxes_rf')
# Metabolic Expression RF
train_100_rf(normalisedMetabolicExpressionTrain, trainingTarget, normalisedMetabolicExpressionTest, testingTarget, 'metabolic_expression_rf')
# Expression only RF
train_100_rf(normalisedExpressionTrainingData, trainingTarget, normalisedExpressionTestingData, testingTarget, 'expression_rf')
# Iterative iRF RF
train_100_rf(normalisediRFTrain, trainingTarget, normalisediRFTest, testingTarget, 'iRF_rf')
# SGL RF
train_100_rf(normalisedSGLTrain, trainingTarget, normalisedSGLTest, testingTarget, 'sgl_rf')
# Concatenated data  RF
train_100_rf(normalisedTrainingData, trainingTarget, normalisedTestData, testingTarget, 'concat_expression_fluxes_rf')


# Late integration RF models 
late_integraion_rf <- function(link1, link2, outlink){
  for (i in c(0:100)){
    first = read.csv(paste('predictions/RF/', link1, '_Predictions_', i, '.csv',sep = ''))
    second = read.csv(paste('predictions/RF/', link2, '_Predictions_', i, '.csv',sep = ''))
    out = (first[,2] + second[,2]) / 2
    rmse <- postResample(pred = out, obs = testingTarget)
    print(rmse)
    write.csv(out, paste('predictions/RF/', outlink, '_Predictions_', i, '.csv',sep = ''))
  }
}

# late integration for RF 
late_integraion_rf('expression_rf', 'fluxes_rf', 'bagged_expression_fluxes_rf')
late_integraion_rf('metabolic_expression_rf', 'fluxes_rf', 'bagged_metabolic_expression_fluxes_rf')


#######################################################################################
# Genetic features ####################################################################

#Random forest for genetic (NSGA-II) features
train_100_rf_genetic <- function(train_data_list, train_target, test_data_list, test_target, data_name){
  for (i in c(9)){
   for (gen in c(11)){
      s <- i*9 + gen - 1 
      train_data <- train_data_list[[gen]]
      test_data <- test_data_list[[gen]]
    #create a list of seed, here change the seed for each resampling
    set.seed(i)
    #length is = (n_repeats*nresampling)+1
    seeds <- vector(mode = "list", length = 16)
    #(25 is the number of tuning parameters)
    for(j in 1:15) seeds[[j]]<- sample.int(n=1000, 25)
    #for the last model
    seeds[[16]]<-sample.int(1000, 1)
    ctr <- trainControl(index=createFolds(train_target), method = "repeatedcv", number = 5, repeats =  3, verboseIter = FALSE, seeds=seeds) # we change the seed to get variation (I hope)
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    param_explore_res <- train(train_data, 
                               train_target,
                               method = "rf",
                               trControl=ctr,
                               tuneGrid = rfgrid)
    stopCluster(cl)
    predictions <- predict(param_explore_res$finalModel, test_data)
    rmse <- postResample(pred = predictions, obs = test_target)
    print(rmse)
    write.csv(predictions, paste('predictions/RF_GEN/', data_name, '_Predictions_', s, '.csv',sep = ''))
    }
  }
}


# SVM for the genetic features
train_100_svm_genetic <- function(train_data_list, train_target, test_data_list, test_target, data_name){
  for (i in c(1:9)){
    for (gen in c(1:11)){
      s <- i*9 + gen - 1 
    train_data <- train_data_list[[gen]]
    test_data <- test_data_list[[gen]]
    #create a list of seed, here change the seed for each resampling
    set.seed(i ** 2)
    #length is = (n_repeats*nresampling)+1
    seeds <- vector(mode = "list", length = 16)
    #(35 is the number of tuning parameters)
    for(j in 1:15) seeds[[j]]<- sample.int(n=1000, 35)
    #for the last model
    seeds[[16]]<-sample.int(1000, 1)
    ctr <- trainControl(index=createFolds(train_target), method = "repeatedcv", number = 5, repeats =  3, verboseIter = FALSE, seeds=seeds) # we change the seed to get variation (I hope)
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    param_explore_res <- train(train_data, 
                               train_target,
                               method = 'svmRadial',
                               trControl=ctr,
                               tuneGrid = svmTune)
    stopCluster(cl)
    predictions <- predict(param_explore_res$finalModel, test_data)
    rmse <- postResample(pred = predictions, obs = test_target)
    print(rmse)
    write.csv(predictions, paste('predictions/SVM_GEN/', data_name, '_Predictions_', s, '.csv',sep = ''))
    }
  }
}


train_100_rf_genetic(normGenTrain, trainingTarget, normGenTest, testingTarget, 'genetic_rf')
train_100_svm_genetic(normGenTrain, trainingTarget, normGenTest, testingTarget, 'genetic_svm')


