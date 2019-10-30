
# A Multi-omic and Mechanism-aware Machine Learning Pipeline Predicts Growth Rate in Saccharomyces cerevisiae

## Objective 

To leverage strain-specific metabolic modelling coupled with gene expression data to produce predictive models for Yeast growth.

**Step one:**

msb145172 dataset (Duibhir et al., 2014) is used with Metrade to create strain-specific metabolic models through altering the reaction upper and lower bound constraints inside a metabolic model. 

**Step two:**

Min-norm flux balance analysis is run to produce a simulated cell network from which the reaction flux rates are extracted

**Step three:**

Feature selection techniques are applied, we explore three state-of-the-art techniques :

* NSGA-II - A genetic algorithm that allows for multiple objectives
* Sparse group lasso - A regularisation technique which promotes group sparsity
* Iterative random forest - Selected for its abilit to account for non-linear interations 

**Step four:**

Machine learning and deep learning techniques are applied to produce predictive models for both the reduced feature data, single-omics and concatenated using early, intermediate and late data integration techniques. We explore: 

* Random forest - The standard setup for the datasets and also a late integration bagged vote model 
* Deep Learning - both a two layer feed forward network and a multi-modal network trained on separate omics 
* Support vector regression - We also explore a the Bayesian Efficient Multiple Kernel Learning algorithm as a multimodal option when combining two single omic datasets. 


## Files

### Metabolic Modelling 

The metabobolic modelling code is run in Matlab (version R2015b) - Code/metabolic_modelling/RUN.m 

This will automatically start the FBA process using the msbData and the yeastmmm.mat metabolic model 

The captured flux levels can then be extracted and used to create a complete dataset (see our experiment example in Code/learning/Data)


### Deep Learning Approach

In order to run the deep learning models two zip files need to be extracted: 

* Code/learning/Data/CompleteDataset.csv

* Code/learning/Data/Reduced_Dataset.csv

#### Requirements:

* Python 3.5

* Tensorflow 2.0 

* Keras 2.0

* Numpy 1.5

* Pandas 0.25

* Seaborn 0.9

* Matplotlib 3.1.1


To run the deep learning models - once these datasets have been extracted - simply run the python file Code/learning/DeepLearningFullSet.py to see the set of results. 


### Tradtional Machine Learning Approaches

The R file Code/learning/machineLearning.R contains the feature selection and traditional machine learning approaches (random forest and support vector regression models) which can be run by first setting your workstation to point its data collection to Code/learning/Data and then running the complete code. 

This code requires that the libraries listed at the top of the code first be installed to ensure it runs to completion.






