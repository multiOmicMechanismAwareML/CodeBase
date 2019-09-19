# Machine Learning Predicts Yeast Growth

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

The deep learning techniques were applied using Python and are found in the Code folder. All other techniques applied were done so using R and can be found in the machineLearning.R file also inside the Code folder.

