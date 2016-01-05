Script for variable selection using the combined data sets of host blood microarray data for predicting TB and OD groups.

# Generating test data
The data sets from various studies were combined using COMBAT as explained in the main paper. The factors are relevelled
so that ATB is 1 and OD is 0. PCA on the samples is performed and some outliers are removed. Infinite and NA values are removed. Genes are further selected, based on if they are differentially expressed between the 2 conditions on a t.test and a wilcoxon test with BH FDR of .1% is used as a cutoff for significantly different genes between ATB and OD. Using this smaller set of genes we set aside 30% of data randomly as test set. Random Forest is classifier is run on the training set to select top 100 genes based on the variable importance (NOTE: this can take several hours to a day to finish, depending on computational resources). This set of 100 variables and the expression values in each sample are provided with the test script.

# Running the script
## Test Data
lData.publish.rds, provided with the repository and is downloaded when the script is run.
## Variable selection library
CCrossValidation.R is downloaded when the script is run from the repository https://github.com/uhkniazi/CCrossValidation  
This is an S4 class R object and details of this class and its functions are provided at the main repository page. 
## Steps for selection
1- Nested random forest to select top 30 variables.  
2- Nested exhaustive subset selection to select combinations of variables that produce a good fit.  
3- Nested 10 fold cross validation using LDA where the model is fit on the training set and prediction is performed on the test set for the 2 to n variable models.
