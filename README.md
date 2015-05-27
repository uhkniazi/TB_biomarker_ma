# TB_biomarker_ma
Selection of public microarray datasets used to select set of genes for class prediction of TB vs Other diseases.

# tb_biomarker_ma_header.R
Header file with libraries, functions and global variables to set. Keep in same directory as other scripts
before running a script.

# berry_et_al_nature.R
Analysis of 'GSE19491' dataset from GEO database (microarray data) of TB and other diseases, whole blood and 
sorted cells gene expression data. The script uses lumi to log2 transform the data, while also adding an offset
to negative values to bring them up above 0 before taking the log. lumi is used to normalise the data and pca
used to check for outliers and removing those samples. The count matrix is scaled (i.e. standardized on sd and
centered to 0 mean) for each sample before further analysis. The genes with a coefficient of variation upto the 40% 
quantile are selected for further analysis. t test and wilcoxon test with BH FDR of 10% is used as a cutoff for 
significantly different genes between PTB and Other conditions. Using this smaller set of genes perform variable reduction,
using random forest (NOTE: this can take several hours to finish). The smaller set of variables are checked again by 
running a simulation of iBoot steps, where the class data is sampled so that proportions are 50:50 for each class, the 
model is fit multiple times in order to assess the error rate and calculate its standard deviation. The results are used
to further reduce the variable size. Fit an lda model with nested k fold cross validation (50:50 for each class) at every
cycle, and fit a ROC curve for the prediction of PTB and cross validation error rate. 

NOTE: further variable reduction can be done maually by choosing variables that are least correlated with each other
by analysing the correlation matrix.

# kaforou_et_al_plosm.R
Analysis of 'GSE37250'. The analysis protocol is similar to the previous berry script.

# long_all_data_test.R
Similar to next script, use second one its latest.

# long_all_data_test_02.R
Analysis of the combined data, provided by long. The annotation for the probes is loaded from Geo Database, while the 
data from long and sample annotation are loaded separately. We remove sarcoidosis and non-active sarcoids. The factors are relevelled
so that ATB is 1 and OD is 0. PCA on the samples is performed and some outliers are removed. Data is z-scale standardized across
the samples before further analysis, and infinite and na values are removed. The genes are subset into groups based on the 
coefficient of variation which shows how precise a value is. Genes are further subset, based on if they are differentially expressed
between the 2 conditions on a t.test and a wilcoxon test. Data is split into training and test set (20%). Random Forest is performed
on the training set to select a set of genes that may be important, based on the variable importance distribution. Using that smaller
set of genes, a nested random forest with equal proportions of both classes, is performed to estimate the importance score, its mean
and standard deviation and coef of variation. Genes with lowest coef of variation are selected. Further genes are removed based
on the colinearity. Another nested variable subset selection step with data split into training and test sets is peformed and 
distribution of test error rates for each number of variables is calculated. Smallest number of variables model with lowest or
more precise error rate is selected. Using those variables a nested 10 fold cross validation is performed on the training data
and a test/training error rate is also calculated for the test data. Summary is plotted as a ROC curve and boxplots.

# long_all_data_clusters.R
Making clusters of related genes and pathways using the expression data saved as a list from the previous script. Choose various 
types of classes (groups) to make clusters. Not a fixed script, needs adjustment according to the question.

# kaforou_et_al_nejm.R
Very similar to previous analysis, however we use the dataset GSE39941, from the transcriptomics study on TB in children with and
without HIV. The PCA plot suggests a few outliers which are removed before normalization using lumi. The data is split into 2 groups
based on all other conditions (OD) and active tb (ATB). The data is split into training and test sets randomly with 20% of data 
as test set. Using the markers selected earlier (in previous scripts) - they are tested using 10 fold cross validation and training
vs test error rates and ROC curves are plotted. The expression values for these genes are also plotted. A pathway analysis is also
performed i.e. clustering based on shared pathways and expression values of the selected 2500 or so markers. The pathway clustering
is very similar to long_all_data_clusters.R and various plots are produced to for significant modules.
