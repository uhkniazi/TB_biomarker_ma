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

NOTE: test variables predicted by one script in the data from the other data set, to see how it performs on the 
new data set.


