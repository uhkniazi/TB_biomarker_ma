# Name: integrate_cgraphclust_results.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 10/08/2015
# Desc: combined tb microarray data from long testing with CGraphClust results

source('tb_biomarker_ma_header.R')


## data loading
# load the data, clean and create factors
dfExp = f_LoadObject(file.choose())

# load the sample annotation
dfSamples = read.csv(file.choose(), header=T)

# sort both the samples and expression data in same order
rownames(dfSamples) = as.character(dfSamples$Sample_ID)
dfSamples = dfSamples[colnames(dfExp),]

# load annotation data
dfAnnotation = f_LoadObject(file.choose())
rownames(dfAnnotation) = as.character(dfAnnotation$ID)

# load the results from CGraphClust
f = file.choose()
dfCluster = read.csv(f, header=T, row.names=1)

# create factors
fGroups = factor(dfSamples$Illness1)
# create a second factor with only 2 levels
# keep ptb at 1 for downstream predictions
fGroups.2 = as.character(dfSamples$Illness)
fGroups.2 = factor(fGroups.2, levels = c('OD', 'ATB'))
dfSamples$fGroups = fGroups
dfSamples$fGroups.2 = fGroups.2

# get the genes of interest
csCluster = 'R-HSA-1280218'

csGenes = as.character(dfCluster[dfCluster$cluster == csCluster,'ENTREZID'])

# get the corresponding probe ids
f = dfAnnotation$Enterez_ID %in% csGenes
csID = as.character(dfAnnotation$ID[f])

mDat = as.matrix(dfExp[csID,])
# phenotypic data
fGroups.2 = dfSamples$fGroups.2

# remove any NA
mDat = na.omit(mDat)
fRepeated.probes = factor(as.character(dfAnnotation[rownames(mDat), 'Enterez_ID']))
# take mean for groups that have a repeated probe
mDat.m = apply(mDat, 2, FUN = function(x) tapply(x, fRepeated.probes, mean))
dfExport = as.data.frame(t(mDat.m))
dfExport$fGroups = dfSamples$fGroups
dfExport$fGroups.2 = dfSamples$fGroups.2 

write.csv(dfExport, file=paste('Temp/', csCluster, '.csv', sep = ''))

# split the data into test and training sets
set.seed(123)
dfImport = read.csv(file.choose(), header=T, row.names=1)
test = sample(1:nrow(dfImport), size = nrow(dfImport) * 0.3, replace = F)
dfData = dfImport[,1:(ncol(dfImport)-2)]

############# variable selection steps
# note: change method to exhaustive or forward
oVar.sub = CVariableSelection.ReduceModel(dfData[-test,], dfImport$fGroups.2[-test], boot.num = 10, cvMethod = 'forward')
plot.var.selection(oVar.sub)
par(mfrow=c(2,2))
# try models of various sizes with CV
for (i in 2:7){
  cvTopGenes.min = CVariableSelection.ReduceModel.getMinModel(oVar.sub, 3, cvMethod = 'forward')
  dfData.train = dfData[-test, cvTopGenes.min]
    
  dfData.test = dfData[test, cvTopGenes.min]
  
  oCV = CCrossValidation.LDA(test.dat = dfData.test, train.dat = dfData.train, test.groups = dfImport$fGroups.2[test],
                             train.groups = dfImport$fGroups.2[-test], level.predict = 'ATB', boot.num = 100)
  
  plot.cv.performance(oCV)
}
par(p.old)

for (i in 2:7) {print(paste('size', i)); print(CVariableSelection.ReduceModel.getMinModel(oVar.sub, i, cvMethod='forward'))}

# # plot using a tree based classifier
# par(mfrow=c(2,2))
# # try models of various sizes with CV
# for (i in 1:9){
#   cvTopGenes.min = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
#   dfData.train = as.data.frame(mDat.sub.train)
#   dfData.train = dfData.train[,colnames(dfData.train) %in% c(cvTopGenes.sub, cvTopGenes.min)]
#   
#   dfData.test = as.data.frame(mDat.sub.test)
#   dfData.test = dfData.test[,colnames(dfData.test) %in% c(cvTopGenes.sub, cvTopGenes.min)]
#   
#   oCV = CCrossValidation.Tree(test.dat = dfData.test, train.dat = dfData.train, test.groups = dfSamples.test$fGroups.2,
#                               train.groups = dfSamples.train$fGroups.2, level.predict = 'ATB', boot.num = 50)
#   
#   plot.cv.performance(oCV)
# }
# par(p.old)
# 
# # testing 2 gene combinations
# gene.cor = cvTopGenes
# par(mfrow=c(2,2))
# for (i in 1:length(gene.cor)){
#   dfData.train = as.data.frame(mDat.sub.train)
#   dfData.train = data.frame(dfData.train[,colnames(dfData.train) %in% c(cvTopGenes.sub, gene.cor[i])])
#   dfData.test = as.data.frame(mDat.sub.test)
#   dfData.test = data.frame(dfData.test[,colnames(dfData.test) %in% c(cvTopGenes.sub, gene.cor[i])])
#   oCV = CCrossValidation.LDA(test.dat = dfData.test, train.dat = dfData.train, test.groups = dfSamples.test$fGroups.2,
#                              train.groups = dfSamples.train$fGroups.2, level.predict = 'ATB', boot.num = 100)
#   
#   plot.cv.performance(oCV)  
# }
# 
# as.data.frame((paste(cvTopGenes.sub, gene.cor)))












