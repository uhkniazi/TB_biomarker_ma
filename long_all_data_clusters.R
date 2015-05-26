# Name: long_all_data_clusters.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 22/05/15
# Desc: pathway and cluster analysis of all combined tb ma data


### libraries
library(annotate)
library(org.Hs.eg.db)
source('~/Dropbox/Home/Data/R/My_Libraries/NGS_functions.R')
library(igraph)
source('../CGraphClust/CGraphClust.R')
library(reactome.db)

### data loading
# data from previous script long_all_data_test.R
lDat = f_LoadObject(file.choose())
lDat$desc

# annotation
dfAnnotation = f_LoadObject(file.choose())

# pheno data
dfSamples = lDat$sample
# expression data
mExp = lDat$expression


### processing
# assign names to the genes from illumina ids to enterez ids
n = colnames(mExp)
i = match(n, dfAnnotation$ID)
# replacement names i.e. enterez ids
n2 = as.character(dfAnnotation$Enterez_ID[i])
colnames(mExp) = n2

# remove any duplicated enterez ids
f = duplicated(n2)
mExp = mExp[,!f]

### create clusters based on various criteria
# difference between hiv positive and negative atb
f = dfSamples$Illness
i = which(f == 'ATB')
dfSamples.atb = dfSamples[i,]
mExp.atb = mExp[i,]

## get the reactome ids
dfGraph = select(reactome.db, keys = colnames(mExp.atb), columns = 'REACTOMEID', keytype = 'ENTREZID')
dfGraph = na.omit(dfGraph)

# remove genes from expression matrix that do not have the 
# matching enterez id and reactome id
n = unique(dfGraph$ENTREZID)
mExp.atb = mExp.atb[,n]

# create correlation matrix
mCor = cor(mExp.atb)

# choose cutoff based on histogram
hist(sample(mCor, 10000))

oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.6)
ig = getFinalGraph(oGr)

# plot results
plot(ig, vertex.size=1, vertex.label=NA, layout=layout.fruchterman.reingold)
plot(getCommunity(oGr), ig, vertex.size=1, vertex.label=NA, layout=layout.fruchterman.reingold)
plot.mean.expressions(oGr, t(mExp.atb), fGroups = dfSamples.atb$HIV_status)

### other criteria
f = dfSamples$HIV_status
i = which(f == 'HIV_NEG')
dfSamples.atb = dfSamples[i,]
mExp.atb = mExp[i,]

f = dfSamples.atb$Illness1
i = which(f %in% c('ATB', 'LTBI', 'Ctrl'))
dfSamples.atb = dfSamples.atb[i,]
mExp.atb = mExp.atb[i,]

## get the reactome ids
dfGraph = select(reactome.db, keys = colnames(mExp.atb), columns = 'REACTOMEID', keytype = 'ENTREZID')
dfGraph = na.omit(dfGraph)

# remove genes from expression matrix that do not have the 
# matching enterez id and reactome id
n = unique(dfGraph$ENTREZID)
mExp.atb = mExp.atb[,n]

# create correlation matrix
mCor = cor(mExp.atb)

# choose cutoff based on histogram
hist(sample(mCor, 10000))

oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.5)
ig = getFinalGraph(oGr)

# plot results
plot(ig, vertex.size=1, vertex.label=NA, layout=layout.fruchterman.reingold)
plot(getCommunity(oGr), ig, vertex.size=1, vertex.label=NA, layout=layout.fruchterman.reingold)
plot.mean.expressions(oGr, t(mExp.atb), fGroups = as.character(dfSamples.atb$Illness1))

# get the dataframe cluster labels
df = getClusterMapping(oGr)
df = df[order(df$type.2),]
write.csv(df, file='Temp/df.csv')

temp = scan(what = character())
i = which(dfAnnotation$Enterez_ID %in% temp)
dfAnnotation[i,]


