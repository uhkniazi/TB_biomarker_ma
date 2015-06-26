# Name: PCR_data_CV.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 26/06/15
# Desc: pcr data 

source('../CCrossValidation//CCrossValidation.R')
library(MASS)
library(ROCR)


dfDat = read.csv(file.choose())

fGroups = dfDat$Diagnosis
gene1 = dfDat[,4]
gene2 = dfDat[,6]

dfDat = data.frame(fGroups, gene1, gene2)
dfDat = na.omit(dfDat)

i = grep('Sarcoid', dfDat$fGroups)
dfDat = dfDat[-i,]

fGroups = rep(NA, length=length(dfDat$fGroups))
i = grep('Active TB$', dfDat$fGroups)
fGroups[i] = 'ATB'
fGroups[-i] = 'OD'

dfDat$fGroups = factor(fGroups)
i = grep('N/A', dfDat$gene2)
dfDat = dfDat[-i,]
dfDat$gene2 = as.numeric(as.character(dfDat$gene2))

i = which(dfDat$gene2 == 0)
dfDat = dfDat[-i,]

test = sample(1:nrow(dfDat), size = nrow(dfDat) * 0.2)
fGroups = dfDat$fGroups
mDat = dfDat[,c(2,3)]

oCV = CCrossValidation.LDA(data.frame(gene=mDat[test,]), data.frame(gene=mDat[-test,]), fGroups[test], fGroups[-test], level.predict = 'ATB', boot.num = 100)
svg('Temp/cv2.svg')
plot.cv.performance(oCV)
dev.off(dev.cur())

