# Name: PCR_data_CV_nov10_2015.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 10/11/2015
# Desc: pcr data cross validation check

source('../CCrossValidation//CCrossValidation.R')


dfDat = read.csv('Data_external/Long/Nov_10/Data for Luis.csv', header=T)

# grouping factor
fGroups = dfDat$Diagoutcome
gene1 = dfDat$
gene2 = dfDat$gene2

dfData.load = dfDat

# create dataframe for model fitting
dfDat = data.frame(fGroups, gene1, gene2)
dim(dfDat)
dfDat = na.omit(dfDat)
dim(dfDat)

# remove sarcoids
i = grep('Sarcoid', dfDat$fGroups)
dfDat = dfDat[-i,]

# relevel and set new factor
fGroups = rep(NA, length=length(dfDat$fGroups))
i = grep('Active TB$', dfDat$fGroups)
fGroups[i] = 'ATB'
fGroups[-i] = 'OD'

dfDat$fGroups = factor(fGroups, levels = c('OD', 'ATB'))

# create a test set vector
set.seed(123)
test = sample(1:nrow(dfDat), size = nrow(dfDat) * 0.2)
fG = dfDat$fGroups
dat = dfDat[,!(colnames(dfDat) %in% 'fGroups')]

cv = CCrossValidation.LDA(dat[test,], train.dat = dat[-test,], fG[test], fG[-test], level.predict = 'ATB', boot.num = 1000)
plot.cv.performance(cv)

# 95% CI for auc
x = cv@oAuc.cv
x = as.numeric(x@y.values)
print(signif(quantile(x, probs = c(0.025, 0.975)), 2))

## repeat the analysis with TST included
dfDat = dfData.load
fGroups = dfDat$Diagoutcome
gene1 = dfDat$gene1
gene2 = dfDat$gene2
tst = dfDat$TSTpos
dfDat = data.frame(fGroups, gene1=log(gene1), gene2=log(gene2), tst)
dim(dfDat)
dfDat = na.omit(dfDat)
dim(dfDat)
i = grep('Sarcoid', dfDat$fGroups)
dfDat = dfDat[-i,]
i = grep('unknown', dfDat$tst)
length(i)
dfDat = dfDat[-i,]
dfDat$tst = factor(as.character(dfDat$tst))
dim(dfDat)
fGroups = rep(NA, length=length(dfDat$fGroups))
i = grep('Active TB$', dfDat$fGroups)
fGroups[i] = 'ATB'
fGroups[-i] = 'OD'
dfDat$fGroups = factor(fGroups, levels = c('OD', 'ATB'))
dim(dfDat)
summary(dfDat)
set.seed(123)
test = sample(1:nrow(dfDat), size = nrow(dfDat) * 0.2)
summary(dfDat[test,])
xtabs(~ tst+fGroups, data=dfDat[test,])
summary(dfDat[-test,])

fG = dfDat$fGroups
dat = dfDat[,!(colnames(dfDat) %in% 'fGroups')]
head(dat)
cv = CCrossValidation.LDA(dat[test,], train.dat = dat[-test,], fG[test], fG[-test], level.predict = 'ATB', boot.num = 1000)
plot.cv.performance(cv)

# 95% CI for auc
x = cv@oAuc.cv
x = as.numeric(x@y.values)
print(signif(quantile(x, probs = c(0.025, 0.975)), 2))

# repeat CV with either gene1 or gene2
# using tst + gene1
fG = dfDat$fGroups
dat = dfDat[,!(colnames(dfDat) %in% c('fGroups', 'gene2'))]
head(dat)
cv = CCrossValidation.LDA(dat[test,], train.dat = dat[-test,], fG[test], fG[-test], level.predict = 'ATB', boot.num = 1000)
plot.cv.performance(cv)

# 95% CI for auc
x = cv@oAuc.cv
x = as.numeric(x@y.values)
print(signif(quantile(x, probs = c(0.025, 0.975)), 2))

# using tst + gene2
fG = dfDat$fGroups
dat = dfDat[,!(colnames(dfDat) %in% c('fGroups', 'gene1'))]
head(dat)
cv = CCrossValidation.LDA(dat[test,], train.dat = dat[-test,], fG[test], fG[-test], level.predict = 'ATB', boot.num = 1000)
plot.cv.performance(cv)

# 95% CI for auc
x = cv@oAuc.cv
x = as.numeric(x@y.values)
print(signif(quantile(x, probs = c(0.025, 0.975)), 2))


# using tst alone
fG = dfDat$fGroups
dat = data.frame(tst=dfDat[,!(colnames(dfDat) %in% c('fGroups', 'gene2', 'gene1'))])
head(dat)
cv = CCrossValidation.LDA(data.frame(g=dat[test,]), train.dat = data.frame(g=dat[-test,]),
                          fG[test], fG[-test], level.predict = 'ATB', boot.num = 100)
plot.cv.performance(cv)

# 95% CI for auc
x = cv@oAuc.cv
x = as.numeric(x@y.values)
print(signif(quantile(x, probs = c(0.025, 0.975)), 2))

