# tb_biomarker_ma_header.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 23/04/2015
# Desc: various functions used in scripts

## libraries and scripts
library(GEOquery)
library(Biobase)
library(annotate)
library(org.Hs.eg.db)
library(lumi)
library(MASS)
library(randomForest)
library(tree)
library(ROCR)
source('~/Dropbox/Home/Data/R/My_Libraries/NGS_functions.R')
library(leaps)
library(caret)
source('../CCrossValidation/CCrossValidation.R')
#library(NMF)

# global variables
p.old = par()

## functions
# Name: f_Plot3DPCA
# Args: mComp = n X 3 matrix with first 3 components as the 3 column vectors
#       color = colours for the points
#       ... additional arguments to the plot function
# Rets: none
# Desc: takes the first 3 components and plots the data in a 3d plot
f_Plot3DPCA = function(mComp, color, ...) {
  x = mComp[,1]
  y = mComp[,2]
  z = mComp[,3]
  if (!require(scatterplot3d)) stop('scatterplot3d library required')
  scatterplot3d(x, y, z, color, ...)
}