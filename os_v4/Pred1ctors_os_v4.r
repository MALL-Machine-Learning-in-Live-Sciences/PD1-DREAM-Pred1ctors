#!/usr/bin/env Rscript
setwd('~/git/PD1-DREAM-Pred1ctors/os_v4/')

# Overall Survival (OS)
# =====================================
# load packages
require(data.table)
require(sva)
require(survival)
require(impute)

# Util function
standarize = function(df, log2 = T, scale = T){

  if (log2 == T){
    print('Applying log2 transformation...')
    df = apply(df, 2, function(x) log2(x + 1))
  }
  
  if (scale == T){
    print('Scaling dataset...')
    df = scale(df)
  }
  
  res = df
  return(res)
}

# get command lines args: args[1] rnas-seq gene level count data; args[2] clinical data; args[3] output file
args = commandArgs(trailingOnly=TRUE)
counts.test = fread('~/git/PD1-DREAM-Pred1ctors/CM_026_formatted_synthetic_data_subset/GRCh37ERCC_refseq105_genes_count.csv',data.table = F); rownames(counts.test) <- counts.test[,1]; counts.test <- counts.test[,-1]; print("done reading in counts")
counts.test = as.data.frame(t(counts.test))
clinical.test = read.csv('~/git/PD1-DREAM-Pred1ctors/CM_026_formatted_synthetic_data_subset/clinical_data.csv', header = T, row.names = 1); print('Done reading clinical data!')

# load ours objects
trainObj = readRDS('models/os_cox_v3.rds')

trainData = trainObj$train
cox = trainObj$model
print('Done load training objects!')

# Select genes
genes = c("TNFAIP3", "ATG7", "CD1A", "LTK", "PTGS2", "VEGFC")
counts.test = counts.test[, match(genes, names(counts.test))]
print('Done extract genes!')

# Standarize and scale
counts.test = standarize(counts.test, log2 = T, scale = T)
print('Done standarization!')

# Impute NAs
impNA = impute.knn(t(counts.test))
counts.test = as.data.frame(t(impNA$data))
print('Done NA imputing!')

# Add Clinical data
counts.test$Drug = 'anti-PD1'
counts.test$Age = round(clinical.test$AAGE, 0)
print('Done add clinical data!')

pred = predict(cox, counts.test)
print('Done prediction!')

res = cbind.data.frame(patientID = rownames(counts.test),
                       prediction = - pred)

print(res)
# write.csv(res, file = args[3], quote = F, row.names = F); print("done writing out signature")