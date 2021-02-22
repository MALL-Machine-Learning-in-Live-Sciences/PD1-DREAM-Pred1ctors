#!/usr/bin/env Rscript
setwd('~/git/PD1-DREAM-Pred1ctors/pfs_v3/')

# Progression Free Survival (PFS)
# =====================================
# load packages
require(data.table)
require(sva)
require(survival)
require(impute)
source('utils.r')

# get command lines args: args[1] rnas-seq gene level count data; args[2] clinical data; args[3] output file
# args = commandArgs(trailingOnly=TRUE)
counts.test = fread('~/git/PD1-DREAM-Pred1ctors/CM_026_formatted_synthetic_data_subset/GRCh37ERCC_refseq105_genes_count.csv',data.table = F); rownames(counts.test) <- counts.test[,1]; counts.test <- counts.test[,-1]; print("done reading in counts")
counts.test = as.data.frame(t(counts.test))
clinical.test = read.csv('~/git/PD1-DREAM-Pred1ctors/CM_026_formatted_synthetic_data_subset/clinical_data.csv', header = T, row.names = 1); print('Done reading clinical data!')

# load ours objects
trainObj = readRDS('models/pfs_biospear_alternative.rds')

pcaModel = trainObj$pca
trainData = trainObj$train
counts.train = trainObj$genes
genes = names(counts.train)
model_PCA = trainObj$model
model_TMB = readRDS('models/pfs_model3.rds')$model
print('Done load training objects!')

# Match genes in test count data
counts.test = counts.test[, match(genes, names(counts.test))]
print('Done match genes!')

# Standarize and scale
counts.test = standarize(counts.test, log2 = T, scale = T)
print('Done standarization!')

# Impute NA�s
impNA = impute.knn(t(counts.test))
counts.test = as.data.frame(t(impNA$data))
print('Done NA imputing!')

# ComBat
meta = data.frame(patID = c(rownames(counts.train), rownames(counts.test)),
                  subset = c(rep('train', nrow(counts.train)), rep('test', nrow(counts.test))))
rownames(meta) = meta$patID
counts = rbind.data.frame(counts.train, counts.test)
m = model.matrix(~1, data = meta)
combat = ComBat(dat = t(counts), batch = meta$subset, mod = m, ref.batch = 'train')
combat = as.data.frame(t(combat))
counts.test = combat[match(rownames(counts.test), rownames(combat)), ]
print('Done ComBat!')

# Predict PCA
test.pca = predict(pcaModel, newdata = counts.test)
test.pca = as.data.frame(test.pca)
print('Done predict PCA�s in test data!')

# Format clinical data
clinical.test$SEX[which(clinical.test$SEX == 'M')] = 2
clinical.test$SEX[which(clinical.test$SEX == 'F')] = 1

clinical.test$CRFHIST[which(clinical.test$CRFHIST == 'NON-SQUAMOUS')] = 1
clinical.test$CRFHIST[which(clinical.test$CRFHIST == 'SQUAMOUS')] = 4

# Add Clinical data
test.pca$Cohort = 4.5
test.pca$Age = round(clinical.test$AAGE, 0)
test.pca$Gender = as.numeric(clinical.test$SEX)
test.pca$tumor.type = as.numeric(clinical.test$CRFHIST)
test.pca$Drug = 0
print('Done add clinical data!')

# Predict module PCA�s
method = 'alassoU'
pred1 = predBiospear(res = model_PCA, method = method, newdata = test.pca, randomPFS = T)[,1]
print('Done prediction in module 1!')

# Predict module 3
clinical.test$TMB = as.numeric(clinical.test$TMB)
TMB_c = ifelse(clinical.test$TMB <= 100, 'low', 
               ifelse(clinical.test$TMB < 242 & clinical.test$TMB > 100, 'medium', 'high'))
TMB_h = ifelse(TMB_c == 'high', 'yes', 'no')

data3 = data.frame(
  Age = round(as.numeric(clinical.test$AAGE), 0),
  TMB = TMB_h,
  Drug = 'anti-PD1',
  row.names = rownames(clinical.test)
)
pred3 = predict(model_TMB, data3)
pred3[which(is.na(pred3))] = 0
print('Done prediction in module 3!')

# Final predictions
pred = pred1 + pred3
print('Predictions summatory!')

res = cbind.data.frame(patientID = rownames(clinical.test),
                       prediction = pred)

print(res)
# write.csv(res, file = args[3], quote = F, row.names = F); print("done writing out signature")
