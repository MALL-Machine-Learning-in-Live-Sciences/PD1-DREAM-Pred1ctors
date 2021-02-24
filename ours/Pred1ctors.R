#!/usr/bin/env Rscript
# simple inflammation signature from Peter Szabo at BMS: reference below
# F. Stephen Hodi, Jedd D. Wolchok, Dirk Schadendorf, James Larkin, Max Qian, Abdel Saci, Tina C. Young, Sujaya Srinivasan, Han Chang, 
# Megan Wind-Rotolo, Jasmine I. Rizzo, Donald G. Jackson, Paolo A. Ascierto. Genomic analyses and immunotherapy in advanced melanoma. 
# In: Proceedings of the American Association for Cancer Research Annual Meeting 2019; 
# 2019 Mar 29-Apr 3; Atlanta, GA. Philadelphia (PA): AACR; Cancer Res 2019;79(13 Suppl):Abstract nr CT037.
# using NOISeq library to get tmm since it is simpler
# to run from command line: Rscript --quiet --vanilla szabo_inflammation_signature.R "input_counts_rna_seq.txt" "output.csv"

# Utils functions
standarize = function(df, log2 = T, scale = T){

  # Details:
  #   Standarize expression matrix (from counts to scale)
  # Arguments:
  #   df    : expression data.frame with patients in rows and genes in columns
  #   log2  : logical value. If TRUE, log2 normalization will be calculated. If FALSE, only scale will be carried out
  #   scale : logical value. If TRUE, scale function will be applied
  # Value:
  #   Data.frame with standarize genes

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

library(data.table)
library(sva)
library(PCAtools)
library(biospear)
library(impute)
source('predRes_modified.r')

# get command line args: args[1] input rna-seq gene level count data, args[2] output file: two column(patient id, inflammation score) .csv file
args        	 <- commandArgs(trailingOnly=TRUE)
counts.test      <- fread(args[1],data.table = F); rownames(counts.test) <- counts.test[,1]; counts.test <- counts.test[,-1]; print("done reading in counts")


# OJO!!! Delete these lines to run with docker!
# ======================================================================================
# Load validation data from local
# counts.test = read.csv('~/git/PD1-DREAM-Pred1ctors/CM_026_formatted_synthetic_data_subset/GRCh37ERCC_refseq105_genes_count.csv', header = T, row.names = 1)
model = readRDS('models/pfs_train_objects.rds')


# ======================================================================================

# Load our models
# ===========================================
# model = readRDS('/models/pfs_train_objects.rds')
# Clinical
clin.train = model$genes[,570:573]
# Expression
counts.train = model$genes[,1:569]
# Genes
genes = names(counts.train)


# Read data test
# ===========================================
# Load RefSeq105 genes raw counts
# Transposed
counts.test = as.data.frame(t(counts.test))
# Match with our genes
counts.test = counts.test[, match(genes, names(counts.test))]
# Apply log2 and scaling
counts.test = standarize(counts.test, log2 = T, scale = T)
# Impute NA´s
print('impute NAs')

res = impute.knn(t(counts.test))
counts.test = as.data.frame(t(res$data))


# Combat
# ===========================================
# Create metadata (with patID and subset variables)
meta = data.frame(patID = c(rownames(counts.train), rownames(counts.test)),
                  subset = c(rep('train', nrow(counts.train)), rep('test', nrow(counts.test))))
# Change rownames in meta in order to match with colnames in counts
rownames(meta) = meta$patID
# Bind rows between train and test counts
counts = rbind.data.frame(counts.train, counts.test)
# Create the model matrix from meta data
m = model.matrix(~1, data = meta)
# Run combat stablishing ref.batch=train
combat = ComBat(dat = t(counts), batch = meta$subset, mod = m, ref.batch = 'train')
# Tranposed combat matrix
combat = as.data.frame(t(combat))
# Retrieve test data counts from combat
print('combat')
counts.test = combat[match(rownames(counts.test), rownames(combat)), ]


# PCA
# ===========================================
# Load our PCA model
pca = model$modelPCA
# Change pca class to prcomp class to make the predict
pca <- list(sdev = pca$sdev,
                 rotation = data.matrix(pca$loadings),
                 x = data.matrix(pca$rotated),
                 center = TRUE, scale = TRUE)
class(pca) <- 'prcomp'
# Make predict
test.pca = predict(pca, newdata = counts.test)
test.pca = as.data.frame(test.pca)

print('normalizing by PCA')

# Add clinical data
# ===========================================
test.pca$Drug = 0
test.pca$Cohort = 1.5


# Predict biospear model
# ===========================================
biospear = readRDS("models/pfs_biospear_model_melanoma_Treatment.rds") #cambiar ruta para correr el docker!!!
print('biospear')

res = predBiospear(res = biospear, method = 'alassoR', newdata = test.pca, randomPFS = T)

predictions = cbind.data.frame(patientID = rownames(res), prediction = res[,1])

# write out inflammation signature to prediciton file
write.csv(predictions, file = args[2], quote = F, row.names = F); print("done writing out signature")



