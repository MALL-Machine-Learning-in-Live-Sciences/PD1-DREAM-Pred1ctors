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


supressMessages(library(sva))
supressMessages(library(PCAtools))

# get command line args: args[1] input rna-seq gene level count data, args[2] output file: two column(patient id, inflammation score) .csv file
args        <- commandArgs(trailingOnly=TRUE)
counts      <- fread(args[1],data.table = F); rownames(counts) <- counts[,1]; counts <- counts[,-1]; print("done reading in counts")


# Data directory of synthetic data
# ===========================================
dirData = '~/git/PD1-DREAM-Pred1ctors/CM_026_formatted_synthetic_data_subset/'
setwd(dirData)


# Load train objects
# ===========================================
model = readRDS('../../../projects/Anti-PD1/data/metacohort/pfs_train_objects.rds')

# Clinical
clin.train = model$genes[,570:573]
# Expression
counts.train = model$genes[,1:569]
# Genes
genes = names(counts.train)


# Read data test
# ===========================================
# Load RefSeq105 genes raw counts
counts.test = read.csv('GRCh37ERCC_refseq105_genes_count.csv', header = T, row.names = 1)
# Transposed
counts.test = as.data.frame(t(counts.test))
# Match with our genes
counts.test = counts.test[, match(genes, names(counts.test))]
# Apply log2 and scaling
counts.test = standarize(counts.test, log2 = T, scale = T)


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
counts.test = combat[match(rownames(counts.test), rownames(combat)), ]


# Prediction of treatment variable
# ===========================================
# ???????


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


# Add clinical data
# ===========================================



# write out inflammation signature to prediciton file
write.csv(inflam_sig, file = args[2], quote = F, row.names = F); print("done writing out signature")

rm(tmm,counts)

