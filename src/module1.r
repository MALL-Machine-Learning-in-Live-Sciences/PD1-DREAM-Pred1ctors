# Module 1 (PCAs)
# =============
# =============
# PCAs calculation
require(PCAtools)
require(biospear)

# PFS
# =================================
pfs = readRDS('~/projects/Anti-PD1/data/metacohort/ProgressionFreeSurvival_train.rds')

meta = as.data.frame(pfs$clinical)
counts = pfs$counts
rownames(meta) = meta$patID

# Calculating PCA큦
pca = pca(counts, meta, removeVar = 0.1, transposed = T, center = F, scale = F)
elbow <- findElbowPoint(pca$variance); print(elbow)

# Creating train dataset (PCA큦 + cvrts)
cvrts = c('Cohort', 'Age', 'Drug', 'PFS', 'PFS.Event')
train = cbind.data.frame(meta[, cvrts], pca$rotated)

# Convert class to PCA object
pca <- list(sdev = pca$sdev,
            rotation = data.matrix(pca$loadings),
            x = data.matrix(pca$rotated),
            center = TRUE, scale = TRUE)
class(pca) <- 'prcomp'

# Select only stageIV patients
train.s4 = train[which(meta$StageIV == 'yes'),]

# Format covariates
train.s4$Cohort = as.numeric(as.factor(train.s4$Cohort)) - 1  # 0 = Gide2019;  1 = Prat2017;  2 = Riaz2017; 3 = TCGA-SKCM
train.s4$Drug = as.numeric(as.factor(train.s4$Drug)) - 1      # 0 = anti-PD1; 1 = Chemo

# Running biospear models
set.seed(1993)
model.pfs = BMsel(train.s4, paste0('PC', 1:15) , c("PFS","PFS.Event"), c("Cohort", "Age"), tt="Drug", inter=T, std.x = TRUE, std.i = FALSE, std.tt = F,
               method = c('alassoR', 'alassoU', 'enet',  'lasso', 'lasso-1se', 'lasso-AIC', 'lasso-BIC',
                          'lasso-HQIC', 'lasso-pct', 'lasso-pcvl','lasso-RIC', 'PCAlasso', 'PLSlasso', 'ridgelasso',  'uniFDR'),
               folds = 10, uni.fdr = 0.05, uni.test = 1, ss.rando = F, ss.nsub = 100,
               ss.fsub = 0.5, ss.fwer = 1, ss.thr = 0.6, dfmax = 70,
               pct.rep = 1, pct.qtl = 0.95, showWarn = TRUE, trace = TRUE)

res = list(pca = pca, train = train.s4, genes = counts, model = model.pfs)
saveRDS(res, file = '~/projects/Anti-PD1/models_v2/pfs_biospear.rds')









# OS
# =================================
os = readRDS('~/projects/Anti-PD1/data/metacohort/overallSurvival_train.rds')

meta = as.data.frame(os$clinical)
counts = os$counts
rownames(meta) = meta$patID

# Calculating PCA큦
pca = pca(counts, meta, removeVar = 0.1, transposed = T, center = F, scale = F)
elbow <- findElbowPoint(pca$variance); print(elbow)

# Creating train dataset (PCA큦 + cvrts)
cvrts = c('Cohort', 'Age', 'Drug', 'OS', 'OS.Event')
train = cbind.data.frame(meta[, cvrts], pca$rotated)

# Convert class to PCA object
pca <- list(sdev = pca$sdev,
            rotation = data.matrix(pca$loadings),
            x = data.matrix(pca$rotated),
            center = TRUE, scale = TRUE)
class(pca) <- 'prcomp'

# Select only stageIV patients
train.s4 = train[which(meta$StageIV == 'yes'),]

# Format covariates
train.s4$Cohort = as.numeric(as.factor(train.s4$Cohort)) - 1
train.s4$Drug = as.numeric(as.factor(train.s4$Drug)) - 1    # 0 = anti-PD1; 1 = Chemo

# Running biospear models
model.os = BMsel(train.s4, paste0('PC', 1:14) , c("OS","OS.Event"), c("Age"), tt="Drug", inter=T, std.x = TRUE, std.i = FALSE, std.tt = F,
                  method = c('alassoR', 'alassoU', 'enet',  'lasso', 'lasso-1se', 'lasso-AIC', 'lasso-BIC',
                             'lasso-HQIC', 'lasso-pct', 'lasso-pcvl','lasso-RIC', 'PCAlasso', 'PLSlasso', 'ridgelasso',  'uniFDR'),
                  folds = 10, uni.fdr = 0.05, uni.test = 1, ss.rando = F, ss.nsub = 100,
                  ss.fsub = 0.5, ss.fwer = 1, ss.thr = 0.6, dfmax = 70,
                  pct.rep = 1, pct.qtl = 0.95, showWarn = TRUE, trace = TRUE)

res = list(pca = pca, train = train.s4, genes = counts, model = model.os)
saveRDS(res, file = '~/projects/Anti-PD1/models_v2/os_biospear.rds')








pred = predRes(res = xx$model,
               method=c('alassoR', 'alassoU', 'enet',  'lasso', 'lasso-1se', 'lasso-AIC', 'lasso-BIC', 'lasso-HQIC',
                        'lasso-pct', 'lasso-pcvl','lasso-RIC', 'PCAlasso', 'PLSlasso', 'ridgelasso'),
               traindata = train.s4,
               int.cv=T,
               int.cv.nfold = 5,
               time=seq(100,1000,100),
               trace = TRUE,
               ncores = 5)
