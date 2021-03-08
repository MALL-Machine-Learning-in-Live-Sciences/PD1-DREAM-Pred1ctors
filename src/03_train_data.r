# See PCA
require(PCAtools)
require(readxl)
require(sva)
require(SummarizedExperiment)
require(survival)
require(rms)
require(survminer)
source('~/git/Anti-PD1/utils.r')

data = readRDS('~/projects/Anti-PD1/data/metacohort/metacohort.rds')
clinical = data$clinical
counts = data$counts

# Select melanoma samples
pats.pd1 = clinical[which(clinical$tumor.type == 'Melanoma' & clinical$Drug == 'anti-PD1'),]$patID
pats.chemo = clinical[which(clinical$tumor.type == 'Melanoma' & clinical$Drug == 'Chemo'),]$patID

pats = c(pats.pd1, pats.chemo)

clinical = clinical[match(pats, clinical$patID),]
counts = counts[match(pats, rownames(counts)), ]

# Biplot pre-combat
rownames(clinical) = rownames(counts)
pca.pre = pca(t(counts), clinical, scale = F, center = F)
biplot(pca.pre, colby = 'Cohort', lab = NULL, legendPosition = 'right')

# ComBat
# =======================================
meta = clinical[, c('Cohort', 'patID')]
rownames(meta) = meta$patID

model = model.matrix(~1, data = meta)
combat = ComBat(dat = t(counts), batch = meta$Cohort, mod = model)
combat = as.data.frame(t(combat))

# Biplot post-combat
rownames(combat) = clinical$patID
pca = pca(t(combat), meta, scale = F, center = F)
biplot(pca, colby = 'Cohort', lab = NULL, legendPosition = 'right')

counts = combat


# Create datasets for OS and PFS
# =======================================
table(is.na(counts))

# OS
os = as.data.frame(clinical[, c('patID', 'OS', 'OS.Event')])
os = os[complete.cases(os),]

os.clin = clinical[match(os$patID, clinical$patID), ]
os.counts = counts[match(os$patID, rownames(counts)), ]

os.clin$OS.Event = as.numeric(os.clin$OS.Event)
saveRDS(list(clinical = os.clin, counts = os.counts), file = '~/projects/Anti-PD1/data/metacohort/overallSurvival_train.rds')



# PFS
pfs = as.data.frame(clinical[, c('patID', 'PFS', 'PFS.Event')])
pfs = pfs[complete.cases(pfs),]

pfs.clin = clinical[match(pfs$patID, clinical$patID), ]
pfs.counts = counts[match(pfs$patID, rownames(counts)), ]

pfs.clin$PFS.Event = as.numeric(pfs.clin$PFS.Event)
saveRDS(list(clinical = pfs.clin, counts = pfs.counts), file = '~/projects/Anti-PD1/data/metacohort/ProgressionFreeSurvival_train.rds')



