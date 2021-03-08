# Create metacohort counts
# =================================
source('~/git/Anti-PD1/utils.r')
require(data.table)
require(GEOquery)

dirData = '~/projects/Anti-PD1/data/'
setwd(dirData)

clinical = readRDS('metacohort/clinical_metacohort.rds')

# Gide 2019
# ==========================
gide2019 = readRDS('processed/Gide2019.rds')
gide2019 = rbind.data.frame(gide2019$antiPD1$counts, gide2019$antiPD1.antiCTLA4$counts)

pats = clinical[which(clinical$Cohort == 'Gide2019'), ]$patID

gide2019 = gide2019[match(intersect(pats, rownames(gide2019)), rownames(gide2019)), ]
names(gide2019) = convertGeneOnto(names(gide2019), from = 'ENTREZID', to = 'SYMBOL')


# GSE52562
# ==========================
gse52562 = readRDS('train/GSE52562/GSE52562_counts.rds')

pats = clinical[which(clinical$Cohort == 'GSE52562'), ]$patID

gse52562 = gse52562[match(intersect(pats, rownames(gse52562)), rownames(gse52562)), ]

annot = getGEO('GPL10558')
annot = annot@dataTable@table
gse52562 = selectSonda(gse52562, annot)


# GSE67501
# ==========================
gse67501 = readRDS('train/GSE67501/GSE67501_counts.rds')

pats = clinical[which(clinical$Cohort == 'GSE67501'), ]$patID

gse67501 = gse67501[match(intersect(pats, rownames(gse67501)), rownames(gse67501)), ]

annot = getGEO('GPL14951')
annot = annot@dataTable@table
gse67501 = selectSonda(gse67501, annot)


# GSE78220 (Hugo 2016)
# ==========================
gse78220 = readRDS('train/GSE78220/GSE78220_counts.rds')
names(gse78220) = gse78220[1,]; gse78220 = gse78220[-1,]

rwn = rownames(gse78220)
patsID = list()
for (i in seq_along(rwn)) {
  patsID[[i]] = strsplit(rwn[i], '\\.')[[1]][1]
}
patsID = unlist(patsID)
rownames(gse78220) = patsID

pats = clinical[which(clinical$Cohort == 'GSE78220'), ]$patID

gse78220 = gse78220[match(intersect(pats, rownames(gse78220)), rownames(gse78220)), ]


# GSE115821
# ==========================
gse115821 = readRDS('train/GSE115821/GSE115821_counts.rds')

pats = clinical[which(clinical$Cohort == 'GSE115821'), ]$patID

rwn = gsub('.bam', '', rownames(gse115821))
rownames(gse115821) = rwn
gse115821 = gse115821[match(intersect(pats, rownames(gse115821)), rownames(gse115821)), ]


# GSE100797 (Lauss et al 2017)
# ==========================
lauss2017 = readRDS('processed/Lauss2017.rds')
lauss2017 = lauss2017$counts

rwn = rownames(lauss2017)
patsID = list()
for (i in seq_along(rwn)) {
  patsID[[i]] = strsplit(rwn[i], '_')[[1]][2]
}
patsID = paste0('p', unlist(patsID))
rownames(lauss2017) = patsID

pats = clinical[which(clinical$Cohort == 'Lauss2017'), ]$patID

lauss2017 = lauss2017[match(intersect(pats, rownames(lauss2017)), rownames(lauss2017)), ]


# Prat2017 (GSE93157)
# ==========================
prat2017 = readRDS('processed/Prat2017.rds')
prat2017 = prat2017$counts

pats = clinical[which(clinical$Cohort == 'Prat2017'), ]$patID

prat2017 = prat2017[match(intersect(pats, rownames(prat2017)), rownames(prat2017)), ]


# Riaz2017 (GSE91061)
# ==========================
riaz2017 = readRDS('train/Riaz2017/GSE91061_counts_pre.rds')

pats = clinical[which(clinical$Cohort == 'Riaz2017'), ]$patID

riaz2017 = riaz2017[match(intersect(pats, rownames(riaz2017)), rownames(riaz2017)), ]
names(riaz2017) = convertGeneOnto(names(riaz2017), from = 'ENTREZID', to = 'SYMBOL')


# TCGA-SKCM
# ==========================
skcm = readRDS('TCGA/TCGA-SKCM-RNASeq_counts.rds')
skcm.o = as.data.frame(t(assay(skcm)))
names(skcm.o) = skcm@rowRanges$external_gene_name

pats = clinical[which(clinical$Cohort == 'TCGA-SKCM'), ]$patID

skcm.o = skcm.o[match(intersect(pats, rownames(skcm.o)), rownames(skcm.o)), ]


# TCGA-LUSC
# ==========================
lusc = readRDS('TCGA/TCGA-LUSC-RNASeq_counts.rds')
lusc.o = as.data.frame(t(assay(lusc)))
names(lusc.o) = lusc@rowRanges$external_gene_name

pats = clinical[which(clinical$Cohort == 'TCGA-LUSC'), ]$patID

lusc.o = lusc.o[match(intersect(pats, rownames(lusc.o)), rownames(lusc.o)), ]


# TCGA-LUAD
# ==========================
luad = readRDS('TCGA/TCGA-LUAD-RNASeq_counts.rds')
luad.o = as.data.frame(t(assay(luad)))
names(luad.o) = luad@rowRanges$external_gene_name

pats = clinical[which(clinical$Cohort == 'TCGA-LUAD'), ]$patID

luad.o = luad.o[match(intersect(pats, rownames(luad.o)), rownames(luad.o)), ]




# Synthetic data
# ========================================================
synthetic = read.csv('synthetic/GRCh37ERCC_refseq105_genes_count.csv', header = T, row.names = 1)
synthetic = as.data.frame(t(synthetic))



rm(list = setdiff(ls(), c('clinical', 'gide2019', 'gse52562', 'gse67501',
                          'gse78220', 'gse115821', 'lauss2017', 'prat2017', 
                          'riaz2017', 'synthetic', 'standarize',
                          'skcm.o', 'lusc.o', 'luad.o')))


# Common genes among datasets and synthetic data
commons = intersect(names(gide2019),
                    intersect(names(gse52562),
                              intersect(names(gse67501),
                                        intersect(names(gse78220),
                                                  intersect(names(lauss2017),
                                                            intersect(names(prat2017),
                                                                      intersect(names(riaz2017),
                                                                                intersect(names(gse115821),
                                                                                          intersect(names(synthetic),
                                                                                                    intersect(names(skcm.o),
                                                                                                              intersect(names(lusc.o),
                                                                                                                        names(luad.o))))))))))))

# Matching according synthetic data annotation
gide2019 = gide2019[, match(commons, names(gide2019))]
gse52562 = gse52562[, match(commons, names(gse52562))]
gse67501 = gse67501[, match(commons, names(gse67501))]
gse78220 = gse78220[, match(commons, names(gse78220))]
lauss2017 = lauss2017[, match(commons, names(lauss2017))]
prat2017 = prat2017[, match(commons, names(prat2017))]
riaz2017 = riaz2017[, match(commons, names(riaz2017))]
gse115821 = gse115821[, match(commons, names(gse115821))]
skcm.o = skcm.o[,match(commons, names(skcm.o))]
luad.o = luad.o[,match(commons, names(luad.o))]
lusc.o = lusc.o[,match(commons, names(lusc.o))]

metaCohort = rbind.data.frame(gide2019,
                 gse52562,
                 gse67501,
                 gse78220,
                 lauss2017,
                 prat2017,
                 riaz2017,
                 gse115821,
                 skcm.o,
                 lusc.o,
                 luad.o)

rwn = rownames(metaCohort)
metaCohort = as.data.frame(sapply(metaCohort, as.numeric))
rownames(metaCohort) = rwn

clinical = clinical[match(rownames(metaCohort), clinical$patID),]

# Standarize cohorts (if necessary)
table(clinical$Cohort)

# Gide2019
pats = clinical[which(clinical$Cohort == 'Gide2019'),]$patID
gide2019 = metaCohort[match(pats, rownames(metaCohort)),]

# GSE115821
pats = clinical[which(clinical$Cohort == 'GSE115821'),]$patID
gse115821 = metaCohort[match(pats, rownames(metaCohort)),]
gse115821 = standarize(gse115821, log2 = T, scale = T)

# GSE52562
pats = clinical[which(clinical$Cohort == 'GSE52562'),]$patID
gse52562 = metaCohort[match(pats, rownames(metaCohort)),]
gse52562 = standarize(gse52562, log2 = T, scale = T)

# GSE67501
pats = clinical[which(clinical$Cohort == 'GSE67501'),]$patID
gse67501 = metaCohort[match(pats, rownames(metaCohort)),]
gse67501 = standarize(gse67501, log2 = F, scale = T)

# GSE78220
pats = clinical[which(clinical$Cohort == 'GSE78220'),]$patID
gse78220 = metaCohort[match(pats, rownames(metaCohort)),]
gse78220 = standarize(gse78220, log2 = F, scale = T)

# Lauss2017
pats = clinical[which(clinical$Cohort == 'Lauss2017'),]$patID
lauss2017 = metaCohort[match(pats, rownames(metaCohort)),]
lauss2017 = standarize(lauss2017, log2 = F, scale = T)

# Prat2017
pats = clinical[which(clinical$Cohort == 'Prat2017'),]$patID
prat2017 = metaCohort[match(pats, rownames(metaCohort)),]
prat2017 = standarize(prat2017, log2 = F, scale = T)

# Riaz2017
pats = clinical[which(clinical$Cohort == 'Riaz2017'),]$patID
riaz2017 = metaCohort[match(pats, rownames(metaCohort)),]
riaz2017 = standarize(riaz2017, log2 = F, scale = T)

# SKCM
pats = clinical[which(clinical$Cohort == 'TCGA-SKCM'),]$patID
skcm.o = metaCohort[match(pats, rownames(metaCohort)),]
skcm.o = standarize(skcm.o, log2 = T, scale = T)

# LUSC
pats = clinical[which(clinical$Cohort == 'TCGA-LUSC'),]$patID
lusc.o = metaCohort[match(pats, rownames(metaCohort)),]
lusc.o = standarize(lusc.o, log2 = T, scale = T)

# LUAD
pats = clinical[which(clinical$Cohort == 'TCGA-LUAD'),]$patID
luad.o = metaCohort[match(pats, rownames(metaCohort)),]
luad.o = standarize(luad.o, log2 = T, scale = T)


metaCohort2 = rbind.data.frame(gide2019,
                              gse52562,
                              gse67501,
                              gse78220,
                              lauss2017,
                              prat2017,
                              riaz2017,
                              gse115821,
                              skcm.o,
                              lusc.o,
                              luad.o)

counts = metaCohort2[ , colSums(is.na(metaCohort2)) == 0]

clinical = clinical[match(rownames(counts), clinical$patID), ]
saveRDS(list(clinical = clinical, counts = counts), file = 'metacohort/metacohort.rds')

# Check!
stopifnot(rownames(counts) == clinical$patID)


