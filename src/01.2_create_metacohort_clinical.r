# Processing clinical data: Selection and labeling of patients by cohort
# Requisitos para añadir cohort:
#   1. Pacientes que tengan disponible datos de expresión antes del tratamientos (se eliminan aquellas muestras secuenciadas post treatment)
#   2. Pacientes que hayan sido tratados con algún tipo de fármaco anti-PD1
#   3. Pacientes que tengan datos de OS, PFS o Response (al menos uno de ellos).

require(data.table)
require(readxl)
dirData = '~/projects/Anti-PD1/data/'
setwd(dirData)


vars = c('patID', 'Cohort', 'tumor.type', 'Treatment', 'OS', 'OS.Event', 'PFS', 'PFS.Event',
         'Response', 'RECIST', 'Age', 'Gender', 'Smoking', 'stage', 'TMB')
res = list()

# Gide 2019
# ==========================
gide2019 = readRDS('processed/Gide2019.rds')
gide2019.c = cbind.data.frame(gide2019$antiPD1$cvrts, gide2019$antiPD1$surv)
gide2019.c1 = cbind.data.frame(gide2019$antiPD1.antiCTLA4$cvrts, gide2019$antiPD1.antiCTLA4$surv)
gide2019.clinical = rbind.data.frame(gide2019.c, gide2019.c1)
gide2019.clinical$patID = paste0('X', rownames(gide2019.clinical))
gide2019.clinical$tumor.type = 'Melanoma'
gide2019.clinical$Cohort = 'Gide2019'
gide2019.clinical$RECIST = gide2019.clinical$Best.RECIST.response
gide2019.clinical$Smoking = NA
gide2019.clinical$TMB = NA
gide2019.clinical$stage = 'stage iv'


clinical = gide2019.clinical 
clinical = clinical[, match(vars, names(clinical))]
clinical$Gender = ifelse(clinical$Gender == 0, 'male', 'female')

res$Gide2019 = clinical
rm(list = setdiff(ls(), c('res', 'dirData', 'vars')))


# GSE52562
# ==========================
gse52562 = readRDS('train/GSE52562/GSE52562_clinical.rds')
gse52562 = gse52562[grep('_pre', gse52562$source_name_ch1),]
os = NA
os.event = NA
pfs = gse52562$`pfs.days:ch1`
pfs.event = gse52562$`pfs.status.censorship:ch1`

clinical = data.frame(patID = gse52562$geo_accession,
           Cohort = 'GSE52562',
           tumor.type = 'Lymphoma',
           Treatment = 'Pidilizumab',
           OS = os,
           OS.Event = os.event,
           PFS = pfs,
           PFS.Event = pfs.event,
           Response = NA,
           RECIST = NA,
           Age = as.numeric(gse52562$`age:ch1`),
           Gender = gse52562$`gender:ch1`,
           Smoking = NA,
           stage = NA,
           TMB = NA)

clinical = clinical[, match(vars, names(clinical))]

res$GSE52562 = clinical
rm(list = setdiff(ls(), c('res', 'dirData', 'vars')))


# GSE67501
# ==========================
gse67501 = readRDS('train/GSE67501/GSE67501_clinical.rds')
response = gse67501$`response to anti-pd-1 (nivolumab) immunotherapy (response or no-response):ch1`
response = ifelse(response == 'response', 1, 0)

clinical = data.frame(
  patID = gse67501$geo_accession,
  Cohort = 'GSE67501',
  tumor.type = 'Renal Cell Carcinoma',
  Treatment = 'Nivolumab',
  OS = NA,
  OS.Event = NA,
  PFS = NA,
  PFS.Event = NA,
  Response = response,
  RECIST = gse67501$`response to anti-pd-1 (nivolumab) immunotherapy (complete response (cr) or partial response (pr) or stable disease (sd) or no response (nr)):ch1`,
  Age = NA,
  Gender = gse67501$`gender:ch1`,
  Smoking = NA,
  stage = gse67501$`primary tumor or metastasis:ch1`, #see supp material in https://cancerimmunolres.aacrjournals.org/content/suppl/2016/07/20/2326-6066.CIR-16-0072.DC1
  TMB = NA
)

clinical = clinical[, match(vars, names(clinical))]

res$GSE67501 = clinical
rm(list = setdiff(ls(), c('res', 'dirData', 'vars')))


# GSE78220 (Hugo 2016)
# ==========================
gse78220 = readRDS('train/GSE78220/GSE78220_clinical.rds')
gse78220 = gse78220[grep('pre-', gse78220$`biopsy time:ch1`),]

os = gse78220$`overall survival (days):ch1`
os.event = gse78220$`vital status:ch1`
response = ''
RECIST = gse78220$`anti-pd-1 response:ch1`

clinical = data.frame(
  patID = gse78220$title,
  Cohort = 'GSE78220',
  tumor.type = 'Melanoma',
  Treatment = 'Pembrolizumab',
  OS = os,
  OS.Event = os.event,
  PFS = NA,
  PFS.Event = NA,
  Response = NA,
  RECIST = RECIST,
  Age = as.numeric(gse78220$`age (yrs):ch1`),
  Gender = gse78220$`gender:ch1`,
  Smoking = NA,
  stage = gse78220$`disease status:ch1`,
  TMB = NA
)

clinical = clinical[, match(vars, names(clinical))]

res$GSE78220 = clinical
rm(list = setdiff(ls(), c('res', 'dirData', 'vars')))


# GSE115821 - for treatment response
# ==========================
gse115821 = readRDS('train/GSE115821/GSE115821_clinical.rds')
gse115821 = gse115821[grep('PRE', gse115821$`treatment state:ch1`), ]

clinical = data.frame(
  patID = gsub('.bam', '',  make.names(gse115821$title)),
  Cohort = 'GSE115821',
  tumor.type = 'Melanoma',
  Treatment = gse115821$`antibody:ch1`,
  OS = NA,
  OS.Event = NA,
  PFS = NA,
  PFS.Event = NA,
  Response = gse115821$`response:ch1`,
  RECIST = NA,
  Age = as.numeric(gse115821$`age at the baseline:ch1`),
  Gender = NA,
  Smoking = NA,
  stage = gse115821$`tumor type:ch1`,
  TMB = NA
)

clinical = clinical[, match(vars, names(clinical))]

res$GSE115821 = clinical
rm(list = setdiff(ls(), c('res', 'dirData', 'vars')))


# GSE100797 (Lauss et al 2017)
# ==========================
lauss2017 = readRDS('processed/Lauss2017.rds')
lauss2017.c = cbind.data.frame(lauss2017$cvrts, lauss2017$surv, lauss2017$meta)

pID = list()
ids = lauss2017.c$title
for (i in seq_along(ids)) {
  pID[[i]] = paste0('p', strsplit(ids[i], '_')[[1]][2])
}
pID = unlist(pID)

clinical = data.frame(
  patID = pID,
  Cohort = 'Lauss2017',
  tumor.type = 'Melanoma',
  Treatment = 'see paper',
  OS = lauss2017.c$`os.time:ch1`,
  OS.Event = lauss2017.c$`os.event:ch1`,
  PFS = lauss2017.c$`pfs.time:ch1`,
  PFS.Event = lauss2017.c$`pfs.event:ch1`,
  Response = NA,
  RECIST = lauss2017.c$`recist:ch1`,
  Age = NA,
  Gender = NA,
  Smoking = NA,
  stage = lauss2017.c$`ajcc.stage:ch1`,
  TMB = NA
)

clinical = clinical[, match(vars, names(clinical))]
clinical$OS = as.numeric(clinical$OS) * 30
clinical$PFS = as.numeric(clinical$PFS) * 30

res$Lauss2017 = clinical
rm(list = setdiff(ls(), c('res', 'dirData', 'vars')))


# Prat2017 (GSE93157)
# ==========================
prat2017 = readRDS('processed/Prat2017.rds')
prat2017 = cbind.data.frame(prat2017$cvrts, prat2017$surv, prat2017$meta)

clinical = data.frame(
  patID = rownames(prat2017),
  Cohort = 'Prat2017',
  tumor.type = prat2017$source_name_ch1,
  Treatment = prat2017$`drug:ch1`,
  OS = NA,
  OS.Event = NA,
  PFS = prat2017$`pfs:ch1`,
  PFS.Event = prat2017$`pfse:ch1`,
  Response = prat2017$`response:ch1`,
  RECIST = prat2017$`best.resp:ch1`,
  Age = as.numeric(prat2017$`age:ch1`),
  Gender = prat2017$`Sex:ch1`,
  Smoking = prat2017$`smoking:ch1`,
  stage = 'stage iv',
  TMB = NA
)

clinical = clinical[, match(vars, names(clinical))]
clinical$PFS = as.numeric(clinical$PFS) * 30

res$Prat2017 = clinical
rm(list = setdiff(ls(), c('res', 'dirData', 'vars')))


# Riaz2017 (GSE91061)
# ==========================
riaz = readRDS('train/Riaz2017/GSE91061_clinical_pre.rds')
clin.sup = as.data.frame(readxl::read_excel('~/projects/Anti-PD1/data/train/Riaz2017/clinical_suplementary_paper.xlsx', skip = 2))
clin.sup = clin.sup[match(riaz$title1, clin.sup$Patient),]
clinical = data.frame(
  patID = riaz$title,
  Cohort = 'Riaz2017',
  tumor.type = 'Melanoma',
  Treatment = 'Nivolumab',
  OS = riaz$OS,
  OS.Event = riaz$OS.Event,
  PFS = riaz$PFS,
  PFS.Event = riaz$PFS.Event,
  Response = riaz$Response,
  RECIST = riaz$RECIST,
  Age = as.numeric(riaz$Age),
  Gender = riaz$Gender,
  Smoking = NA,
  stage = clin.sup$`M Stage`,
  TMB = riaz$Mutation
)

clinical = clinical[, match(vars, names(clinical))]
clinical$Gender = ifelse(clinical$Gender == 0, 'male', 'female')

res$Riaz2017 = clinical
rm(list = setdiff(ls(), c('res', 'dirData', 'vars')))


# TCGA-SKCM
# ==========================
skcm = readRDS('TCGA/TCGA-SKCM-RNASeq_counts.rds')
skcm.c = skcm@colData
skcm.c = skcm.c[which(skcm.c$sample_type_id == '06' | skcm.c$sample_type_id == '01'),]
clinTCGA = read_excel('TCGA/clinical.xlsx')
clinTCGA = clinTCGA[which(clinTCGA$type == 'SKCM'),]

pats = skcm.c$patient
skcm.c = skcm.c[match(pats, skcm.c$patient),]
clinTCGA = clinTCGA[match(pats,clinTCGA$bcr_patient_barcode),]

tcga.tmb = read.delim2('~/git/PD1-DREAM-Pred1ctors/clinic_TCGA_TMB.tsv', header = T, sep = '\t')

clinical = data.frame(
  patID = skcm.c$barcode,
  Cohort = 'TCGA-SKCM',
  tumor.type = 'Melanoma',
  Treatment = 'Chemo',
  OS = clinTCGA$OS.time,
  OS.Event = clinTCGA$OS,
  PFS = clinTCGA$PFI.time,
  PFS.Event = clinTCGA$PFI,
  Response = NA,
  RECIST = NA,
  Age = skcm.c$age_at_diagnosis / 365,
  Gender = skcm.c$gender,
  Smoking = NA,
  stage = skcm.c$tumor_stage,
  TMB = tcga.tmb[match(skcm.c$patient, tcga.tmb$bcr_patient_barcode),]$TMB
)

clinical = clinical[, match(vars, names(clinical))]
res$TCGA.SKCM = clinical

rm(list = setdiff(ls(), c('res', 'dirData', 'vars')))


# TCGA-LUSC
# ==========================
lusc = readRDS('TCGA/TCGA-LUSC-RNASeq_counts.rds')
lusc.c = lusc@colData
lusc.c = lusc.c[which(lusc.c$tumor_stage == 'stage iii' |
                        lusc.c$tumor_stage == 'stage iiia' |
                        lusc.c$tumor_stage == 'stage iiib' |
                        lusc.c$tumor_stage == 'stage iv'),]
clinTCGA = read_excel('TCGA/clinical.xlsx')
clinTCGA = clinTCGA[which(clinTCGA$type == 'LUSC'),]

pats = lusc.c$patient
lusc.c = lusc.c[match(pats, lusc.c$patient),]
clinTCGA = clinTCGA[match(pats,clinTCGA$bcr_patient_barcode),]

tcga.tmb = read.delim2('~/git/PD1-DREAM-Pred1ctors/clinic_TCGA_TMB.tsv', header = T, sep = '\t')

clinical = data.frame(
  patID = lusc.c$barcode,
  Cohort = 'TCGA-LUSC',
  tumor.type = 'SQUAMOUS LUNG CANCER',
  Treatment = 'Chemo',
  OS = clinTCGA$OS.time,
  OS.Event = clinTCGA$OS,
  PFS = clinTCGA$PFI.time,
  PFS.Event = clinTCGA$PFI,
  Response = NA,
  RECIST = NA,
  Age = lusc.c$age_at_diagnosis / 365,
  Gender = lusc.c$gender,
  Smoking = NA,
  stage = lusc.c$tumor_stage,
  TMB = tcga.tmb[match(lusc.c$patient, tcga.tmb$bcr_patient_barcode),]$TMB
)

clinical = clinical[, match(vars, names(clinical))]
res$TCGA.LUSC = clinical

rm(list = setdiff(ls(), c('res', 'dirData', 'vars')))


# TCGA-LUAD
# ==========================
luad = readRDS('TCGA/TCGA-LUAD-RNASeq_counts.rds')
luad.c = luad@colData
luad.c = luad.c[which(luad.c$tumor_stage == 'stage iiia' |
                        luad.c$tumor_stage == 'stage iiib' |
                        luad.c$tumor_stage == 'stage iv'),]
clinTCGA = read_excel('TCGA/clinical.xlsx')
clinTCGA = clinTCGA[which(clinTCGA$type == 'LUAD'),]

pats = luad.c$patient
luad.c = luad.c[match(pats, luad.c$patient),]
clinTCGA = clinTCGA[match(pats,clinTCGA$bcr_patient_barcode),]

tcga.tmb = read.delim2('~/git/PD1-DREAM-Pred1ctors/clinic_TCGA_TMB.tsv', header = T, sep = '\t')

clinical = data.frame(
  patID = luad.c$barcode,
  Cohort = 'TCGA-LUAD',
  tumor.type = 'LUNG NON-SQUAMOUS CANCER',
  Treatment = 'Chemo',
  OS = clinTCGA$OS.time,
  OS.Event = clinTCGA$OS,
  PFS = clinTCGA$PFI.time,
  PFS.Event = clinTCGA$PFI,
  Response = NA,
  RECIST = NA,
  Age = luad.c$age_at_diagnosis / 365,
  Gender = luad.c$gender,
  Smoking = NA,
  stage = luad.c$tumor_stage,
  TMB = tcga.tmb[match(luad.c$patient, tcga.tmb$bcr_patient_barcode),]$TMB
)

clinical = clinical[, match(vars, names(clinical))]
res$TCGA.LUAD = clinical

rm(list = setdiff(ls(), c('res', 'dirData', 'vars')))





# Clinical Data of metacohort
# =======================================================================
source('~/git/Anti-PD1/utils.r')
res = rbindlist(res)

res$tumor.type = replace(res$tumor.type,
                         res$tumor.type == 'MELANOMA',
                         'Melanoma')

res$Treatment = replace(res$Treatment, 
                        res$Treatment == 'NIVOLUMAB',
                        'Nivolumab')

res$Treatment = replace(res$Treatment, 
                        res$Treatment == 'PEMBROLIZUMAB',
                        'Pembrolizumab')

res$Drug = ifelse(res$Treatment == 'anti-PD-1' |
                    res$Treatment == 'anti-PD-1+anti-CTLA-4' |
                    res$Treatment == 'Ipilimumab + Nivolumab' |
                    res$Treatment == 'Ipilimumab + Pembrolizumab' |
                    res$Treatment == 'Nivolumab' |
                    res$Treatment == 'Pembrolizumab' |
                    res$Treatment == 'Pidilizumab', 'anti-PD1',
                  ifelse(res$Treatment == 'Chemo', 'Chemo', 
                         ifelse(res$Treatment == 'see paper', 'see paper', 'no anti-PD-1')))

res$OS = as.numeric(res$OS)
res$PFS = as.numeric(res$PFS)

res$OS.Event = replace(res$OS.Event,
                       res$OS.Event == 'Alive',
                       0)
res$OS.Event = replace(res$OS.Event,
                       res$OS.Event == 'Dead',
                       1)

# res = labelData(5, res)

res$Age = round(res$Age, 0)

res$Gender = replace(res$Gender,
                     res$Gender == 'F',
                     'female')
res$Gender = replace(res$Gender,
                     res$Gender == 'M',
                     'male')

res$StageIV = ifelse(res$stage == 'M1a' | res$stage == 'M1a' | res$stage == 'M1A' | res$stage == 'M1b' | res$stage == 'M1B' | res$stage == 'M1c' | res$stage == 'M1C' |
         res$stage == 'stage iv' | res$stage == 'Metastasis' | res$stage == 'metastatic melanoma', 'yes', 'no')


saveRDS(res, file = 'metacohort/clinical_metacohort2.rds')
