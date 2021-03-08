require(readxl)
luad = readRDS('~/projects/Anti-PD1/data/TCGA/TCGA-LUAD-RNASeq_counts.rds')
lusc = readRDS('~/projects/Anti-PD1/data/TCGA/TCGA-LUSC-RNASeq_counts.rds')
clinTCGA = read_excel('~/projects/Anti-PD1/data/TCGA/clinical.xlsx')
clinTCGA = clinTCGA[which(clinTCGA$type == 'LUAD' | clinTCGA$type == 'LUSC'),]

luad = luad[match(luad$patient, luad$patient),]
lusc = lusc[match(lusc$patient, lusc$patient),]

luad.surv = clinTCGA[match(luad$patient, clinTCGA$bcr_patient_barcode),]
lusc.surv = clinTCGA[match(lusc$patient, clinTCGA$bcr_patient_barcode),]

luad.c = data.frame(
  age = round(luad$age_at_diagnosis / 365, 0),
  gender = luad$gender,
  tumorType = 'NON-SQUAMOUS',
  PFS = luad.surv$PFI.time,
  PFS.Event = luad.surv$PFI,
  OS = luad.surv$OS.time,
  OS.Event = luad.surv$OS,
  row.names = luad$barcode
)

lusc.c = data.frame(
  age = round(lusc$age_at_diagnosis / 365, 0),
  gender = lusc$gender,
  tumorType = 'SQUAMOUS',
  PFS = lusc.surv$PFI.time,
  PFS.Event = lusc.surv$PFI,
  OS = lusc.surv$OS.time,
  OS.Event = lusc.surv$OS,
  row.names = lusc$barcode
)

data = rbind.data.frame(luad.c, lusc.c)

data$gender[which(data$gender == 'female')] = 'F'
data$gender[which(data$gender == 'male')] = 'M'
ggsurvplot(survfit(Surv(PFS, PFS.Event) ~ tumorType, data = data), data = data, risk.table = T)

# PFS
# =================
cox = coxph(Surv(PFS, PFS.Event) ~ age + tumorType + gender, data = data, x = T)
summary(cox)
saveRDS(list(train = data, model = cox), file = '~/projects/Anti-PD1/models_v2/pfs_model2.rds')


# OS
# =================
cox = coxph(Surv(OS, OS.Event) ~ age + tumorType + gender, data = data, x = T)
summary(cox)
saveRDS(list(train = data, model = cox), file = '~/projects/Anti-PD1/models_v2/os_model2.rds')
