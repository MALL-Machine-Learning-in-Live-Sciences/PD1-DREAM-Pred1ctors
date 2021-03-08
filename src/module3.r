# Module 3 (TMB)
# =============
# =============
require(survival)
require(rms)

# PFS
# ============================
data = readRDS('~/projects/Anti-PD1/data/metacohort/metacohort.rds')
clinical = data$clinical

# Filter patients
tmb = clinical[-which(is.na(clinical$TMB)),]
cc = tmb
cc = cc[which(cc$StageIV == 'yes'),]
cc = cc[which(cc$tumor.type == 'Melanoma'),]
print(table(cc$tumor.type, cc$Drug))
cc$TMB_c = ifelse(cc$TMB <= 100, 'low', 
                  ifelse(cc$TMB < 242 & cc$TMB > 100, 'medium', 'high'))
cc$TMB_h = ifelse(cc$TMB_c == 'high', 'yes', 'no')

cc$PFS.Event = as.numeric(cc$PFS.Event)
cc = cc[-which(is.na(cc$PFS)),]

cc$Gender[which(cc$Gender == 'female')] = 'F'
cc$Gender[which(cc$Gender == 'male')] = 'M'

cox = coxph(Surv(PFS, PFS.Event) ~ Age + Gender + TMB_h*Drug, data = cc, x = T, y = T)
summary(cox)

saveRDS(list(model = cox, train = cc), file = '~/projects/Anti-PD1/models_v2/pfs_model3.rds')


# OS
# ============================
data = readRDS('~/projects/Anti-PD1/data/metacohort/metacohort.rds')
clinical = data$clinical

# Filter patients
tmb = clinical[-which(is.na(clinical$TMB)),]
cc = tmb
cc = cc[which(cc$StageIV == 'yes'),]
print(table(cc$tumor.type, cc$Drug))

cc$TMB_c = ifelse(cc$TMB <= 100, 'low', 
                  ifelse(cc$TMB < 242 & cc$TMB > 100, 'medium', 'high'))
cc$TMB_h = ifelse(cc$TMB_c == 'high', 'yes', 'no')

cc$OS.Event = as.numeric(cc$OS.Event)
cc = cc[-which(is.na(cc$OS)),]

cox = coxph(Surv(OS, OS.Event) ~ Age + TMB_h*Drug, data = cc, x = T, y = T)
summary(cox)

saveRDS(list(model = cox, train = cc), file = '~/projects/Anti-PD1/models_v2/os_model3.rds')

