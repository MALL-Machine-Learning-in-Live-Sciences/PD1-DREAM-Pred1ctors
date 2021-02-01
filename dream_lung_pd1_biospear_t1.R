# biospear in antipd1

mel_data = readRDS("melanoma_PCAs.rds")


mel_data$Cohort_num = as.numeric(as.factor(mel_data$Cohort))-1
mel_data$Drug_num = as.numeric(as.factor(mel_data$Drug))-1
mel_data$PFS.Event_num = ifelse(mel_data$PFS.Event=="Progression",1,0)
mel_data$PFS[which(mel_data$PFS==0)]=1


library(biospear)


# con 15 , con 50
set.seed(1980)
model1 = BMsel(mel_data,6:20 , c("PFS","PFS.Event_num"), c("Cohort_num"),tt="Drug_num", inter=T, std.x = TRUE, std.i = FALSE, std.tt = F,
               #method = c('alassoL', 'alassoR', 'alassoU', 'enet', 'gboost', 'glasso', 'lasso', 'lasso-1se', 'lasso-AIC', 'lasso-BIC',
               #           'lasso-HQIC', 'lasso-pct', 'lasso-pcvl','lasso-RIC', 'modCov', 'PCAlasso', 'PLSlasso', 'ridge', 'ridgelasso', 'stabSel', 'uniFDR'),
               method = c('alassoR', 'alassoU', 'enet',  'lasso', 'lasso-1se', 'lasso-AIC', 'lasso-BIC',
                          'lasso-HQIC', 'lasso-pct', 'lasso-pcvl','lasso-RIC', 'PCAlasso', 'PLSlasso', 'ridge', 'ridgelasso',  'uniFDR'),
              #method = c('lasso',  'uniFDR'),
               folds = 10, uni.fdr = 0.05, uni.test = 1, ss.rando = F, ss.nsub = 100,
               ss.fsub = 0.5, ss.fwer = 1, ss.thr = 0.6, dfmax = 70,
               pct.rep = 1, pct.qtl = 0.95, showWarn = TRUE, trace = TRUE)

model2 = BMsel(mel_data,6:55 , c("PFS","PFS.Event_num"), c("Cohort_num"),tt="Drug_num", inter=T, std.x = TRUE, std.i = FALSE, std.tt = F,
               #method = c('alassoL', 'alassoR', 'alassoU', 'enet', 'gboost', 'glasso', 'lasso', 'lasso-1se', 'lasso-AIC', 'lasso-BIC',
               #           'lasso-HQIC', 'lasso-pct', 'lasso-pcvl','lasso-RIC', 'modCov', 'PCAlasso', 'PLSlasso', 'ridge', 'ridgelasso', 'stabSel', 'uniFDR'),
               method = c('alassoR', 'alassoU', 'enet',  'lasso', 'lasso-1se', 'lasso-AIC', 'lasso-BIC',
                          'lasso-HQIC', 'lasso-pct', 'lasso-pcvl','lasso-RIC', 'PCAlasso', 'PLSlasso', 'ridge', 'ridgelasso',  'uniFDR'),
               #method = c('lasso',  'uniFDR'),
               folds = 10, uni.fdr = 0.05, uni.test = 1, ss.rando = F, ss.nsub = 100,
               ss.fsub = 0.5, ss.fwer = 1, ss.thr = 0.6, dfmax = 70,
               pct.rep = 1, pct.qtl = 0.95, showWarn = TRUE, trace = TRUE)

save(model1,model2,file="models_melanoma.RData")


pred1 = predRes(model1, method=c('alassoR', 'alassoU', 'enet',  'lasso', 'lasso-1se', 'lasso-AIC', 'lasso-BIC', 'lasso-HQIC', 'lasso-pct', 'lasso-pcvl','lasso-RIC', 'PCAlasso', 'PLSlasso',  'ridgelasso',  'uniFDR'), 
                mel_data,  int.cv=T, int.cv.nfold = 5, time=seq(100,1500,100),
                trace = TRUE, ncores = 5)

pred1_50 = predRes(model2, method=c('alassoR', 'alassoU', 'enet',  'lasso', 'lasso-1se', 'lasso-AIC', 'lasso-BIC', 'lasso-HQIC', 'lasso-pct', 'lasso-pcvl','lasso-RIC', 'PCAlasso', 'PLSlasso',  'ridgelasso',  'uniFDR'), 
                mel_data,  int.cv=T, int.cv.nfold = 5, time=seq(100,1500,100),
                trace = TRUE, ncores = 20)

save(pred1,pred1_50,file="model_results_melanoma.RData")

pred2 = predRes(model1, method=c('lasso-1se', 'lasso-AIC'), 
                mel_data,  int.cv=T, int.cv.nfold = 5, time=c(500,800,1000,1200),
                trace = TRUE, ncores = 5)
