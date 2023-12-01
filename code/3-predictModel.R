# Predict model
source('./Function.R')
dt_score <- openxlsx::read.xlsx('./dt_score.xlsx')
use_dt_case <- openxlsx::read.xlsx('./use_SGA-627case-prognosis.xlsx',sheet='use_dt_case')
use_dt_ctr <- openxlsx::read.xlsx('./use_SGA-627case-prognosis.xlsx',sheet='use_dt_ctr')
dt_vali <- rio::import('./dt_vali.xlsx') 
dt <- plyr::rbind.fill(use_dt_case,use_dt_ctr)
dt <- merge(dt,dt_score,by.x='sample',by.y='Sample_name',all.x=T)
dt <- dt %>% mutate(score=ifelse(is.na(score),0,score))
dt_vali <- dt_vali %>% mutate(score=ifelse(is.na(score),0,score))
rownames(dt) <- dt$sample
rownames(dt_vali) <- dt_vali$sample
dt <- dt[,c("censor","gender","GA_week","BW","Central_nervous_system_abnormalities",     
            "Cardiovascular_abnormalities","Evidence_of_metabolic_disease",
            "Skeletal_abnormalities","Digestive_system_abnormalities",
            "Respiratory_system_abnormalities","Urinary_or_reproductive_system_abnormalities",
            "Hematologic_abnormalities","Infections_and_immune_abnormalities",
            "Craniofacial_anomalies","Skin",
            "score")]
dt_vali <- dt_vali[,c("censor","gender","GA_week","BW","Central_nervous_system_abnormalities",     
                      "Cardiovascular_abnormalities","Evidence_of_metabolic_disease",
                      "Skeletal_abnormalities","Digestive_system_abnormalities",
                      "Respiratory_system_abnormalities","Urinary_or_reproductive_system_abnormalities",
                      "Hematologic_abnormalities","Infections_and_immune_abnormalities",
                      "Craniofacial_anomalies","Skin",
                      "score")]
dt[is.na(dt)] <- 'N'
dt_vali[is.na(dt_vali)] <- 'N'
cont_var <- c('GA_week','BW','score')
cat_var <- setdiff(colnames(dt),c(cont_var,"censor"))
dt$gender <- factor(dt$gender,levels = c('Male','Female'))
dt_vali$gender <- factor(dt_vali$gender,levels = c('male','female'))
dt <- dt %>% mutate_at(cat_var,function(x)as.factor(x))
dt_vali <- dt_vali %>% mutate_at(cat_var,function(x)as.factor(x))


use_var <- c("censor",'GA_week','BW',"gender","score","Central_nervous_system_abnormalities",
             "Evidence_of_metabolic_disease","Skeletal_abnormalities",
             "Digestive_system_abnormalities","Respiratory_system_abnormalities",
             "Urinary_or_reproductive_system_abnormalities","Infections_and_immune_abnormalities",
             "Craniofacial_anomalies","Skin")
cat_var <- c("censor","gender","Central_nervous_system_abnormalities",
             "Evidence_of_metabolic_disease","Skeletal_abnormalities",
             "Digestive_system_abnormalities","Respiratory_system_abnormalities",
             "Urinary_or_reproductive_system_abnormalities","Infections_and_immune_abnormalities",
             "Craniofacial_anomalies","Skin")
cont_var <- c('GA_week','BW','score')
dt <- dt[,use_var]
dt_vali <- dt_vali[,use_var]
flag <- as.numeric(as.factor(dt$censor))
dt <- dt %>% mutate_at(cat_var,function(x)as.numeric(as.factor(x)))
dt_vali <- dt_vali %>% mutate_at(cat_var,function(x)as.numeric(as.factor(x)))
set.seed(42)
train_index <- sample(rownames(dt),0.7*nrow(dt))
test_index <- setdiff(rownames(dt),train_index)
dt_train <- dt[train_index,]
dt_test <- dt[test_index,]
dt_train_flag <- dt_train$censor
dt_test_flag <- dt_test$censor
dt_vali_flag <- dt_vali$censor
preProcValues <- preProcess(dt_train[,-1],method = c("center", "scale"))
dt_train <- predict(preProcValues, dt_train[,-1])
dt_test <- predict(preProcValues, dt_test[,-1])
dt_vali <- predict(preProcValues, dt_vali[,-1])
dt_train$flag <- dt_train_flag
dt_test$flag <- dt_test_flag
dt_vali$flag <- dt_vali_flag

nfolds <- 10
auc_list_set <- list()
dat <- rbind(dt_train,dt_test)
pred_col <- 'flag'
fit_list <- list()

if(T){
  factor_variable <- setdiff(cat_var,pred_col)
  cont_variable <- cont_var
  model_var1 <- c('GA_week','BW',"Central_nervous_system_abnormalities",
                  "Evidence_of_metabolic_disease","Respiratory_system_abnormalities",
                  "Infections_and_immune_abnormalities",
                  "Craniofacial_anomalies")
  model_var2 <- c(model_var1,"score")
  trainx1 <- dt_train[,model_var1] %>% as.matrix()
  trainx2 <- dt_train[,model_var2] %>% as.matrix()
  trainy <- dt_train[,pred_col] %>% as.factor()
  testx1 <- dt_test[,model_var1] %>% as.matrix()
  testx2 <- dt_test[,model_var2] %>% as.matrix()
  testy <- dt_test[,pred_col] %>% as.factor()
  valix1 <- dt_vali[,model_var1] %>% as.matrix()
  valix2 <- dt_vali[,model_var2] %>% as.matrix()
  valiy <- dt_vali[,pred_col] %>% as.factor()
  
  
  fitControl <- trainControl(method="repeatedcv", number=10, repeats=3, returnResamp="all")
  set.seed(42)
  fit_list[["gbm1"]] <- fit_gbm1 <- train(trainx1, trainy, method="gbm", trControl=fitControl, verbose=F, metric = 'Accuracy')
  fit_list[["gbm2"]] <- fit_gbm2 <- train(trainx2, trainy, method="gbm", trControl=fitControl, verbose=F, metric = 'Accuracy')
  
  Y_test <- as.numeric(as.character(testy))
  Y_vali <- as.numeric(as.character(valiy))
}


# test set
## AUC
auc_list <- list()
for (i in c("gbm1","gbm2")) {
  # i <- 'gbm1'
  # print(i)
  fit <- fit_list[[i]]
  if(str_detect(i,'1')){
    pred_Y1 <- predict(fit,newdata=trainx1,type = "prob")
    pred_Y_test1 <- predict(fit, newdata=testx1, type = 'prob')
  }else if(str_detect(i,'2')){
    pred_Y1 <- predict(fit,newdata=trainx2,type = "prob")
    pred_Y_test1 <- predict(fit, newdata=testx2, type = 'prob')
  }else if(str_detect(i,'3')){
    pred_Y1 <- predict(fit,newdata=trainx3,type = "prob")
    pred_Y_test1 <- predict(fit, newdata=testx3, type = 'prob')
  }
  
  ### ROC Curve
  Y_train <- as.numeric(as.character(trainy))
  pred_Y1 <- as.numeric(pred_Y1[,2])
  Y_test <- as.numeric(as.character(testy))
  pred_Y_test1 <- as.numeric(pred_Y_test1[,2])
  
  pdf(sprintf('./ROC_%s.pdf',i))
  rocobj <- plot.roc(Y_train,pred_Y1,print.auc=TRUE,#smooth=T,
                     auc.polygon=TRUE,print.thres='best',print.thres.best.method='youden',
                     ci=T,
                     max.auc.polygon=TRUE,main='Trainingset');
  ciobj <- ci.se(rocobj, # CI of sensitivity
                 specificities=seq(0, 1, 0.1)) # over a select set of specificities
  plot(ciobj, type="shape", col="#1c61b6AA") # plot as a blue shape
  plot(ci(rocobj, of="thresholds", thresholds="best")) # add one threshold
  
  auc_list[[i]] <- rocobj <- plot.roc(Y_test,pred_Y_test1,print.auc=TRUE,#smooth=T,
                                      auc.polygon=TRUE,print.thres='best',print.thres.best.method='youden',
                                      ci=T,
                                      max.auc.polygon=TRUE,main='Testingset');
  ciobj <- ci.se(rocobj, # CI of sensitivity
                 specificities=seq(0, 1, 0.1)) # over a select set of specificities
  plot(ciobj, type="shape", col="#1c61b6AA") # plot as a blue shape
  plot(ci(rocobj, of="thresholds", thresholds="best")) # add one threshold
  dev.off()
}
auc_list_set[[pred_col]] <- auc_list


auc_df_list <- list()
for (i in names(auc_list_set)) {
  # i <- names(auc_list_set)[1]
  se_tmp <- lapply(auc_list_set[[i]], function(x){
    se <- x[["sensitivities"]]
    sp <- x[["specificities"]]
    tmp <- se+sp
    max_index <- which(tmp==max(tmp))
    return(data.frame(Se=se[max_index],Sp=sp[max_index]))
  })
  se_sp_df <- do.call('rbind',se_tmp)
  auc_tmp <- lapply(auc_list_set[[i]], ci) %>% as.data.frame() %>% t() %>% as.data.frame()
  colnames(auc_tmp) <- c("lower","AUC","upper")
  auc_tmp$pred <- i
  auc_tmp$method <- rownames(auc_tmp)
  auc_se <- cbind(auc_tmp,se_sp_df)
  auc_df_list[[i]] <- auc_se
}
auc_df <- do.call('rbind',auc_df_list)
auc_df$`AUC(95%CI)` <- str_c(round(auc_df$AUC,2)," [",round(auc_df$lower,2),", ",round(auc_df$upper,2),"]")
auc_df$Se <- round(auc_df$Se,2)
auc_df$Sp <- round(auc_df$Sp,2)
auc_df <- auc_df %>% dplyr::select(pred,method,`AUC(95%CI)`,Se,Sp)
auc_df
rio::export(auc_df,file = "./AUC_model.xlsx")


gbm_fit1 <- fit_list[[1]]
gbm_fit2 <- fit_list[[2]]
pred_Y_test1 <- as.numeric(predict(gbm_fit1, newdata=testx1, type = 'prob')[,2])
pred_Y_test2 <- as.numeric(predict(gbm_fit2, newdata=testx2, type = 'prob')[,2])

# Combine two ROC curve
if(T){
  pdf("./ROC_combine.pdf",width = 5,height = 5)
  opar <- par(no.readonly = TRUE)
  par(pin=c(12,12),oma=c(1,1,1,1))
  rocobj1 <- plot.roc(Y_test,pred_Y_test1,
                      main="ROC curves",percent=TRUE, col="#1c61b6")
  rocobj2 <- lines.roc(Y_test,pred_Y_test2, percent=TRUE, col="#008600")
  testobj1 <- roc.test(rocobj1, rocobj2)
  text(60, 40, labels=paste("P =", format.pval(testobj1$p.value,digits = 3,scientific=T)), adj=c(0, .5),cex=0.8)
  legend("bottomright", 
         legend=c("Clinical factors","Clinical factors + Genetic factors"),
         col=c("#1c61b6", "#008600","#00BFFF","#FFC0CB"),bty="n",lty = c(1,1,1,2),lwd=2,cex=0.8,y.intersp=1,title="Area under the ROC curve")
  par(opar)
  dev.off()
}

# validation set 
## AUC
auc_list <- list()
for (i in c("gbm1","gbm2")) {
  fit <- fit_list[[i]]
  if(str_detect(i,'1')){
    pred_Y_vali1 <- predict(fit, newdata=valix1, type = 'prob')
  }else if(str_detect(i,'2')){
    pred_Y_vali1 <- predict(fit, newdata=valix2, type = 'prob')
  }
  
  ### ROC Curve
  Y_vali <- as.numeric(as.character(valiy))
  pred_Y_vali1 <- as.numeric(pred_Y_vali1[,2])
  
  
  auc_list[[i]] <- rocobj <- plot.roc(Y_vali,pred_Y_vali1,print.auc=TRUE,#smooth=T,
                                      auc.polygon=TRUE,print.thres='best',print.thres.best.method='youden',
                                      ci=T,
                                      max.auc.polygon=TRUE,main='Testingset');
}
auc_list_set[[pred_col]] <- auc_list


auc_df_list <- list()
for (i in names(auc_list_set)) {
  # i <- names(auc_list_set)[1]
  se_tmp <- lapply(auc_list_set[[i]], function(x){
    se <- x[["sensitivities"]]
    sp <- x[["specificities"]]
    tmp <- se+sp
    max_index <- which(tmp==max(tmp))
    return(data.frame(Se=se[max_index],Sp=sp[max_index]))
  })
  se_sp_df <- do.call('rbind',se_tmp)
  auc_tmp <- lapply(auc_list_set[[i]], ci) %>% as.data.frame() %>% t() %>% as.data.frame()
  colnames(auc_tmp) <- c("lower","AUC","upper")
  auc_tmp$pred <- i
  auc_tmp$method <- rownames(auc_tmp)
  auc_se <- cbind(auc_tmp,se_sp_df)
  auc_df_list[[i]] <- auc_se
}
auc_df <- do.call('rbind',auc_df_list)
auc_df$`AUC(95%CI)` <- str_c(round(auc_df$AUC,2)," [",round(auc_df$lower,2),", ",round(auc_df$upper,2),"]")
auc_df$Se <- round(auc_df$Se,2)
auc_df$Sp <- round(auc_df$Sp,2)
auc_df <- auc_df %>% dplyr::select(pred,method,`AUC(95%CI)`,Se,Sp)
auc_df
rio::export(auc_df,file = "./AUC_model_vali.xlsx")


gbm_fit1 <- fit_list[[1]]
gbm_fit2 <- fit_list[[2]]
pred_Y_vali1 <- as.numeric(predict(gbm_fit1, newdata=valix1, type = 'prob')[,2])
pred_Y_vali2 <- as.numeric(predict(gbm_fit2, newdata=valix2, type = 'prob')[,2])

# Combine two ROC curve
if(T){
  pdf("./ROC_combine_vali.pdf",width = 5,height = 5)
  opar <- par(no.readonly = TRUE)
  par(pin=c(12,12),oma=c(1,1,1,1))
  rocobj1 <- plot.roc(Y_vali,pred_Y_vali1,
                      main="ROC curves",percent=TRUE, col="#1c61b6")
  rocobj2 <- lines.roc(Y_vali,pred_Y_vali2, percent=TRUE, col="#008600")
  testobj1 <- roc.test(rocobj1, rocobj2)
  text(60, 40, labels=paste("P =", format.pval(testobj1$p.value,digits = 3,scientific=T)), adj=c(0, .5),cex=0.8)
  legend("bottomright", 
         legend=c("Clinical factors","Clinical factors + Genetic factors"),
         col=c("#1c61b6", "#008600","#00BFFF","#FFC0CB"),bty="n",lty = c(1,1,1,2),lwd=2,cex=0.8,y.intersp=1,title="Area under the ROC curve")
  par(opar)
  dev.off()
}

