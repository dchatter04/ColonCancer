rm(list=ls())
#ROC curve for testing dataset
#set working directory
#setwd("./path/to/Colon_data")
setwd("/Users/dchatter/Documents/Frontiers_submission_Xiuxiu/Colon_data/Colon_data")

library(pROC)
library(ggplot2)
library(MLmetrics)

#GSE17538 ----
d1<- read.csv('testing dataset_GSE17538.csv',header=T)
d2 <- read.csv('8genes_GSE17538_dataALL.csv',header=T)
d1$eightgene_riskscore <- d2$score
d1<- d1[which(d1$RFS_IND!=''),]
d1$event<- ifelse(d1$RFS>36,0,d1$RFS_IND)
event <- d1$event
XGBoost=d1$XG_riskscore
RandomForest=d1$RF_riskscore
LASSO=d1$LASSO_riskscore
Cox_PH=d1$Cox_riskscore
eightgene = d1$eightgene_riskscore
res <- roc(event ~ XGBoost +RandomForest + LASSO +Cox_PH+eightgene_riskscore,data=d1,
           aur=TRUE)
p <-ggroc(res,legacy.axes = TRUE)+
  geom_segment(aes(x=0,xend=1,y=0,yend=1),
               color="darkgrey",linetype=5)+
  scale_color_discrete(labels=paste(c('XGBoost','RandomForest','LASSO','Cox_PH','8-gene'), 
                                    rep(" AUC=", 5), 
                                    c(round(res$`XGBoost`$auc,3),round(res$`RandomForest`$auc,3),round(res$`LASSO`$auc,3),round(res$`Cox_PH`$auc,3),round(res$`eightgene_riskscore`$auc,3)))) + 
  theme_bw()+ ggtitle('ROC_GSE17538') +
  theme(legend.title=element_blank(),
        legend.position = c(0.75, 0.25),
        legend.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())  # Remove minor grid lines
p
ggsave(p,file="/Users/dchatter/Documents/Frontiers_submission_Xiuxiu/Travis_codes/Travis_codes_corrected/ROCcurve_GSE17538_Nov4.pdf",height=6,width=6)

# xbg----
a1 <- roc(event,XGBoost)
auc<- auc(a1)
prauc <- PRAUC(XGBoost,event)
sens <- Sensitivity(event,ifelse(XGBoost>median(XGBoost),1,0),positive = "1")
spec <- Specificity(event,ifelse(XGBoost>median(XGBoost),1,0),positive = "1")
f1 <- F1_Score(event,ifelse(XGBoost>median(XGBoost),1,0),positive="1")
# Print results vertically
cat("AUC=", round(auc, 3), "\n",
    "PRAUC=", round(prauc, 3), "\n",
    "Sensitivity=", round(sens, 3), "\n",
    "Specificity=", round(spec, 3), "\n",
    "F1=", round(f1, 3), "\n")

#rf----
a1 <- roc(event,RandomForest)
auc<- auc(a1)
prauc <- PRAUC(RandomForest,event)
sens <- Sensitivity(event,ifelse(RandomForest>median(RandomForest),1,0),positive = "1")
spec <- Specificity(event,ifelse(RandomForest>median(RandomForest),1,0),positive = "1")
f1 <- F1_Score(event,ifelse(RandomForest>median(RandomForest),1,0),positive="1")
# Print results vertically
cat("AUC=", round(auc, 3), "\n",
    "PRAUC=", round(prauc, 3), "\n",
    "Sensitivity=", round(sens, 3), "\n",
    "Specificity=", round(spec, 3), "\n",
    "F1=", round(f1, 3), "\n")
#lasso----
a1 <- roc(event,LASSO)
auc<- auc(a1)
prauc <- PRAUC(LASSO,event)
sens <- Sensitivity(event,ifelse(LASSO>median(LASSO),1,0),positive = "1")
spec <- Specificity(event,ifelse(LASSO>median(LASSO),1,0),positive = "1")
f1 <- F1_Score(event,ifelse(LASSO>median(LASSO),1,0),positive="1")
# Print results vertically
cat("AUC=", round(auc, 3), "\n",
    "PRAUC=", round(prauc, 3), "\n",
    "Sensitivity=", round(sens, 3), "\n",
    "Specificity=", round(spec, 3), "\n",
    "F1=", round(f1, 3), "\n")
#coxph----
a1 <- roc(event,Cox_PH)
auc<- auc(a1)
prauc <- PRAUC(Cox_PH,event)
sens <- Sensitivity(event,ifelse(Cox_PH>median(Cox_PH),1,0),positive = "1")
spec <- Specificity(event,ifelse(Cox_PH>median(Cox_PH),1,0),positive = "1")
f1 <- F1_Score(event,ifelse(Cox_PH>median(Cox_PH),1,0),positive="1")
# Print results vertically
cat("AUC=", round(auc, 3), "\n",
    "PRAUC=", round(prauc, 3), "\n",
    "Sensitivity=", round(sens, 3), "\n",
    "Specificity=", round(spec, 3), "\n",
    "F1=", round(f1, 3), "\n")

#8gene----
a1 <- roc(event,eightgene)
auc<- auc(a1)
prauc <- PRAUC(eightgene,event)
sens <- Sensitivity(event,ifelse(eightgene>median(eightgene),1,0),positive = "1")
spec <- Specificity(event,ifelse(eightgene>median(eightgene),1,0),positive = "1")
f1 <- F1_Score(event,ifelse(eightgene>median(eightgene),1,0),positive="1")
# Print results vertically
cat("AUC=", round(auc, 3), "\n",
    "PRAUC=", round(prauc, 3), "\n",
    "Sensitivity=", round(sens, 3), "\n",
    "Specificity=", round(spec, 3), "\n",
    "F1=", round(f1, 3), "\n")
######################################################
#GSE39582 ----
d1<- read.csv('testing dataset_GSE39582.csv',header=T)
d2 <- read.csv('8genes_GSE39582_dataALL.csv',header=T)
d1$eightgene_riskscore <- d2$score
d1<- d1[which(d1$RFS_IND!=''),]
d1$event<- ifelse(d1$RFS>36,0,d1$RFS_IND)
event <- d1$event
XGBoost=d1$XG_riskscore
RandomForest=d1$RF_riskscore
LASSO=d1$LASSO_riskscore
Cox_PH=d1$Cox_riskscore
eightgene = d1$eightgene_riskscore
res <- roc(event ~ XGBoost +RandomForest + LASSO +Cox_PH+eightgene_riskscore,data=d1,
           aur=TRUE)
p1 <-ggroc(res,legacy.axes = TRUE)+
  geom_segment(aes(x=0,xend=1,y=0,yend=1),
               color="darkgrey",linetype=5)+
  scale_color_discrete(labels=paste(c('XGBoost','RandomForest','LASSO','Cox_PH','8-gene'), 
                                    rep(" AUC=", 5), 
                                    c(round(res$`XGBoost`$auc,3),round(res$`RandomForest`$auc,3),round(res$`LASSO`$auc,3),round(res$`Cox_PH`$auc,3),round(res$`eightgene_riskscore`$auc,3)))) + 
  theme_bw()+ ggtitle('ROC_GSE39582') +
  theme(legend.title=element_blank(),
        legend.position = c(0.75, 0.25),
        legend.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())  # Remove minor grid lines
p1
ggsave(p1,file="/Users/dchatter/Documents/Frontiers_submission_Xiuxiu/Travis_codes/Travis_codes_corrected/ROCcurve_GSE39582_Nov4.pdf",height=6,width=6)
#AUC(XGBoost,event)
#AUC(RandomForest,event)
#AUC(LASSO,event)
#AUC(Cox_PH,event)
#AUC(eightgene,event)
# xbg----
a1 <- roc(event,XGBoost)
auc<- auc(a1)
prauc <- PRAUC(XGBoost,event)
sens <- Sensitivity(event,ifelse(XGBoost>median(XGBoost),1,0),positive = "1")
spec <- Specificity(event,ifelse(XGBoost>median(XGBoost),1,0),positive = "1")
f1 <- F1_Score(event,ifelse(XGBoost>median(XGBoost),1,0),positive="1")
# Print results vertically
cat("AUC=", round(auc, 3), "\n",
    "PRAUC=", round(prauc, 3), "\n",
    "Sensitivity=", round(sens, 3), "\n",
    "Specificity=", round(spec, 3), "\n",
    "F1=", round(f1, 3), "\n")

#rf----
a1 <- roc(event,RandomForest)
auc<- auc(a1)
prauc <- PRAUC(RandomForest,event)
sens <- Sensitivity(event,ifelse(RandomForest>median(RandomForest),1,0),positive = "1")
spec <- Specificity(event,ifelse(RandomForest>median(RandomForest),1,0),positive = "1")
f1 <- F1_Score(event,ifelse(RandomForest>median(RandomForest),1,0),positive="1")
# Print results vertically
cat("AUC=", round(auc, 3), "\n",
    "PRAUC=", round(prauc, 3), "\n",
    "Sensitivity=", round(sens, 3), "\n",
    "Specificity=", round(spec, 3), "\n",
    "F1=", round(f1, 3), "\n")
#lasso----
a1 <- roc(event,LASSO)
auc<- auc(a1)
prauc <- PRAUC(LASSO,event)
sens <- Sensitivity(event,ifelse(LASSO>median(LASSO),1,0),positive = "1")
spec <- Specificity(event,ifelse(LASSO>median(LASSO),1,0),positive = "1")
f1 <- F1_Score(event,ifelse(LASSO>median(LASSO),1,0),positive="1")
# Print results vertically
cat("AUC=", round(auc, 3), "\n",
    "PRAUC=", round(prauc, 3), "\n",
    "Sensitivity=", round(sens, 3), "\n",
    "Specificity=", round(spec, 3), "\n",
    "F1=", round(f1, 3), "\n")
#coxph----
a1 <- roc(event,Cox_PH)
auc<- auc(a1)
prauc <- PRAUC(Cox_PH,event)
sens <- Sensitivity(event,ifelse(Cox_PH>median(Cox_PH),1,0),positive = "1")
spec <- Specificity(event,ifelse(Cox_PH>median(Cox_PH),1,0),positive = "1")
f1 <- F1_Score(event,ifelse(Cox_PH>median(Cox_PH),1,0),positive="1")
# Print results vertically
cat("AUC=", round(auc, 3), "\n",
    "PRAUC=", round(prauc, 3), "\n",
    "Sensitivity=", round(sens, 3), "\n",
    "Specificity=", round(spec, 3), "\n",
    "F1=", round(f1, 3), "\n")

#8gene----
a1 <- roc(event,eightgene)
auc<- auc(a1)
prauc <- PRAUC(eightgene,event)
sens <- Sensitivity(event,ifelse(eightgene>median(eightgene),1,0),positive = "1")
spec <- Specificity(event,ifelse(eightgene>median(eightgene),1,0),positive = "1")
f1 <- F1_Score(event,ifelse(eightgene>median(eightgene),1,0),positive="1")
# Print results vertically
cat("AUC=", round(auc, 3), "\n",
    "PRAUC=", round(prauc, 3), "\n",
    "Sensitivity=", round(sens, 3), "\n",
    "Specificity=", round(spec, 3), "\n",
    "F1=", round(f1, 3), "\n")


###############################################
#GSE33113 ----
d1<- read.csv('testing dataset_GSE33113.csv',header=T)
d2 <- read.csv('8genes_GSE33113_dataALL.csv',header=T)
d1$eightgene_riskscore <- d2$score
d1<- d1[which(d1$RFS_IND!=''),]
d1$event<- ifelse(d1$RFS>36,0,d1$RFS_IND)
event <- d1$event
XGBoost=d1$XG_riskscore
RandomForest=d1$RF_riskscore
LASSO=d1$LASSO_riskscore
Cox_PH=d1$Cox_riskscore
eightgene = d1$eightgene_riskscore
res <- roc(event ~ XGBoost +RandomForest + LASSO +Cox_PH+eightgene_riskscore,data=d1,
           aur=TRUE)
p2 <-ggroc(res,legacy.axes = TRUE)+
  geom_segment(aes(x=0,xend=1,y=0,yend=1),
               color="darkgrey",linetype=5)+
  scale_color_discrete(labels=paste(c('XGBoost','RandomForest','LASSO','Cox_PH','8-gene'), 
                                    rep(" AUC=", 5), 
                                    c(round(res$`XGBoost`$auc,3),round(res$`RandomForest`$auc,3),round(res$`LASSO`$auc,3),round(res$`Cox_PH`$auc,3),round(res$`eightgene_riskscore`$auc,3)))) + 
  theme_bw()+ ggtitle('ROC_GSE33113') +
  theme(legend.title=element_blank(),
        legend.position = c(0.75, 0.25),
        legend.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())  # Remove minor grid lines
p2
ggsave(p2,file="/Users/dchatter/Documents/Frontiers_submission_Xiuxiu/Travis_codes/Travis_codes_corrected/ROCcurve_GSE33113_Nov4.pdf",height=6,width=6)
#AUC(XGBoost,event)
#AUC(RandomForest,event)
#AUC(LASSO,event)
#AUC(Cox_PH,event)
#AUC(eightgene,event)
# xbg----
a1 <- roc(event,XGBoost)
auc<- auc(a1)
prauc <- PRAUC(XGBoost,event)
sens <- Sensitivity(event,ifelse(XGBoost>median(XGBoost),1,0),positive = "1")
spec <- Specificity(event,ifelse(XGBoost>median(XGBoost),1,0),positive = "1")
f1 <- F1_Score(event,ifelse(XGBoost>median(XGBoost),1,0),positive="1")
# Print results vertically
cat("AUC=", round(auc, 3), "\n",
    "PRAUC=", round(prauc, 3), "\n",
    "Sensitivity=", round(sens, 3), "\n",
    "Specificity=", round(spec, 3), "\n",
    "F1=", round(f1, 3), "\n")

#rf----
a1 <- roc(event,RandomForest)
auc<- auc(a1)
prauc <- PRAUC(RandomForest,event)
sens <- Sensitivity(event,ifelse(RandomForest>median(RandomForest),1,0),positive = "1")
spec <- Specificity(event,ifelse(RandomForest>median(RandomForest),1,0),positive = "1")
f1 <- F1_Score(event,ifelse(RandomForest>median(RandomForest),1,0),positive="1")
# Print results vertically
cat("AUC=", round(auc, 3), "\n",
    "PRAUC=", round(prauc, 3), "\n",
    "Sensitivity=", round(sens, 3), "\n",
    "Specificity=", round(spec, 3), "\n",
    "F1=", round(f1, 3), "\n")
#lasso----
a1 <- roc(event,LASSO)
auc<- auc(a1)
prauc <- PRAUC(LASSO,event)
sens <- Sensitivity(event,ifelse(LASSO>median(LASSO),1,0),positive = "1")
spec <- Specificity(event,ifelse(LASSO>median(LASSO),1,0),positive = "1")
f1 <- F1_Score(event,ifelse(LASSO>median(LASSO),1,0),positive="1")
# Print results vertically
cat("AUC=", round(auc, 3), "\n",
    "PRAUC=", round(prauc, 3), "\n",
    "Sensitivity=", round(sens, 3), "\n",
    "Specificity=", round(spec, 3), "\n",
    "F1=", round(f1, 3), "\n")
#coxph----
a1 <- roc(event,Cox_PH)
auc<- auc(a1)
prauc <- PRAUC(Cox_PH,event)
sens <- Sensitivity(event,ifelse(Cox_PH>median(Cox_PH),1,0),positive = "1")
spec <- Specificity(event,ifelse(Cox_PH>median(Cox_PH),1,0),positive = "1")
f1 <- F1_Score(event,ifelse(Cox_PH>median(Cox_PH),1,0),positive="1")
# Print results vertically
cat("AUC=", round(auc, 3), "\n",
    "PRAUC=", round(prauc, 3), "\n",
    "Sensitivity=", round(sens, 3), "\n",
    "Specificity=", round(spec, 3), "\n",
    "F1=", round(f1, 3), "\n")

#8gene----
a1 <- roc(event,eightgene)
auc<- auc(a1)
prauc <- PRAUC(eightgene,event)
sens <- Sensitivity(event,ifelse(eightgene>median(eightgene),1,0),positive = "1")
spec <- Specificity(event,ifelse(eightgene>median(eightgene),1,0),positive = "1")
f1 <- F1_Score(event,ifelse(eightgene>median(eightgene),1,0),positive="1")
# Print results vertically
cat("AUC=", round(auc, 3), "\n",
    "PRAUC=", round(prauc, 3), "\n",
    "Sensitivity=", round(sens, 3), "\n",
    "Specificity=", round(spec, 3), "\n",
    "F1=", round(f1, 3), "\n")
