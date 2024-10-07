#ROC curve for testing dataset
#set working directory
setwd("./path/to/Colon_data")
library(pROC)
library(ggplot2)

#GSE17538
d1<- read.csv('testing dataset_GSE17538.csv',header=T)
d1<- d1[which(d1$RFS_IND!=''),]
d1$event<- ifelse(d1$RFS>36,0,d1$RFS_IND)
event <- d1$event
XGBoost=d1$XG_riskscore
RandomForest=d1$RF_riskscore
LASSO=d1$LASSO_riskscore
Cox_PH=d1$Cox_riskscore
res <- roc(event ~ XGBoost +RandomForest + LASSO +Cox_PH,data=d1,
           aur=TRUE)
#pdf("ROCcurve_GSE17538.pdf")
p <-ggroc(res,legacy.axes = TRUE)+
  geom_segment(aes(x=0,xend=1,y=0,yend=1),
               color="darkgrey",linetype=5)+
  scale_color_discrete(labels=paste(c('XGBoost','RandomForest','LASSO','Cox_PH'), 
                                    rep(" AUC=", 4), 
                                    c(round(res$`XGBoost`$auc,3),round(res$`RandomForest`$auc,3),round(res$`LASSO`$auc,3),round(res$`Cox_PH`$auc,3)))) + 
  theme_bw()+ ggtitle('ROC_GSE17538') +
  theme(legend.title=element_blank(),
        legend.position = c(0.75, 0.25),
        legend.background = element_rect(fill = "white", color = "white"))
p

#GSE39582
d1<- read.csv('testing dataset_GSE39582.csv',header=T)
d1<- d1[which(d1$RFS_IND!=''),]
d1$event<- ifelse(d1$RFS>36,0,d1$RFS_IND)
event <- d1$event
XGBoost=d1$XG_riskscore
RandomForest=d1$RF_riskscore
LASSO=d1$LASSO_riskscore
Cox_PH=d1$Cox_riskscore

res <- roc(event ~ XGBoost +RandomForest + LASSO +Cox_PH,data=d1,
           aur=TRUE)
#pdf("ROCcurve_GSE39582.pdf")
p1 <-ggroc(res,legacy.axes = TRUE)+
  geom_segment(aes(x=0,xend=1,y=0,yend=1),
               color="darkgrey",linetype=5)+
  scale_color_discrete(labels=paste(c('XGBoost','RandomForest','LASSO','Cox_PH'), 
                                    rep(" AUC=", 4), 
                                    c(round(res$`XGBoost`$auc,3),round(res$`RandomForest`$auc,3),round(res$`LASSO`$auc,3),round(res$`Cox_PH`$auc,3)))) + 
  theme_bw()+ ggtitle('ROC_GSE39582') + 
  theme(legend.title=element_blank(),
        legend.position = c(0.75, 0.25),
        legend.background = element_rect(fill = "white", color = "white"))

p1


#GSE33113
d1<- read.csv('testing dataset_GSE33113.csv',header=T)
d1<- d1[which(d1$RFS_IND!=''),]
d1$event<- ifelse(d1$RFS>36,0,d1$RFS_IND)
event <- d1$event
XGBoost=d1$XG_riskscore
RandomForest=d1$RF_riskscore
LASSO=d1$LASSO_riskscore
Cox_PH=d1$Cox_riskscore

res <- roc(event ~ XGBoost +RandomForest + LASSO +Cox_PH,data=d1,
           aur=TRUE)
#pdf("ROCcurve_36M_GSE33113_new.pdf")
p2 <-ggroc(res,legacy.axes = TRUE)+
  geom_segment(aes(x=0,xend=1,y=0,yend=1),
               color="darkgrey",linetype=5)+
  scale_color_discrete(labels=paste(c('XGBoost','RandomForest','LASSO','Cox_PH'), 
                                    rep(" AUC=", 4), 
                                    c(round(res$`XGBoost`$auc,3),round(res$`RandomForest`$auc,3),round(res$`LASSO`$auc,3),round(res$`Cox_PH`$auc,3)))) + 
  theme_bw()+ ggtitle('ROC_GSE33113') + 
  theme(legend.title=element_blank(),
        legend.position = c(0.75, 0.25),
        legend.background = element_rect(fill = "white", color = "white"))

  p2