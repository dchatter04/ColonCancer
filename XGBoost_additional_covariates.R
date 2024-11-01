#set working directory
rm(list=ls())
setwd("~/Colon_data")
library('survival')
library("xgboost")
library('gbm')
library("caret")
library("datasets")
library("Hmisc")
library('ggplot2')
library(survcomp)
library(MLmetrics)

df <-read.csv("TCGA_COAD_clinic_clean.csv",header=T,row.names=1)
df1 <- df[which(df$DFS !="NA" & (df$DFS_IND !="NA")),] 
df1$T_sign <-ifelse(df1$DFS_IND==1,df1$DFS,-df1$DFS)
train=df1
rownames.train <- rownames(train)
train1<- data.frame(mapply(train, FUN=as.numeric))
rownames(train1)<- rownames.train
dtrain <- xgb.DMatrix(as.matrix(train1[,1:75]), label=train1$T_sign)

#including the age and disease atage information in the XGBoost model.
#following are the set of genes identified by the main XGBoost model
top = c("UGT2A3","SLC27A5","ATP1A1","GNAS","ATP1A4","BAAT","ABCC2","ATP1B1","KCNN2","ADCY5","SLC2A1","ABCC4","SCTR","HMGCR","AQP9")
train2 = train[,c(top,"Diagnosis.Age")]
train2$Sex = ifelse(train$Sex=="Male",1,0)
train2$Stage = train$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code
train2$Stage[train2$Stage %in% c("Stage I", "Stage IA","Stage IB")] = "1"
train2$Stage[train2$Stage %in% c("Stage II","Stage IIA","Stage IIB","Stage IIC")] = "2"
train2$Stage[train2$Stage %in% c("Stage III","Stage IIIA","Stage IIIB","Stage IIIC")] = "3"
train2$Stage[train2$Stage %in% c("Stage IV","Stage IVA","Stage IVB")] = "4"
train2 =na.omit(train2)
colnames(train2) = c(top,"Age","Gender","Stage")
train2 = train2[,c(top,"Age","Stage")]
train2 <- data.frame(mapply(train2, FUN=as.numeric))
train2$Age = log2(train2$Age)
label = train1$T_sign[!is.na(train$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code)]
dtrainAll <- xgb.DMatrix(as.matrix(train2),label=label)

set.seed(123)
rnds = 60
param <- list(max_depth=6,verbosity=0,objective = "survival:cox",
              eval_metric = "cox-nloglik",subsample=0.5,gamma=0)
modelAll <- xgb.train(param, dtrainAll, nrounds=rnds)
#saveRDS(modelAll,file="modelAll_final.rds")
#modelAll = readRDS("modelAll.rds")

# GSE17538
gse_id<- "GSE17538"
used_gene <-top
dataGene<- read.csv(paste(gse_id,"_exprs.txt", sep=""), header=T,
                    sep='\t', row.names=1,stringsAsFactors=FALSE)

names(dataGene)<- gsub('.CEL.gz', '', names(dataGene))
names(dataGene)<- as.character(data.frame(do.call('rbind', 
                                                  strsplit(as.character(names(dataGene)),'_',fixed=TRUE)))[,1])
dataGene<- dataGene[used_gene,]
dataClinic<- read.csv(paste(gse_id,"_clinic.csv", sep=""), header=T,
                      row.names=1,stringsAsFactors=FALSE)
dataAll<- merge(data.frame(t(dataGene)),dataClinic, by='row.names', all.x=T)
row.names(dataAll) <-dataAll[,1]
dataAll <-dataAll[,-1]
test_top <- dataAll[,(names(dataAll) %in% c(top,"Age","Gender","Stage","RFS","RFS_IND","DFS","DFS_IND"))]
test_top$Gender = ifelse(test_top$Gender == "MALE",1,0)
test_top2 = test_top[,c(top,"Age","Stage")]
test_top2$Age = log2(test_top2$Age)
test_top2 = data.frame(mapply(test_top2, FUN=as.numeric))
test_top2_xgb = xgb.DMatrix(as.matrix(test_top2))
test_top2_predsAll = predict(modelAll,test_top2_xgb)
ci =concordance.index(test_top2_predsAll[!is.na(test_top$RFS_IND)], as.numeric(test_top$RFS)[!is.na(test_top$RFS_IND)], as.numeric(test_top$RFS_IND)[!is.na(test_top$RFS_IND)])
ci$c.index


event = rep(NA,dim(test_top)[1])
event[test_top$RFS>36] = 0
event[test_top$RFS<=36 & test_top$RFS_IND==1] = 1
AUC(test_top2_predsAll[!is.na(event)],event[!is.na(event)])
PRAUC(test_top2_predsAll[!is.na(event)],event[!is.na(event)])
Sensitivity(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0))
Specificity(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0))
F1_Score(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0),positive=1)
roc36_GSE17538 <- roc(event[!is.na(event)], test_top2_predsAll[!is.na(event)])


event = rep(NA,dim(test_top)[1])
event[test_top$RFS>60] = 0
event[test_top$RFS<=60 & test_top$RFS_IND==1] = 1
AUC(test_top2_predsAll[!is.na(event)],event[!is.na(event)])
PRAUC(test_top2_predsAll[!is.na(event)],event[!is.na(event)])
Sensitivity(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0))
Specificity(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0))
F1_Score(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0),positive=1)
roc60_GSE17538 <- roc(event[!is.na(event)], test_top2_predsAll[!is.na(event)])

event = rep(NA,dim(test_top)[1])
event[test_top$RFS>120] = 0
event[test_top$RFS<=120 & test_top$RFS_IND==1] = 1
AUC(test_top2_predsAll[!is.na(event)],event[!is.na(event)])
PRAUC(test_top2_predsAll[!is.na(event)],event[!is.na(event)])
Sensitivity(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0))
Specificity(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0))
F1_Score(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0),positive=1)

outGSE17538 = data.frame(patients=rownames(test_top),hazards=test_top2_predsAll)
#write.csv(outGSE17538,file="GSE17538_preds_final.csv")

# GSE33113
gse_id<- "GSE33113"
used_gene <-top
dataGene<- read.csv(paste(gse_id,"_exprs.txt", sep=""), header=T,
                    sep='\t', row.names=1,stringsAsFactors=FALSE)

names(dataGene)<- gsub('.CEL.gz', '', names(dataGene))
names(dataGene)<- as.character(data.frame(do.call('rbind', 
                                                  strsplit(as.character(names(dataGene)),'_',fixed=TRUE)))[,1])
dataGene<- dataGene[used_gene,]
dataClinic<- read.csv(paste(gse_id,"_clinic.csv", sep=""), header=T,
                      row.names=1,stringsAsFactors=FALSE)
dataAll<- merge(data.frame(t(dataGene)),dataClinic, by='row.names', all.x=T)
row.names(dataAll) <-dataAll[,1]
dataAll <-dataAll[,-1]
test_top <- dataAll[,(names(dataAll) %in% c(top,"age","sex","RFS","RFS_IND","DFS","DFS_IND"))]
test_top$RFS = round(test_top$RFS/30,2)
test_top$RFS_IND=1
test_top2 = test_top[,c(top,"age","sex")]
test_top2$sex = ifelse(test_top2$sex=="m","MALE","FEMALE")
colnames(test_top2) = c(top,"Age","Gender")
test_top2$Gender = ifelse(test_top2$Gender == "MALE",1,0)
test_top2 = test_top2[,c(top,"Age")]
#test_top2 = test_top2[,c(top,"Age","Gender")]
test_top2$Stage = 2
test_top2 = data.frame(mapply(test_top2, FUN=as.numeric))
test_top2$Age = log2(test_top2$Age)
test_top2_xgb = xgb.DMatrix(as.matrix(test_top2))
test_top2_predsAll = predict(modelAll,test_top2_xgb)
ci =concordance.index(test_top2_predsAll[!is.na(test_top$RFS_IND)], as.numeric(test_top$RFS), as.numeric(test_top$RFS_IND))
ci$c.index
event = rep(NA,dim(test_top)[1])
event[test_top$RFS>36] = 0
event[test_top$RFS<=36 & test_top$RFS_IND==1] = 1
AUC(test_top2_predsAll[!is.na(event)],event[!is.na(event)])
PRAUC(test_top2_predsAll[!is.na(event)],event[!is.na(event)])
Sensitivity(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0))
Specificity(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0))
F1_Score(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0),positive=1)
roc36_GSE33113 <- roc(event[!is.na(event)], test_top2_predsAll[!is.na(event)])

event = rep(NA,dim(test_top)[1])
event[test_top$RFS>60] = 0
event[test_top$RFS<=60 & test_top$RFS_IND==1] = 1
AUC(test_top2_predsAll[!is.na(event)],event[!is.na(event)])
PRAUC(test_top2_predsAll[!is.na(event)],event[!is.na(event)])
Sensitivity(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0))
Specificity(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0))
F1_Score(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0),positive=1)
roc60_GSE33113 <- roc(event[!is.na(event)], test_top2_predsAll[!is.na(event)])

outGSE33113 = data.frame(patients=rownames(test_top),hazards=test_top2_predsAll)
#write.csv(outGSE33113,file="GSE33113_preds_final.csv")


# GSE39582
gse_id<- "GSE39582"
used_gene <-top
dataGene<- read.csv(paste(gse_id,"_exprs.txt", sep=""), header=T,
                    sep='\t', row.names=1,stringsAsFactors=FALSE)

names(dataGene)<- gsub('.CEL.gz', '', names(dataGene))
names(dataGene)<- as.character(data.frame(do.call('rbind', 
                                                  strsplit(as.character(names(dataGene)),'_',fixed=TRUE)))[,1])
dataGene<- dataGene[used_gene,]
dataClinic<- read.csv(paste(gse_id,"_clinic.csv", sep=""), header=T,
                      row.names=1,stringsAsFactors=FALSE)
dataAll<- merge(data.frame(t(dataGene)),dataClinic, by='row.names', all.x=T)
row.names(dataAll) <-dataAll[,1]
dataAll <-dataAll[,-1]
test_top <- dataAll[,(names(dataAll) %in% c(top,"Age","Gender","Stage","RFS","RFS_IND","DFS","DFS_IND"))]
test_top2 = test_top[,c(top,"Age","Gender","Stage")]
test_top2$Gender = ifelse(test_top$Gender == "MALE",1,0)
test_top2$Age = log2(test_top2$Age)
test_top2 = test_top2[,c(top,"Age","Stage")]
test_top2 = data.frame(mapply(test_top2, FUN=as.numeric))
test_top2_xgb = xgb.DMatrix(as.matrix(test_top2))
test_top2_predsAll = predict(modelAll,test_top2_xgb)
ci =concordance.index(test_top2_predsAll[!is.na(test_top$RFS_IND)], as.numeric(test_top$RFS)[!is.na(test_top$RFS_IND)], as.numeric(test_top$RFS_IND)[!is.na(test_top$RFS_IND)])
ci$c.index
event = rep(NA,dim(test_top)[1])
event[test_top$RFS>36] = 0
event[test_top$RFS<=36 & test_top$RFS_IND==1] = 1
AUC(test_top2_predsAll[!is.na(event)],event[!is.na(event)])
PRAUC(test_top2_predsAll[!is.na(event)],event[!is.na(event)])
Sensitivity(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0))
Specificity(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0))
F1_Score(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0),positive=1)
roc36_GSE39582 <- roc(event[!is.na(event)], test_top2_predsAll[!is.na(event)])

event = rep(NA,dim(test_top)[1])
event[test_top$RFS>60] = 0
event[test_top$RFS<=60 & test_top$RFS_IND==1] = 1
AUC(test_top2_predsAll[!is.na(event)],event[!is.na(event)])
PRAUC(test_top2_predsAll[!is.na(event)],event[!is.na(event)])
Sensitivity(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0))
Specificity(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0))
F1_Score(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0),positive=1)
roc60_GSE39582 <- roc(event[!is.na(event)], test_top2_predsAll[!is.na(event)])

event = rep(NA,dim(test_top)[1])
event[test_top$RFS>120] = 0
event[test_top$RFS<=120 & test_top$RFS_IND==1] = 1
AUC(test_top2_predsAll[!is.na(event)],event[!is.na(event)])
PRAUC(test_top2_predsAll[!is.na(event)],event[!is.na(event)])
Sensitivity(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0))
Specificity(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0))
F1_Score(event[!is.na(event)],ifelse(test_top2_predsAll[!is.na(event)]>median(test_top2_predsAll[!is.na(event)]),1,0),positive=1)

outGSE39582 = data.frame(patients=rownames(test_top),hazards=test_top2_predsAll)
#write.csv(outGSE39582,file="GSE39582_preds_final.csv")


combined_labels <- c(paste0("GSE17538 (3yrs) [AUC=", round(roc36_GSE17538$auc, 3), "]"), 
                     paste0("GSE33113 (3yrs) [AUC=", round(roc36_GSE33113$auc, 3), "]"),
                     paste0("GSE39582 (3yrs) [AUC=", round(roc36_GSE39582$auc, 3), "]"),
                     paste0("GSE17538 (5yrs) [AUC=", round(roc60_GSE17538$auc, 3), "]"), 
                     paste0("GSE33113 (5yrs) [AUC=", round(roc60_GSE33113$auc, 3), "]"),
                     paste0("GSE39582 (5yrs) [AUC=", round(roc60_GSE39582$auc, 3), "]"))
colors <- c("GSE17538 (3yrs) [AUC=0.734]" = "#00AFBB", 
            "GSE33113 (3yrs) [AUC=0.586]" = "salmon", 
            "GSE39582 (3yrs) [AUC=0.631]" = "#52854C", 
            "GSE17538 (5yrs) [AUC=0.747]" = "#00AFBB", 
            "GSE33113 (5yrs) [AUC=0.72]" = "salmon", 
            "GSE39582 (5yrs) [AUC=0.625]" = "#52854C")
line_types <- c("GSE17538 (3yrs) [AUC=0.734]" = "solid", 
                "GSE33113 (3yrs) [AUC=0.586]" = "solid", 
                "GSE39582 (3yrs) [AUC=0.631]" = "solid", 
                "GSE17538 (5yrs) [AUC=0.747]" = "dashed", 
                "GSE33113 (5yrs) [AUC=0.72]" = "dashed", 
                "GSE39582 (5yrs) [AUC=0.625]" = "dashed")


plot.df <- list("GSE17538 (3yrs)" = roc36_GSE17538, 
                "GSE33113 (3yrs)" = roc36_GSE33113, 
                "GSE39582 (3yrs)"= roc36_GSE39582, 
                "GSE17538 (5yrs)" = roc60_GSE17538, 
                "GSE33113 (5yrs)" = roc60_GSE33113, 
                "GSE39582 (5yrs)"= roc60_GSE39582)
names(plot.df) <- combined_labels

p <- ggroc(plot.df, 
           aes = c("color", "linetype"), legacy.axes = TRUE) +
  scale_linetype_manual(values = line_types) + 
  scale_color_manual(values = colors) + 
  geom_segment(aes(x=0,xend=1,y=0,yend=1),
               color="darkgrey",linetype=5) +
  theme_minimal() +
  ggtitle('ROC for XBGoost including additional covariates') +
  theme(legend.title=element_blank(),
        legend.position = c(0.75, 0.25),
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add a border around the plot panel
        #legend.background = element_rect(color = "black", fill = NA),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())  # Remove minor grid lines
#p
ggsave(p, file="ROC_XGboost_stages.pdf",height=6,width=6)

