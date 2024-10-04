#set working directory
setwd("./path/to/Colon_data")
library('survival')
library("xgboost")
library('gbm')
library("caret")
library("datasets")
library("Hmisc")
library('ggplot2')

df <-read.csv("TCGA_COAD_clinic_clean.csv",header=T,row.names=1)
df1 <- df[which(df$DFS !="NA" & (df$DFS_IND !="NA")),] 
df1$T_sign <-ifelse(df1$DFS_IND==1,df1$DFS,-df1$DFS)
#randomize the clinical data to set up the training and test datasets
random.index <-sample(c(1:dim(df1)[1]),size=round(dim(df1)[1]))
train=df1[random.index,]
rownames.train <- rownames(train)
train1<- data.frame(mapply(train, FUN=as.numeric))
rownames(train1)<- rownames.train
dtrain <- xgb.DMatrix(as.matrix(train1[,1:75]), label=train1$T_sign)

#set a random seed and a set of parameters
set.seed(123)
param <- list(max_depth=6,verbosity=0,objective = "survival:cox",
              eval_metric = "cox-nloglik",subsample=0.5,gamma=0)
model <- xgb.train(param, dtrain, nrounds=60)

pred.train <-log(predict(model,dtrain))
Hmisc::rcorr.cens(-pred.train, Surv(train1$DFS, train1$DFS_IND))
time.interest <-sort(unique(train1$DFS[train1$DFS_IND==1]))
basehaz.cum <-basehaz.gbm(train1$DFS,train1$DFS_IND,pred.train,cumulative=TRUE)
surv.rate <-exp(-exp(pred.train)*basehaz.cum[10])
importance_matrix=xgb.importance(colnames(dtrain),model=model)

top_result <- data.frame(matrix(NA,nrow=0,ncol=16))
names(top_result) <- c("Top feature","train_pvalue","train_hr","train_concor",rep(c("test_dataset","test_pvalue","test_hr","test_concor"),3))
top_XGBoost <-data.frame(importance_matrix[1:15,])
top <-top_XGBoost$Feature
top_clinic <-cbind(top,colnames(df1[,c(76:83)]))
top_expr_clinic <- data.frame(df1[,which((colnames(df1) %in% top_clinic))])
top_clinic <- top_expr_clinic[,(names(top_expr_clinic) %in% c(top,"RFS","RFS_IND","DFS","DFS_IND"))]
cox_top <-coxph(Surv(DFS,DFS_IND)~.,data=top_clinic)
coeffs_top <- round(as.numeric(summary(cox_top)$coefficient[,1]),3)
coef_list<- summary(cox_top)$coefficients[,1]
cox_coef <- as.numeric(coef_list[top])

score<- data.frame(sample=rownames(top_clinic),stringsAsFactors=FALSE)

for(j in 1:length(top)){
  score1 <- top_clinic[,top[j]]*cox_coef[j]
  score <- cbind(score1, score)
}
score$sample<-NULL
top_expr_clinic$riskscore<-apply(score, 1, sum)
rr<-summary(top_expr_clinic$riskscore)
median <- as.numeric(rr[3])
top_expr_clinic$medianB <- as.numeric(as.numeric(top_expr_clinic$riskscore>=median))
out1 <-top_expr_clinic
ss1 <- survdiff(Surv(DFS,DFS_IND) ~ medianB,data=top_expr_clinic)
pValue1 <- 1-pchisq(ss1$chisq, 1)
cox_run1<- coxph(Surv(DFS,DFS_IND) ~ medianB,,data=top_expr_clinic)
hr1<- as.numeric(summary(cox_run1)$coefficients[2])
concor1 <- cox_run1$concordance[6]
top_result[1,1] <- length(top)
top_result[1,2] <- pValue1
top_result[1,3] <- hr1
top_result[1,4] <- concor1
##validation
dataset<- data.frame(set=c("GSE39582","GSE17538","GSE33113"))
used_gene <-top
out_val<- c()
for ( k in 1:dim(dataset)[1]){
  gse_id<- as.character(dataset[k,])
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
  
  test_top <- dataAll[,(names(dataAll) %in% c(top,"RFS","RFS_IND","DFS","DFS_IND"))]
  
  coef_list<- summary(cox_top)$coefficients[,1]
  cox_coef <- as.numeric(coef_list[top])
  score<- data.frame(sample=rownames(test_top),stringsAsFactors=FALSE)
  for(i in 1:length(used_gene)){
    score1 <- as.numeric(test_top[,used_gene[i]])*cox_coef[i]
    score <- cbind(score1, score)
  }
  score$sample<-NULL
  test_top$riskscore<-apply(score, 1, sum)
  rr <-summary(as.numeric(test_top$riskscore))
  median <- as.numeric(rr[3])
  test_top$medianB <-as.numeric(as.numeric(test_top$riskscore)>=median)
  ss <-survdiff(Surv(RFS,RFS_IND)~medianB,data=test_top)
  pValue <- 1-pchisq(ss$chisq,1)
  cox_run<- coxph(Surv(RFS,RFS_IND) ~ medianB,,data=test_top)
  hr<- as.numeric(summary(cox_run)$coefficients[2])
  concor <-cox_run$concordance[6]
  out_val<- c(out_val,c(dataset[k,],pValue,hr,concor) )
}
top_result[1,5:16]<- out_val
write.csv(top_result,"Result_XGBoostmodel.csv")