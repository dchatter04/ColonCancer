#set working directory 
setwd("./path/to/Colon_data")
library(glmnet)
library(rms)
library(VIM)
library(survival)
library(Hmisc)
library(data.table)

data <- read.csv("TCGA_COAD_clinic_clean.csv",header=T,row.names=1)
data1 <- data[,c(1:75,78:79)]
data2 <- na.omit(data1)
#Extract the disease free survival information (DFS) and the disease free survival status indicator (DFS_IND as 0 and 1, where 1 indicates survival)
data2$DFS <- as.numeric(data2$DFS)
data3 <- data2[which(data2$DFS!=0),]
data3$DFS <- as.double(data3$DFS)
data3$DFS_IND <- as.double(data3$DFS_IND)

x <- data.matrix(data3[,c(1:75)]) #The gene expression matrix for the TCGA data
y <- data.matrix(Surv(time=data3$DFS,event=data3$DFS_IND)) #The disease free survival information matrix
fit <- glmnet(x,y,family="cox",alpha=1)
cv_fit <- cv.glmnet(x, y,family="cox",nfolds=10,type.measure="C")
coef <- coef(fit,s=cv_fit$lambda.1se)
active.Index <- which(coef !=0)
length(active.Index)
active.coef <- coef[active.Index]
get_coe <- function(fit,s){
  coef <- coef(fit, s =cv_fit$lambda.1se)
  active.Index <- which(coef != 0)
  active.coef <- coef[active.Index]
  re <- data.frame(rownames(coef)[active.Index],active.coef)
  re <- data.table('var_names'=rownames(coef)[active.Index],
                   'coef'=active.coef)
  re$expcoef <- exp(re$coef)
  return(re[order(expcoef)])
}
b <- get_coe(cv_fit,cv_fit$lambda.1se)
b1 <- b[order(b$coef,decreasing=T),]
b2 <- head(b1,20)
var_select <- b2$var_names
top_result <- data.frame(matrix(NA,nrow=0,ncol=16))
names(top_result) <- c("variable","train_pvalue","train_hr","train_concor",rep(c("test_dataset","test_pvalue","test_hr","test_concor"),3))
  top <-var_select[1:15]
  top_exp <- data.frame(data3[,which((colnames(data3) %in% c(top,"DFS","DFS_IND")))])
  cox_top<- coxph(Surv(DFS,DFS_IND) ~.,data=top_exp) 
  coeffs_top <- round(as.numeric(summary(cox_top)$coefficient[,1]),4)
  coef_list<- summary(cox_top)$coefficients[,1]
  cox_coef <- as.numeric(coef_list[top])
  score<- data.frame(sample=rownames(top_exp),stringsAsFactors=FALSE)
  for(j in 1:length(top)){
    score1 <- top_exp[,top[j]]*cox_coef[j]
    score <- cbind(score1, score)
  }
  score$sample<-NULL
  top_exp$riskscore<-apply(score, 1, sum)
  rr<-summary(top_exp$riskscore)
  median <- as.numeric(rr[3])
  top_exp$medianB <- as.numeric(as.numeric(top_exp$riskscore>=median))
  ss1 <- survdiff(Surv(DFS,DFS_IND) ~ medianB,data=top_exp)
  pValue1 <- 1-pchisq(ss1$chisq, 1)
  cox_run1<- coxph(Surv(DFS,DFS_IND) ~ medianB,,data=top_exp)
  hr1<- as.numeric(summary(cox_run1)$coefficients[2])
  concor1 <-cox_run1$concordance[6]
  top_result[1,1] <- length(top)
  top_result[1,2] <- pValue1
  top_result[1,3] <- hr1
  top_result[1,4] <- concor1
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
  

write.csv(top_result,"Result_Lassomodel.csv")