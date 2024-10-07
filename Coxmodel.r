#set working directory
setwd("./path/to/Colon_data")

library("survival")
TCGA_COAD<-read.csv("TCGA_COAD_clinic_clean.csv")
Used_genes <- data.frame(colnames(TCGA_COAD)[2:76])
DFS<- TCGA_COAD$DFS
DFS_IND<- TCGA_COAD$DFS_IND
results<- data.frame(matrix(ncol =5, nrow = 0)) 
colnames(results)<- c("gene","cox_p","hr","lower","upper")
for (i in 1:length(Used_genes[,1])){
  gene<- Used_genes[i,1]
  genExp<- TCGA_COAD[,gene]
  
  cc<- coxph(Surv(DFS, DFS_IND) ~ genExp)
  cox_p<- round(as.numeric(summary(cc)$coefficient[5]),4)
  hr <- round(as.numeric(summary(cc)$coefficient[2]),4)
  lower<- round(as.numeric(summary(cc)$conf.int[3]),4)
  upper<- round(as.numeric(summary(cc)$conf.int[4]),4)
  
  results[i,1] <- gene
  results[i,2] <- cox_p
  results[i,3] <- hr
  results[i,4] <- lower
  results[i,5] <- upper
}
Sig_res <- results[which(results$cox_p < 0.05),]

top_result <- data.frame(matrix(NA,nrow=0,ncol=16))
names(top_result) <- c("variable","train_pvalue","train_hr","train_concor",rep(c("test_dataset","test_pvalue","test_hr","test_concor"),3))
top <-Sig_res$gene
top_exp <- data.frame(TCGA_COAD[,which((colnames(TCGA_COAD) %in% c(top,"DFS","DFS_IND")))])
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

write.csv(top_result,"Result_TCGA_COAD.csv")
