#8-eight_survival curve
#set working directory
setwd("C:\\Users\\xiuxi\\Desktop\\survplot")

library(survival)

train_data <- read.csv("8genes_TCGA_dataALL.csv")
ss <-survdiff(Surv(DFS,DFS_IND)~medianB,data=train_data)
pValue <- 1-pchisq(ss$chisq,1)
pValue <-round(pValue,6)
fit <-survfit(Surv(DFS,DFS_IND)~medianB,data=train_data)

pdf(paste("8genes_","TCGA","_survplot.pdf"))  

plot(fit,lty=1,lwd=2,col=c("blue","red"),mark=3,yscale=100,las=1,
     xlab="time(months)",ylab="survival probability")
legend(0,0.25,c("low-risk","high-risk"),lty=1,lwd=2,col=c("blue","red"),
       box.lty=0)
text(x=22,y=0.05,labels=paste("log-rank P =",pValue))

dev.off()


#test datasets for survplot
dataset <- data.frame(set=c("GSE33113","GSE17538","GSE39582"))

for (k in 1:dim(dataset)[1]){
  gse_id<- as.character(dataset[k,])
  test_data<- read.csv(paste("8genes_",gse_id,"_dataALL.csv",sep=""),header=T,row.names=1,
                       stringsAsFactors=FALSE)
  
  pdf(paste("8genes_",gse_id,"_survplot.pdf"))  
  
  ss1 <-survdiff(Surv(RFS,RFS_IND)~medianB,data=test_data)
  pValue1 <- 1-pchisq(ss1$chisq,1)
  pValue1 <-round(pValue1,3)
  fit1 <-survfit(Surv(RFS,RFS_IND)~medianB,data=test_data)
  
  p3 <- plot(fit1, lty = 1, lwd=2, col=c('blue','red'), mark=3, yscale=100,las=1,
             xlab='time(months)',ylab='survival probability')
  
  legend(0,0.25,c("Low risk","High risk"), lty=1, lwd=2, col=c('blue','red'), box.lty=0)
  text(x= 22, y= 0.05, labels=paste("Log -rank P =", pValue1))
  
  dev.off()
}