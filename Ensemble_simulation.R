
install.packages("mvtnorm")
install.packages("glmnet")
install.packages("SIS")
install.packages("spls")
install.packages("MASS")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("qvalue", version = "3.8")



source("rankprodbounds.R")
library(mvtnorm)
library(glmnet)
library(SIS)
library(spls)
library(parallel)
library(qvalue)
library(MASS)

n=500
p=1000
nn=100


beta1=rep(1,30)
beta0=rep(0,p-length(beta1))
beta=c(beta1,beta0)
beta1 <- which(beta!=0)   ##index variable
beta0 <-which(beta==0)   ####

lSp=rep(0,nn)     
lSn=rep(0,nn)    
lnum=rep(0,nn) 

onum=rep(0,nn)    # num of outliers
power=rep(0,nn) 
Ierror=rep(0,nn) 
Yorden=rep(0,nn)  # Yoden

data=list()
data1=list()

file=list()
file1=list()

txtfile=paste(1:nn,"train.txt",sep="")
txtfile1=paste(1:nn,"test.txt",sep="")   

for(m in 1:nn)
       {
       file[[m]]=paste("D:/2019/simulation/n500p1000eps5%/",txtfile[m],sep="")
       data[[m]]=read.table(file=file[[m]],header=TRUE)
       y=as.matrix(data[[m]][1])
       outlier_indi<-as.matrix(data[[m]][2])
       routliers<-which(outlier_indi==1)
       rinliers<-which(outlier_indi==0)

       x=as.matrix(data[[m]][-c(1,2)])
       file1[[m]]=paste("D:/2019/simulation/n500p1000eps5%/",txtfile1[m],sep="")

       data1[[m]]=read.table(file=file1[[m]],header=TRUE)

       y1=as.matrix(data1[[m]][1])
       x1=as.matrix(data1[[m]][-1])

      ##elastic net

      my.alpha <- seq(0.1,0.9,0.1)     
      nvar.selected.EN <- matrix(0,1,length(my.alpha))
      pred.EN <- matrix(0,dim(x)[1],length(my.alpha))
      MSE.EN <- matrix(0,1,length(my.alpha))

      set.seed(2010)
      foldid.tnbc <- sample(1:10,size=length(y),replace=TRUE)

      for (j in 1:length(my.alpha)){
 
       # Logistic model fitting with 10-fold cross-validation for glmnet:
       fit.EN.cv <- cv.glmnet(x,y,family="binomial",foldid=foldid.tnbc,alpha=my.alpha[j])
 
       var.selected.EN <- which(fit.EN.cv$glmnet.fit$beta[,which(fit.EN.cv$cvm==min(fit.EN.cv$cvm))]!= 0)
       nvar.selected.EN[j] <- length(var.selected.EN)
 
       # Predictions obtained by model i
       pred.EN[,j] <- predict(fit.EN.cv,newx=x[1:n,],s="lambda.min",type="response")
 
       # Mean squared error of prediction (MSE)
      MSE.EN[j] <- mean((y-pred.EN[,j])^2)
     }

     my.alpha[which(MSE.EN == min(MSE.EN))]  

     ##three models for en, SPLS-DA,SGPLS
     n.models <- 3
     nvar.selected <- matrix(0,1,n.models)
     MSE <- matrix(0,1,n.models)    ##predciton error for the train datasets
     MSE1 <- matrix(0,1,n.models)   ##predciton error for the test datasets
     Miscl <- matrix(0,dim(x)[1],n.models) ##MS rate for the train datasets
     Miscl1 <- matrix(0,dim(x)[1],n.models)   ##MS rate for the test datasets
     CookD <- matrix(0,dim(x)[1],n.models)

     i <- 1
      my.alpha.opt <- my.alpha[which(MSE.EN == min(MSE.EN))]
 
      fit.EN.cv <- cv.glmnet(x,y,family="binomial",foldid=foldid.tnbc,alpha=my.alpha[which(MSE.EN == min(MSE.EN))])

      var.selected.EN <- which(fit.EN.cv$glmnet.fit$beta[,which(fit.EN.cv$cvm == min(fit.EN.cv$cvm))] != 0)
      nvar.selected[i] <- length(var.selected.EN)
       assign(paste("var.selected", i, sep = ""),var.selected.EN)
  
      # Predictions obtained by model i:
      pred.EN <- predict(fit.EN.cv,x[1:n,],s="lambda.min",type="response")       
      pred.EN1 <- predict(fit.EN.cv,x1[1:n,],s="lambda.min",type="response")     

      # Mean squared error of prediction (MSE):                                
       MSE[i] <- mean((y-pred.EN)^2)  
       MSE1[i] <- mean((y1-pred.EN1)^2)                                

      # Misclassified individuals based on model i:
       Miscl[,i] <- mean(abs(y- round(pred.EN)))                      
       table(y,round(pred.EN))
       table(y1,round(pred.EN1))
 
      # The Cook’s distance for each individual based on model i:
        V <- diag(as.vector(sqrt(pred.EN*(1-pred.EN))))
        H <- V %*% x[,var.selected.EN] %*% (ginv(t(x[,var.selected.EN]) %*% V %*% x[,var.selected.EN])) %*% t(x[,var.selected.EN]) %*% V             # the hat matrix

      CookD[,i] <-(y- pred.EN)^2*diag(H)/((pred.EN*(-pred.EN+1))*((-diag(H)+1))^2)
 
       # SPLSDA 
      i <- 2
 
      library(spls)
      ## Sparse Partial Least Squares (SPLS) Regression and
      ## Classification (version 2.2-1)
      # Optimizing K and eta by cross-validation
   
      # Optimizing K and eta by cross-validation
       fit.splsda.LOGIT.cv <- cv.splsda(x,y, fold=10, K = c(1:5), eta = c(0.9,0.8,0.7), kappa=0.5, classifier="logistic", scale.x=FALSE,n.core=10)
 
      # fit.splsda.LOGIT.cv$K.opt
      # fit.splsda.LOGIT.cv$eta.opt
 
      # fixing the optimum K and eta obtained by cross-validation (for reproducibility)


      fit.splsda.LOGIT <- splsda(x,y, K = fit.splsda.LOGIT.cv$K.opt, eta=fit.splsda.LOGIT.cv$eta.opt, classifier="logistic", scale.x=FALSE)
      var.selected.splsda.LOGIT <- fit.splsda.LOGIT$A
      nvar.selected[i] <- length(var.selected.splsda.LOGIT)
       pred.splsda.LOGIT <- predict(fit.splsda.LOGIT, x,type = "fit","coefficient", fit.type = "response")  ##orginal samples
       pred.splsda.LOGIT1 <- predict(fit.splsda.LOGIT, x1,type = "fit","coefficient", fit.type = "response")   ##test samples
     

      table(y,round(pred.splsda.LOGIT))
      MSE[i] <- mean((y-pred.splsda.LOGIT)^2)     
      MSE1[i] <- mean((y1-pred.splsda.LOGIT1)^2)  

      Miscl[,i] <- abs(y - round(pred.splsda.LOGIT))  
      Miscl1[,i] <- abs(y1 - round(pred.splsda.LOGIT1))  

      V <- diag(as.vector(sqrt(pred.splsda.LOGIT*(1-pred.splsda.LOGIT))))
      H <- V %*% (fit.splsda.LOGIT$T) %*% (solve(t(fit.splsda.LOGIT$T) %*% V %*% fit.splsda.LOGIT$T)) %*% t(fit.splsda.LOGIT$T) %*% V # the hat matrix
 
      CookD[,i] <- (y - pred.splsda.LOGIT)^2 * diag(H) / ((pred.splsda.LOGIT*(-pred.splsda.LOGIT+1))*((-diag(H)+1))^2)

   
      i <- 3
   
      library(spls)
      library(parallel)
 
     # Optimizing K and eta by cross-validation
      fit.sgpls.LOGIT.cv <- cv.sgpls(x,y, fold=10, K = c(1:5), eta = c(0.9,0.8,0.7), scale.x=FALSE,n.core=10)
      #fit.sgpls.LOGIT.cv$K.opt
      #fit.sgpls.LOGIT.cv$eta.opt
      # fixing the optimum K and eta obtained by cross-validation (for reproducibility)
      fit.sgpls.LOGIT <- sgpls(x,y, K = fit.sgpls.LOGIT.cv$K.opt, eta=  fit.sgpls.LOGIT.cv$eta.opt, scale.x=FALSE)
  
      var.selected.sgpls.LOGIT <- fit.sgpls.LOGIT$A
      nvar.selected[i] <- length(var.selected.sgpls.LOGIT)
       assign(paste("var.selected", i, sep = ""), var.selected.sgpls.LOGIT)
 
      pred.sgpls.LOGIT <- predict(fit.sgpls.LOGIT, x,type="fit",fit.type="response")
      pred.sgpls.LOGIT1 <- predict(fit.sgpls.LOGIT, x1,type="fit",fit.type="response")

      table(y,round(pred.sgpls.LOGIT))
      table(y1,round(pred.sgpls.LOGIT1))

       # Misclassified individuals based on model i:
       Miscl[,i] <- abs(y - round(pred.sgpls.LOGIT))
       Miscl1[,i] <- abs(y1 - round(pred.sgpls.LOGIT1))

       # The Cook’s distance for each individual based on model i:
       V <- diag(as.vector(sqrt(pred.sgpls.LOGIT*(1-pred.sgpls.LOGIT))))
       H <- V %*% (x[,fit.sgpls.LOGIT$A]%*%fit.sgpls.LOGIT$W) %*% (solve(t(x[,fit.sgpls.LOGIT$A]%*%fit.sgpls.LOGIT$W) %*% V %*% x[,fit.sgpls.LOGIT$A]%*              %fit.sgpls.LOGIT$W)) %*% t(x[,fit.sgpls.LOGIT$A]%*%fit.sgpls.LOGIT$W) %*% V
       CookD[,i] <- (y - pred.sgpls.LOGIT)^2 * diag(H) / ((pred.sgpls.LOGIT*(-pred.sgpls.LOGIT+1))*((-diag(H)+1))^2)

       # The rank product:
       rank.matrix <- apply(CookD,2, function(CookD) rank(-(CookD),ties.method = "first"))
       rank.product <- apply(rank.matrix, 1, prod)

       rho <-rank.product
       n <- dim(x)[1]
       k <- dim(rank.matrix)[2]

       # The p-values:
       pvalues <- as.vector(rankprodbounds(rho,n,k,Delta ='geometric'))

        # The q-values:
       library(qvalue)
       qobj <- qvalue(pvalues)
       qvalues <- qobj$qvalues
 
      # Misclassifications (%) by the three models
       Miscl.percent <- apply(Miscl,1, function(Miscl) round((sum(Miscl)/k)*100))

      id <- as.vector(seq(1:dim(x)[1]))
      outliers.rank <- as.data.frame(cbind(id, rank.matrix,rank.product,pvalues,qvalues,Miscl.percent))

      # influent patients
      TNBC.influent <- outliers.rank[which(outliers.rank[,7]<0.05),1]
       onum[m]=length(TNBC.influent)   ###异常点个数

      ###outlier detection
      odent=TNBC.influent
      ident=(1:n)[-odent]
      TPresult<- length(intersect(odent,routliers)) #TP
      FPresult<- length(intersect(odent,rinliers)) #FP
      TNresult<- length(intersect(ident,rinliers)) #TN
      FNresult<- length(intersect(ident,routliers)) #FN

      power[m]=TPresult/length(routliers)   
      Ierror[m]=FPresult/length(rinliers) 
      Yorden[m]=power[m]-Ierror[m]

      ###the intersection of variables selection 

      a=intersect(var.selected.sgpls.LOGIT,var.selected.splsda.LOGIT )
      var.selected.intersect=intersect(a,as.vector(var.selected.EN))
      var.unselected=(1:p)[-var.selected.intersect]

     lnum[m]=length(var.selected.intersect)  

     shi1=which(beta!=0)
     shi2=which(beta==0)   
     
     TP<- length(intersect(shi1, var.selected.intersect)) #TP
     FP<- length(intersect(shi2,var.selected.intersect)) #FP
     TN<- length(intersect(shi2,var.unselected)) #TN
     FN<- length(intersect(shi1,var.unselected)) #FN
  
     lSn[m]=TP/(TP+FN)
         
     if((TP+FP)==0)
        {lSp[m]=0} else
      {lSp[m]=TP/(TP+FP)}

 }
 
   lSnSp=sqrt(lSn*lSp)  #GM
 
paste("Marta的结果：")
paste("mean(onum) is",mean(onum))
paste("mean(power) is",mean(power))
paste("mean(Ierror) is",mean(Ierror))
paste("mean(Yorden) is",mean(Yorden))
paste("mean(num) is",mean(lnum))
paste("mean(fdr) is",mean(1-lSp))
paste("mean(Sn) is",mean(lSn))
paste("mean(SnSp) is",mean(lSnSp))

paste("sd(onum) is",sd(onum))
paste("sd(power) is",sd(power))
paste("sd(Ierror) is",sd(Ierror))
paste("sd(Yorden) is",sd(Yorden))
paste("sd(num) is",sd(lnum))
paste("sd(fdr) is",sd(1-lSp))
paste("sd(Sn) is",sd(lSn))
paste("sd(SnSp) is",sd(lSnSp))

write(onum,"D:\\2019\\simulation\\n500p1000eps5%\\m_onum.xls") 
write(power,"D:\\2019\\simulation\\n500p1000eps5%\\m_power.xls") 
write(Ierror,"D:\\2019\\simulation\\n500p1000eps5%\\m_Ierror.xls") 
write(Yorden,"D:\\2019\\simulation\\n500p1000eps5%\\m_Yorden.xls") 
write(lnum,"D:\\2019\\simulation\\n500p1000eps5%\\m_lnum.xls")
write(1-lSp,"D:\\2019\\simulation\\n500p1000eps5%\\m_fdr.xls")
write(lSn,"D:\\2019\\simulation\\n500p1000eps5%\\m_lSn.xls")
write(lSnSp,"D:\\2019\\simulation\\n500p1000eps5%\\m_lSnSp.xls")


    
   

