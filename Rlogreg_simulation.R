
install.packages("mvtnorm")
install.packages("glmnet")

install.packages("R.matlab")
install.packages("SIS")

library(mvtnorm)
library(glmnet)
library(SIS)
library(R.matlab)

n=500
p=1000
nn=100

beta1=rep(1,30)
beta0=rep(0,p-length(beta1))
beta=c(beta1,beta0)
beta1 <- which(beta!=0)##index variable
beta0 <-which(beta==0)  

lSp=rep(0,nn)   
lSn=rep(0,nn)     
lnum=rep(0,nn)  
MSE1=rep(0,nn)   
Miscl1=rep(0,nn) 

onum=rep(0,nn)   
power=rep(0,nn)  
Ierror=rep(0,nn)  
Yorden=rep(0,nn)  

data=list()
data1=list()

file=list()
file1=list()

txtfile=paste(1:nn,"train.txt",sep="")
txtfile1=paste(1:nn,"test.txt",sep="")   

w<- readMat("D:/2019/simulation/n500p1000eps5%/w.mat")
w<-as.matrix(w$w)

    
for(m in 1:nn)
   {
    
   file[[m]]=paste("D:/2019/simulation/n500p1000eps5%/",txtfile[m],sep="")
   data[[m]]=read.table(file=file[[m]],header=TRUE)

    y=as.matrix(data[[m]][1])
    outlier_indi<-as.matrix(data[[m]][2])

    routliers<-which(outlier_indi==1)
    rinliers<-which(outlier_indi==0)
    
    x=as.matrix(data[[m]][-c(1,2)])
    xsd=standardize(x)  

    ###Accuracy of variables selection
    
    ww=w[-1,m]
    wpredict=exp(xsd%*%ww)/(1+exp(xsd%*%ww))

    lnum[m]=length(which(ww!=0))   #model size of lasso
     
    TP=length(which(ww[beta1]!=0)) 
    FN=length(which(ww[beta1]==0)) 
    TN=length(which(ww[beta0]==0))    
    FP=length(which(ww[beta0]!=0))  
    lSn[m]=TP/(TP+FN)
      
    if((TP+FP)==0)
    {lSp[m]=0} else
    {lSp[m]=TP/(TP+FP)}

     #####misclassification#############

    file1[[m]]=paste("D:/2019/simulation/n500p1000eps5%/",txtfile1[m],sep="")   ###test data
 
    data1[[m]]=read.table(file=file1[[m]],header=TRUE)   

    y1=as.matrix(data1[[m]][1])
    x1=as.matrix(data1[[m]][-1]) 
    xsd1=standardize(x1)                   

    wpredict1=exp(xsd1%*%ww)/(1+exp(xsd1%*%ww))   

    # Mean squared error of prediction (MSE):                                
      
    MSE1[m] <- mean((y1-wpredict1)^2)                                

    #misclassification rate

    Miscl1[m] <- mean(abs(y1- round(wpredict1)));                      
 
    ######outliers detection
      
    yhat<-rep(0,dim(x)[1])       
    yhat[which(wpredict>=0.5)]=1           
    outliers<-which(y!=yhat)      ##outliers are the mislabelled sample.

    onum[m]=length(outliers)   ###Number of outliers

    ##########outliers detection########################

    odent=outliers
    ident=(1:n)[-odent]
    TPresult<- length(intersect(odent,routliers)) #TP
    FPresult<- length(intersect(odent,rinliers)) #FP
    TNresult<- length(intersect(ident,rinliers)) #TN
    FNresult<- length(intersect(ident,routliers)) #FN

    power[m]=TPresult/length(routliers)
    Ierror[m]=FPresult/length(rinliers)  
    Yorden[m]=power[m]-Ierror[m] 

}

lSnSp=sqrt(lSn*lSp)

paste("n100p200eps5%,RlogregµÄ½á¹û£º")
paste("mean(onum) is",mean(onum))
paste("mean(power) is",mean(power))
paste("mean(Ierror) is",mean(Ierror))
paste("mean(Yorden) is",mean(Yorden))
paste("mean(lnum) is",mean(lnum))
paste("mean(fdr) is",mean(1-lSp))
paste("mean(lSn) is",mean(lSn))
paste("mean(lSnSp) is",mean(lSnSp))
paste("mean(MSE1) is",mean(MSE1))
paste("mean(Miscl1) is",mean(Miscl1))


paste("sd(onum) is",sd(onum))
paste("sd(power) is",sd(power))
paste("sd(Ierror) is",sd(Ierror))
paste("sd(Yorden) is",sd(Yorden))
paste("sd(lnum) is",sd(lnum))
paste("sd(fdr) is",sd(1-lSp))
paste("sd(lSn) is",sd(lSn))
paste("sd(lSnSp) is",sd(lSnSp))
paste("sd(MSE1) is",sd(MSE1))
paste("sd(Miscl1) is",sd(Miscl1))

  
write(onum,"D:\\2019\\simulation\\n500p1000eps5%\\r_onum.xls") 
write(power,"D:\\2019\\simulation\\n500p1000eps5%\\r_power.xls") 
write(Ierror,"D:\\2019\\simulation\\n500p1000eps5%\\r_Ierror.xls") 
write(Yorden,"D:\\2019\\simulation\\n500p1000eps5%\\r_Yorden.xls") 
write(lnum,"D:\\2019\\simulation\\n500p1000eps5%\\r_lnum.xls")
write(1-lSp,"D:\\2019\\simulation\\n500p1000eps5%\\r_fdr.xls")
write(lSn,"D:\\2019\\simulation\\n500p1000eps5%\\r_lSn.xls")
write(lSnSp,"D:\\2019\\simulation\\n500p1000eps5%\\r_lSnSp.xls")
write(MSE1,"D:\\2019\\simulation\\n500p1000eps5%\\r_MSE1.xls")
write(Miscl1,"D:\\2019\\simulation\\n500p1000eps5%\\r_Miscl1.xls")  



