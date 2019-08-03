#n150,p200,epso0.15,xy异常

##注意 ，用的还是qvlue, 弹性网，用的是Pvalue


install.packages("mvtnorm")
install.packages("glmnet")
install.packages("SIS")
install.packages("spls")
install.packages("enetLTS")

library(mvtnorm)
library(glmnet)
library(SIS)
library(enetLTS)

n=500
p=1000
nn=100 #repetition

beta1=rep(1,30)
beta0=rep(0,p-length(beta1))
beta=c(beta1,beta0)
beta1 <- which(beta!=0)
beta0 <-which(beta==0)  

lt2Sp=rep(0,nn)      
lt2Sn=rep(0,nn)    
lt2num=rep(0,nn)  
lt2error=rep(0,nn)  ##MS rate
lt2rmspe=rep(0,nn)  ##MSE

onum=rep(0,nn)    #num of outliers
power=rep(0,nn)  #Sn
Ierror=rep(0,nn)  #FPR
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
       outlier_indi=as.vector(outlier_indi)
       routliers<-which(outlier_indi==1)
       rinliers<-which(outlier_indi==0)
    
      x=as.matrix(data[[m]][-c(1,2)])
      file1[[m]]=paste("D:/2019/simulation/n500p1000eps5%/",txtfile1[m],sep="")
      data1[[m]]=read.table(file=file1[[m]],header=TRUE)

      y1=as.matrix(data1[[m]][1])
      x1=as.matrix(data1[[m]][-1])

      xsd=standardize(x)  
      xsd1=standardize(x1)  
   
      ####enetLTS######

      fit <- enetLTS(x,y,family="binomial",plot=FALSE) 
      ll2=fit$coefficients  
      lt2predict2=as.matrix(fit$fitted.values)   ### for train dataset
      lt2predict=exp(xsd1%*%ll2)/(1+exp(xsd1%*%ll2))  ###for test data
      lt2num[m]=length(which(ll2!=0))   
   
      TP=length(which(ll2[beta1]!=0))    
      FN=length(which(ll2[beta1]==0))   
      TN=length(which(ll2[beta0]==0))    
      FP=length(which(ll2[beta0]!=0)) 
      lt2Sn[m]=TP/(TP+FN)
      if((TP+FP)==0)
        {lt2Sp[m]=0} else
       {lt2Sp[m]=TP/(TP+FP)}
       
      lt2rmspe[m]=sqrt(sum((y1-lt2predict)^2)/n) #prediction error
 
      lt2yhat=rep(0,n)
      lt2yhat[which(lt2predict>=0.5)]=1
      lt2error[m]= (n-length(which(y1==lt2yhat)))/n    #MS rate 

      ##outliers detection
      
      outliers=(1:n)[fit$wt==0]   
      onum[m]=length(outliers) 

      outlier_ident=rep(0,n)
      outlier_ident[fit$wt==0]=1          

      odent=outliers
      ident=(1:n)[-odent]

      TPresult<- length(intersect(odent,routliers)) #TP*
      FPresult<- length(intersect(odent,rinliers)) #FP*
      TNresult<- length(intersect(ident,rinliers)) #TN*
      FNresult<- length(intersect(ident,routliers)) #FN*

       power[m]=TPresult/length(routliers)
      Ierror[m]=FPresult/length(rinliers)  
       Yorden[m]=power[m]-Ierror[m]
   
 }

 lt2SnSp=sqrt(lt2Sn*lt2Sp) #GM

paste("enetLTS的结果：")
paste("mean(onum) is",mean(onum))
paste("mean(power) is",mean(power))
paste("mean(Ierror) is",mean(Ierror))
paste("mean(Yorden) is",mean(Yorden))
paste("mean(num) is",mean(lt2num))
paste("mean(fdr) is",mean(1-lt2Sp))
paste("mean(Sn) is",mean(lt2Sn))
paste("mean(SnSp) is",mean(lt2SnSp))
paste("mean(MSE) is",mean(lt2rmspe))
paste("mean(Miscl) is",mean(lt2error))

paste("sd(onum) is",sd(onum))
paste("sd(power) is",sd(power))
paste("sd(Ierror) is",sd(Ierror))
paste("sd(Yorden) is",sd(Yorden))
paste("sd(num) is",sd(lt2num))
paste("sd(fdr) is",sd(1-lt2Sp))
paste("sd(Sn) is",sd(lt2Sn))
paste("sd(SnSp) is",sd(lt2SnSp))
paste("sd(MSE) is",sd(lt2rmspe))
paste("sd(Miscl) is",sd(lt2error))

write(onum,"D:\\2019\\simulation\\n500p1000eps5%\\k_onum.xls") 
write(power,"D:\\2019\\simulation\\n500p1000eps5%\\k_power.xls") 
write(Ierror,"D:\\2019\\simulation\\n500p1000eps5%\\k_Ierror.xls") 
write(Yorden,"D:\\2019\\simulation\\n500p1000eps5%\\k_Yorden.xls") 
write(lt2num,"D:\\2019\\simulation\\n500p1000eps5%\\k_lnum.xls")
write(1-lt2Sp,"D:\\2019\\simulation\\n500p1000eps5%\\k_fdr.xls")
write(lt2Sn,"D:\\2019\\simulation\\n500p1000eps5%\\k_lSn.xls")
write(lt2SnSp,"D:\\2019\\simulation\\n500p1000eps5%\\k_lSnSp.xls")
write(lt2rmspe^2,"D:\\2019\\simulation\\n500p1000eps5%\\k_MSE1.xls")
write(lt2error,"D:\\2019\\simulation\\n500p1000eps5%\\k_Miscl.xls")  



