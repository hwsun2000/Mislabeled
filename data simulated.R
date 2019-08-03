install.packages("mvtnorm")
install.packages("glmnet")
install.packages("SIS")
install.packages("parallel")
install.packages("R.matlab")

library(mvtnorm)
library(glmnet)
library(SIS)
library(R.matlab)


set.seed(101) # for reproducibility

p=1000  
n=500   #sample size
epsilon=0.05  #proprotion of outliers
epsilon1=epsilon*(1/3)  #y=0 misclassified to y=1
epsilon2=epsilon*(2/3)  #y=1 misclassified to y=0

nn=100  #replication

#coeficient

beta1=rep(1,30)
beta0=rep(0,p-length(beta1))
beta=c(beta1,beta0)

ssigma =0.9   #corelation coeficient
sigma =ssigma^t(sapply(1:p, function(i, j) abs(i-j), 1:p))#covariance matrix



outPath <- "D:/2019/simulation/n500p1000eps5%" 
out_fileName <- sapply(names(data),function(x){
                    paste(x, ".txt", sep='')}) 
out_filePath  <- sapply(out_fileName, function(x){
                     paste(outPath,x,sep='/')}) 

data<-list()

for(i in 1:nn)
    {
       
       x=rmvnorm(n,sigma=sigma) # generate x 

       feta = x%*%beta; 
       fprob = exp(feta)/(1+exp(feta))
       y=rbinom(n,1,fprob)       
      
       index1=which(y==1) 
       index0=which(y==0)
       n1=length(index1)
       n0=length(index0)

       k0=sample(index0,floor(n*epsilon1))  
       k1=sample(index1,floor(n*epsilon2))   

       y[k0]=1 
       y[k1]=0 

       outlier=c(k0,k1)
        ## x[outlier]=x[outlier]+20

       outlier_indi<-rep(0,n)
       outlier_indi[outlier]=1
  
     
       filnm <- paste("D:/2019/simulation/n500p1000eps5%/x",i,".mat",sep="")
       writeMat(filnm, x = x)

        filnm <- paste("D:/2019/simulation/n500p1000eps5%/y",i,".mat",sep="")
       writeMat(filnm, y = y)

         filnm <- paste("D:/2019/simulation/n500p1000eps5%/fd",i,".mat",sep="")
       writeMat(filnm, fd = outlier_indi)

       data[[i]]=as.data.frame(cbind(y,outlier_indi,x)) 

   }

      file=list()
      txtfile=paste(1:nn,"train.txt",sep="")


for(i in 1:nn)
   {
       file[[i]]=paste("D:/2019/simulation/n500p1000eps5%/",txtfile[i],sep="")
       write.table(data[[i]],file=file[[i]])
    }


      ##test data
for(i in 1:nn)
   {
      x1=rmvnorm(n, sigma=sigma)
      feta = x1%*%beta; 
      fprob = exp(feta)/(1+exp(feta))
      y1=rbinom(n, 1, fprob)     
      data[[i]]=as.data.frame(cbind(y1,x1)) 
    }


file=list()
txtfile=paste(1:nn,"test.txt",sep="")

for(i in 1:nn)
    {
    file[[i]]=paste("D:/2019/simulation/n500p1000eps5%/",txtfile[i],sep="")
    write.table(data[[i]],file=file[[i]])
    }









  