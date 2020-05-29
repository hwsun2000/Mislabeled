 #Ensemble outlier detection
 #TNBC original data

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("qvalue", version = "3.8")


source("rankprodbounds.R")
library(mvtnorm)
library(glmnet)
library(spls)
library(parallel)
library(qvalue)
        
data.plus.race.age2<-read.table("data.plus.race.age2.txt",header=TRUE)
x_origin924<-as.matrix(data.plus.race.age2[,1:19690])
y<-as.matrix(data.plus.race.age2[,19691])

ER_PR_HER<-read.table("ER_PR_HER.txt",header=TRUE)


set.seed(2010)
 
my.alpha <- seq(0.1,0.9,0.1)
 
## Optimizing alpha and lambda for logistic regression with Elastic net regularization
 
nvar.selected.EN <- matrix(0,1,length(my.alpha))
pred.EN <- matrix(0,dim(x_origin924)[1],length(my.alpha))
MSE.EN <- matrix(0,1,length(my.alpha))
 
# assigning samples to folds, to be used in cross-validation when tunning alpha
set.seed(2010)
foldid.tnbc <- sample(1:10,size=length(y),replace=TRUE)
 
for (j in 1:length(my.alpha)){
 
  # Logistic model fitting with 10-fold cross-validation for glmnet:
  fit.EN.cv <- cv.glmnet(as.matrix(x_origin924),as.factor(y),family="binomial",foldid=foldid.tnbc,alpha=my.alpha[j])
 
  var.selected.EN <- which(fit.EN.cv$glmnet.fit$beta[,which(fit.EN.cv$cvm == min(fit.EN.cv$cvm))] != 0)
  nvar.selected.EN[j] <- length(var.selected.EN)
 
  # Predictions obtained by model i
  pred.EN[,j] <- predict(fit.EN.cv,as.matrix(x_origin924),s="lambda.min",type="response")
 
  # Mean squared error of prediction (MSE)
  MSE.EN[j] <- mean((y-pred.EN[,j])^2)
}
 
MSE.EN   ##sun
 [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
[1,] 0.0256318 0.02474779 0.02777295 0.02860052 0.02968559 0.02947821
           [,7]       [,8]       [,9]
[1,] 0.02927494 0.02782152 0.02763202


nvar.selected.EN  #sun
 [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
[1,]  389  248  175  140  116   99   88   86   77


my.alpha[which(MSE.EN == min(MSE.EN))]  #sun
[1] 0.2
## Model ensemble
 
n.models <- 3
nvar.selected <- matrix(0,1,n.models)
MSE <- matrix(0,1,n.models)
Miscl <- matrix(0,dim(x_origin924)[1],n.models)
CookD <- matrix(0,dim(x_origin924)[1],n.models)
 
# logistic regression with Elastic net regularization (LOGIT-EN)
 
i <- 1
 
my.alpha.opt <- my.alpha[which(MSE.EN == min(MSE.EN))]
fit.EN.cv <- cv.glmnet(as.matrix(x_origin924),as.factor(y),family="binomial",foldid=foldid.tnbc,alpha=my.alpha[which(MSE.EN == min(MSE.EN))])
 
var.selected.EN <- which(fit.EN.cv$glmnet.fit$beta[,which(fit.EN.cv$cvm == min(fit.EN.cv$cvm))] != 0)
nvar.selected[i] <- length(var.selected.EN)
 
assign(paste("var.selected", i, sep = ""), var.selected.EN)
 
# Predictions obtained by model i:
pred.EN <- predict(fit.EN.cv,as.matrix(x_origin924),s="lambda.min",type="response")
 
# Mean squared error of prediction (MSE):
MSE[i] <- mean((y-pred.EN)^2)
 
# Misclassified individuals based on model i:
Miscl[,i] <- abs(y - round(pred.EN))
 

table(y,round(pred.EN))  ##Sun
y     0   1
  0 756  15
  1  11 142
> 




# The Cook’s distance for each individual based on model i:
V <- diag(as.vector(sqrt(pred.EN*(1-pred.EN))))
H <- V %*% x_origin924[,var.selected.EN] %*% (solve(t(x_origin924[,var.selected.EN]) %*% V %*% x_origin924[,var.selected.EN])) %*% t(x_origin924[,var.selected.EN]) %*% V # the hat matrix
 
CookD[,i] <- (y - pred.EN)^2 * diag(H) / ((pred.EN*(-pred.EN+1))*((-diag(H)+1))^2)
 
 
# SPLSDA
 
i <- 2
 
library(spls)
## Sparse Partial Least Squares (SPLS) Regression and
## Classification (version 2.2-1)
# Optimizing K and eta by cross-validation

fit.splsda.LOGIT.cv <- cv.splsda(as.matrix(x_origin924),as.factor(y), fold=10, K = c(1:5), eta = c(0.9,0.8,0.7), kappa=0.5, classifier="logistic", scale.x=FALSE,n.core=10)
 
#fit.splsda.LOGIT.cv$K.opt
#fit.splsda.LOGIT.cv$eta.opt
 
##sun   optimal parameters:  eta = 0.9, K = 5

# fixing the optimum K and eta obtainedby cross-validation (for reproducibility)
fit.splsda.LOGIT <- splsda(as.matrix(x_origin924),as.factor(y), K = fit.splsda.LOGIT.cv$K.opt, eta=fit.splsda.LOGIT.cv$eta.opt, classifier="logistic", scale.x=FALSE)
 
var.selected.splsda.LOGIT <- fit.splsda.LOGIT$A
nvar.selected[i] <- length(var.selected.splsda.LOGIT)
 
assign(paste("var.selected", i, sep = ""), var.selected.splsda.LOGIT)
 
# Predictions obtained by model i:
pred.splsda.LOGIT <- predict(fit.splsda.LOGIT, as.matrix(x_origin924),type = "fit","coefficient", fit.type = "response")
 
table(y,round(pred.splsda.LOGIT))  #sun
y     0   1
  0 733  38
  1  22 131

# Mean squared error of prediction (MSE):
MSE[i] <- mean((y-pred.splsda.LOGIT)^2)
 
# Misclassified individuals based on model i:
Miscl[,i] <- abs(y - round(pred.splsda.LOGIT))
 
# The Cook’s distance for each individual based on model i:
V <- diag(as.vector(sqrt(pred.splsda.LOGIT*(1-pred.splsda.LOGIT))))
H <- V %*% (fit.splsda.LOGIT$T) %*% (solve(t(fit.splsda.LOGIT$T) %*% V %*% fit.splsda.LOGIT$T)) %*% t(fit.splsda.LOGIT$T) %*% V # the hat matrix
 
CookD[,i] <- (y - pred.splsda.LOGIT)^2 * diag(H) / ((pred.splsda.LOGIT*(-pred.splsda.LOGIT+1))*((-diag(H)+1))^2)
 
# SGPLS
 
i <- 3
 
library(spls)
library(parallel)
 
# Optimizing K and eta by cross-validation
fit.sgpls.LOGIT.cv <- cv.sgpls(as.matrix(x_origin924),y, fold=10, K = c(1:5), eta = c(0.9,0.8,0.7), scale.x=FALSE,n.core=10)
 
# fixing the optimum K and eta obtained by cross-validation (for reproducibility)
  fit.sgpls.LOGIT <- sgpls(as.matrix(x_origin924),y, K = fit.sgpls.LOGIT.cv$K.opt, eta=  fit.sgpls.LOGIT.cv$eta.opt, scale.x=FALSE)
 
var.selected.sgpls.LOGIT <- fit.sgpls.LOGIT$A
nvar.selected[i] <- length(var.selected.sgpls.LOGIT)
 
assign(paste("var.selected", i, sep = ""), var.selected.sgpls.LOGIT)
 
# Predictions obtained by model i:
pred.sgpls.LOGIT <- predict(fit.sgpls.LOGIT, as.matrix(x_origin924),type="fit",fit.type="response")
 
table(y,round(pred.sgpls.LOGIT))  #sun

##y     0   1
  0 741  30
  1  24 129


# Mean squared error of prediction (MSE):
MSE[i] <- mean((y-pred.sgpls.LOGIT)^2)
 
# Misclassified individuals based on model i:
Miscl[,i] <- abs(y - round(pred.sgpls.LOGIT))
 
# The Cook’s distance for each individual based on model i:
V <- diag(as.vector(sqrt(pred.sgpls.LOGIT*(1-pred.sgpls.LOGIT))))
H <- V %*% (x_origin924[,fit.sgpls.LOGIT$A]%*%fit.sgpls.LOGIT$W) %*% (solve(t(x_origin924[,fit.sgpls.LOGIT$A]%*%fit.sgpls.LOGIT$W) %*% V %*% x_origin924[,fit.sgpls.LOGIT$A]%*%fit.sgpls.LOGIT$W)) %*% t(x_origin924[,fit.sgpls.LOGIT$A]%*%fit.sgpls.LOGIT$W) %*% V
 
CookD[,i] <- (y - pred.sgpls.LOGIT)^2 * diag(H) / ((pred.sgpls.LOGIT*(-pred.sgpls.LOGIT+1))*((-diag(H)+1))^2)
 
# Mean squared error for the three models:
MSE  #sun
     [,1]       [,2]      [,3]
[1,] 0.02474779 0.04929222 0.1029172



# Misclassifications for the three models:
apply(Miscl,2,sum)  #sun
[1] 11 60 55

# Number of variables selected by the three models:
nvar.selected   #sun

   [,1] [,2] [,3]
[1,]  175   22   33


# Variables selected by the three models:
var.selected=list(var.selected1,var.selected2,var.selected3) 

# The rank product:
rank.matrix <- apply(CookD,2, function(CookD) rank(-(CookD),ties.method = "first"))
rank.product <- apply(rank.matrix, 1, prod)
 
rho <-rank.product
n <- dim(x_origin924)[1]
k <- dim(rank.matrix)[2]

 #####到这一步了
##source("~/Heskes_pvalues.R")

source("rankprodbounds.R")
 
# The p-values:
pvalues <- as.vector(rankprodbounds(rho,n,k,Delta ='geometric'))
 
# The q-values:
library(qvalue)
qobj <- qvalue(pvalues)
qvalues <- qobj$qvalues
 
# Misclassifications (%) by the three models
Miscl.percent <- apply(Miscl,1, function(Miscl) round((sum(Miscl)/k)*100))



id <- as.vector(seq(1:dim(x_origin924)[1]))
 
#the list of identified outliers, where miscl.percent is the percentage of misclassification in the three models.

outliers.rank <- as.data.frame(cbind(id, rank.matrix,rank.product,pvalues,qvalues,Miscl.percent))
 


outliers.rank <- outliers.rank[order(qvalues),]
 
# influent patients.outliers.rank[,7] is the qvalue.



TNBC.influent <- outliers.rank[which(outliers.rank[,7] < 0.05),1]

outliers.rank[which(outliers.rank[,7] < 0.05),]
  id  V2 V3  V4 rank.product      pvalues     qvalues Miscl.percent
654 654   6 26   2          312 8.903763e-06 0.004113539           100
869 869   2 50   3          300 8.461619e-06 0.004113539           100
874 874  24  1  21          504 1.651119e-05 0.005085448            67
808 808   7  9  12          756 2.763930e-05 0.006384678           100
25   25  12 18   8         1728 7.428429e-05 0.007626520            67
35   35  21 51   1         1071 4.243684e-05 0.007626520           100
59   59  30  3  18         1620 6.894211e-05 0.007626520            67
705 705   9  6  31         1674 7.160895e-05 0.007626520           100
804 804  23 12   6         1656 7.071904e-05 0.007626520           100
64   64  13 13  19         3211 1.500056e-04 0.011481640            67
65   65   4 33  26         3432 1.615382e-04 0.011481640           100
442 442   5 22  30         3300 1.546445e-04 0.011481640           100
566 566  26 14   9         3276 1.533928e-04 0.011481640            33
716 716  10 21  22         4620 2.241546e-04 0.013807925           100
887 887  49  4  23         4508 2.182156e-04 0.013807925            67
36   36  25 39   7         6825 3.419914e-04 0.014363637            67
73   73  17  8  40         5440 2.678014e-04 0.014363637           100
131 131  36  5  36         6480 3.234696e-04 0.014363637            67
325 325  40 35   4         5600 2.763464e-04 0.014363637           100
378 378 109 11   5         5995 2.974741e-04 0.014363637             0
500 500  39  7  25         6825 3.419914e-04 0.014363637            67
551 551   3 10 210         6300 3.138155e-04 0.014363637           100
549 549  32 24  10         7680 3.879766e-04 0.015586537            67
372 372   8 43  29         9976 5.118029e-04 0.019704413           100
137 137   1 57 190        10830 5.579120e-04 0.020620426            67
54   54  37 20  16        11840 6.124434e-04 0.021765296           100
358 358  38 17  27        17442 9.142112e-04 0.031286339           100
223 223  15 34  42        21420 1.127257e-03 0.037199480           100
777 777  20 72  20        28800 1.519046e-03 0.048399936            33
351 351  80 27  14        30240 1.594944e-03 0.049124268            33
> 



length(TNBC.influent) ##sun
30

# brca type of the influent patients


outliers_ensemble=rownames(x_origin924)[TNBC.influent]  ##sun
outliers_ensemble
[1] "TCGA-E9-A1ND-01A-11R-A144-07" "TCGA-AR-A1AJ-01A-21R-A12P-07"
 [3] "TCGA-A2-A04U-01A-11R-A115-07" "TCGA-E9-A22G-01A-11R-A157-07"
 [5] "TCGA-OL-A97C-01A-32R-A41B-07" "TCGA-AC-A62X-01A-11R-A29R-07"
 [7] "TCGA-A7-A13D-01A-13R-A277-07" "TCGA-BH-A42U-01A-12R-A24H-07"
 [9] "TCGA-OL-A5S0-01A-11R-A28M-07" "TCGA-E2-A1II-01A-11R-A144-07"
[11] "TCGA-A2-A0EQ-01A-11R-A034-07" "TCGA-B6-A0IJ-01A-11R-A034-07"
[13] "TCGA-C8-A26Y-01A-11R-A16F-07" "TCGA-A2-A1G6-01A-11R-A13Q-07"
[15] "TCGA-BH-A5IZ-01A-11R-A27Q-07" "TCGA-A2-A0YJ-01A-11R-A109-07"
[17] "TCGA-AR-A1AH-01A-11R-A12D-07" "TCGA-BH-A0DL-01A-11R-A115-07"
[19] "TCGA-E9-A1NC-01A-21R-A26B-07" "TCGA-AO-A03U-01B-21R-A10J-07"
[21] "TCGA-S3-AA0Z-01A-11R-A41B-07" "TCGA-C8-A3M7-01A-12R-A21T-07"
[23] "TCGA-A7-A13E-01A-11R-A12P-07" "TCGA-A2-A3Y0-01A-11R-A239-07"
[25] "TCGA-EW-A1OV-01A-11R-A144-07" "TCGA-LL-A5YP-01A-21R-A28M-07"
[27] "TCGA-LL-A6FR-01A-12R-A31O-07" "TCGA-BH-A1EW-01A-11R-A137-07"
[29] "TCGA-D8-A1JF-01A-11R-A13Q-07" "TCGA-A2-A1G1-01A-21R-A13Q-07"




sum(y[TNBC.influent]) #sun  10 TNBC and 20 non-TNBC patients

## [1] 10

sub_outliers_ensemble<- substr(outliers_ensemble,1,12)  
sub_outliers_ensemble
[1] "TCGA-E9-A1ND" "TCGA-AR-A1AJ" "TCGA-A2-A04U" "TCGA-E9-A22G"
 [5] "TCGA-OL-A97C" "TCGA-AC-A62X" "TCGA-A7-A13D" "TCGA-BH-A42U"
 [9] "TCGA-OL-A5S0" "TCGA-E2-A1II" "TCGA-A2-A0EQ" "TCGA-B6-A0IJ"
[13] "TCGA-C8-A26Y" "TCGA-A2-A1G6" "TCGA-BH-A5IZ" "TCGA-A2-A0YJ"
[17] "TCGA-AR-A1AH" "TCGA-BH-A0DL" "TCGA-E9-A1NC" "TCGA-AO-A03U"
[21] "TCGA-S3-AA0Z" "TCGA-C8-A3M7" "TCGA-A7-A13E" "TCGA-A2-A3Y0"
[25] "TCGA-EW-A1OV" "TCGA-LL-A5YP" "TCGA-LL-A6FR" "TCGA-BH-A1EW"
[29] "TCGA-D8-A1JF" "TCGA-A2-A1G1"



 write.table(sub_outliers_ensemble,file="C:/real data analysis_TNBC20200429/ensemble/sub_outliers_ensemble.txt") 

######

suspect<-read.table("28suspect.txt")
suspect<-as.data.frame(suspect)
suspect28<-as.character(suspect$V1)
##  Coincident with 28 suspected individuals
sub_outliers_ensemble[which(sub_outliers_ensemble %in% suspect28)]

[1] "TCGA-A2-A04U" "TCGA-A2-A0EQ" "TCGA-LL-A5YP"



ER_PR_HER<-read.table(file="ER_PR_HER.txt")
a<-rownames(ER_PR_HER)

ER_PR_HER<-cbind(a,ER_PR_HER)
a=outliers_ensemble
idd<-1:length(a)
        
 a=as.data.frame(cbind(a,idd)) 
   
   ##The individual corresponding to outliers_ensemble in ER_PR_HER was screened by Merge function

   outliers_expression<-merge(a,ER_PR_HER,by="a",all=F)
  outliers_expression<-as.matrix(outliers_expression)
 
    ###按照Outliers_ensemble的排序方式
    outliers_expression[order(as.numeric(outliers_expression[,2])),]
    
  
    a                              idd  ENSG00000091831 ENSG00000082175 ENSG00000141736 y  
 [1,] "TCGA-E9-A1ND-01A-11R-A144-07" "1"  " 1.43733608"   " 0.05212077"   " 13.054248"    "0"
 [2,] "TCGA-AR-A1AJ-01A-21R-A12P-07" "2"  " 1.46653566"   " 0.06680532"   "  9.742704"    "0"
 [3,] "TCGA-A2-A04U-01A-11R-A115-07" "3"  " 0.02393686"   " 0.02283666"   "  9.636197"    "0"
 [4,] "TCGA-E9-A22G-01A-11R-A157-07" "4"  " 0.44347936"   " 0.01983464"   " 15.317028"    "0"
 [5,] "TCGA-OL-A97C-01A-32R-A41B-07" "5"  "16.24660389"   " 8.56199213"   " 24.037033"    "1"
 [6,] "TCGA-AC-A62X-01A-11R-A29R-07" "6"  " 0.18607614"   " 0.02166699"   " 28.531373"    "0"
 [7,] "TCGA-A7-A13D-01A-13R-A277-07" "7"  " 0.52408030"   " 0.80807584"   " 42.281452"    "0"
 [8,] "TCGA-BH-A42U-01A-12R-A24H-07" "8"  " 9.18726081"   " 1.82763551"   " 38.367228"    "1"
 [9,] "TCGA-OL-A5S0-01A-11R-A28M-07" "9"  " 0.09174435"   " 0.06319600"   " 31.919459"    "0"
[10,] "TCGA-E2-A1II-01A-11R-A144-07" "10" " 0.14270408"   " 0.19413278"   " 10.728130"    "0"
[11,] "TCGA-A2-A0EQ-01A-11R-A034-07" "11" " 2.12986768"   " 0.03993924"   " 30.146140"    "1"
[12,] "TCGA-B6-A0IJ-01A-11R-A034-07" "12" " 1.18315709"   " 0.45871168"   " 11.124741"    "0"
[13,] "TCGA-C8-A26Y-01A-11R-A16F-07" "13" " 0.11584087"   " 0.05476926"   " 22.916953"    "1"
[14,] "TCGA-A2-A1G6-01A-11R-A13Q-07" "14" "23.89538282"   "21.45268111"   " 29.740007"    "1"
[15,] "TCGA-BH-A5IZ-01A-11R-A27Q-07" "15" " 5.11519233"   " 0.03021609"   " 28.078545"    "0"
[16,] "TCGA-A2-A0YJ-01A-11R-A109-07" "16" " 0.08534338"   " 0.03114457"   "240.243961"    "0"
[17,] "TCGA-AR-A1AH-01A-11R-A12D-07" "17" " 0.03165734"   " 0.03035718"   " 34.123225"    "0"
[18,] "TCGA-BH-A0DL-01A-11R-A115-07" "18" " 6.98732049"   " 0.04055316"   "  9.921791"    "0"
[19,] "TCGA-E9-A1NC-01A-21R-A26B-07" "19" " 0.10552295"   " 0.07491512"   " 15.905961"    "0"
[20,] "TCGA-AO-A03U-01B-21R-A10J-07" "20" " 0.55896368"   " 0.12362745"   " 17.063622"    "1"
[21,] "TCGA-S3-AA0Z-01A-11R-A41B-07" "21" "16.66754167"   " 0.06623186"   " 33.066854"    "0"
[22,] "TCGA-C8-A3M7-01A-12R-A21T-07" "22" " 4.26609694"   " 0.75780988"   " 25.470280"    "1"
[23,] "TCGA-A7-A13E-01A-11R-A12P-07" "23" " 0.82330363"   " 0.05741937"   " 46.079485"    "0"
[24,] "TCGA-A2-A3Y0-01A-11R-A239-07" "24" " 2.18089099"   " 0.02724167"   " 11.336328"    "0"
[25,] "TCGA-EW-A1OV-01A-11R-A144-07" "25" " 0.22802146"   " 0.02717836"   " 28.908449"    "1"
[26,] "TCGA-LL-A5YP-01A-21R-A28M-07" "26" " 0.15598544"   " 0.05139112"   " 15.096317"    "0"
[27,] "TCGA-LL-A6FR-01A-12R-A31O-07" "27" " 0.32548582"   " 0.03895484"   " 32.131599"    "0"
[28,] "TCGA-BH-A1EW-01A-11R-A137-07" "28" "29.97892747"   "18.90331520"   " 42.465462"    "1"
[29,] "TCGA-D8-A1JF-01A-11R-A13Q-07" "29" " 1.25673293"   " 0.11023371"   " 32.930532"    "1"
[30,] "TCGA-A2-A1G1-01A-21R-A13Q-07" "30" " 0.53237004"   " 0.17102961"   "819.759153"    "0"


   

clinical.tnbc.aux2<-read.table("clinical.tnbc.aux2.txt",header=TRUE)

dim(clinical.tnbc.aux2)

###[1] 1090    6


#####sub_outliers_ensemble


b<-rownames(clinical.tnbc.aux2)


clinical.tnbc.aux2<-cbind(b,clinical.tnbc.aux2)
b=sub_outliers_ensemble
idd<-1:length(b)
        
 b=as.data.frame(cbind(b,idd)) 
   
   ##The individual corresponding to sub_outliers_ensemble in ER_PR_HER was screened by Merge function


   outliers_IHC_FISH<-merge(b,clinical.tnbc.aux2,by="b",all=F)
  outliers_IHC_FISH<-as.matrix(outliers_IHC_FISH)
 
    ###按照Outliers_ensemble的排序方式
    outliers_IHC_FISH[order(as.numeric(outliers_IHC_FISH[,2])),]
    

b              idd  ESR1       PGR        HER2_level HER2_status HER2_FISH  my_HER2_level  
 [1,] "TCGA-E9-A1ND" "1"  "Negative" "Negative" ""         "Positive"  ""         ""             
 [2,] "TCGA-AR-A1AJ" "2"  "Positive" "Negative" ""         "Negative"  ""         ""             
 [3,] "TCGA-A2-A04U" "3"  "Negative" "Negative" "1+"       "Negative"  "Positive" "Negative"     
 [4,] "TCGA-E9-A22G" "4"  "Negative" "Negative" ""         "Positive"  ""         ""             
 [5,] "TCGA-OL-A97C" "5"  "Negative" "Negative" ""         ""          "Negative" ""             
 [6,] "TCGA-AC-A62X" "6"  "Positive" "Negative" ""         ""          ""         ""             
 [7,] "TCGA-A7-A13D" "7"  "Negative" "Positive" "2+"       "Equivocal" "Negative" "Indeterminate"
 [8,] "TCGA-BH-A42U" "8"  "Negative" "Negative" ""         "Negative"  ""         ""             
 [9,] "TCGA-OL-A5S0" "9"  "Positive" "Negative" ""         ""          "Positive" ""             
[10,] "TCGA-E2-A1II" "10" "Negative" "Positive" "1+"       "Negative"  ""         "Negative"     
[11,] "TCGA-A2-A0EQ" "11" "Negative" "Negative" "3+"       "Positive"  "Negative" "Positive"     
[12,] "TCGA-B6-A0IJ" "12" "Positive" "Positive" ""         ""          ""         ""             
[13,] "TCGA-C8-A26Y" "13" "Negative" "Negative" "1+"       "Negative"  ""         "Negative"     
[14,] "TCGA-A2-A1G6" "14" "Negative" "Negative" "1+"       "Negative"  ""         "Negative"     
[15,] "TCGA-BH-A5IZ" "15" "Positive" "Negative" ""         "Negative"  "Negative" ""             
[16,] "TCGA-A2-A0YJ" "16" "Positive" "Negative" "0"        "Negative"  ""         "Negative"     
[17,] "TCGA-AR-A1AH" "17" "Positive" "Negative" ""         "Negative"  ""         ""             
[18,] "TCGA-BH-A0DL" "18" "Positive" "Negative" ""         "Negative"  ""         ""             
[19,] "TCGA-E9-A1NC" "19" "Negative" "Positive" ""         "Positive"  ""         ""             
[20,] "TCGA-AO-A03U" "20" "Negative" "Negative" "0"        "Negative"  "Negative" "Negative"     
[21,] "TCGA-S3-AA0Z" "21" "Positive" "Positive" "1+"       "Equivocal" "Negative" "Negative"     
[22,] "TCGA-C8-A3M7" "22" "Negative" "Negative" ""         "Negative"  ""         ""             
[23,] "TCGA-A7-A13E" "23" "Positive" "Negative" "2+"       "Equivocal" "Negative" "Indeterminate"
[24,] "TCGA-A2-A3Y0" "24" "Positive" "Negative" "1+"       "Negative"  ""         "Negative"     
[25,] "TCGA-EW-A1OV" "25" "Negative" "Negative" ""         "Negative"  "Negative" ""             
[26,] "TCGA-LL-A5YP" "26" "Positive" "Negative" "1+"       "Negative"  "Positive" "Negative"     
[27,] "TCGA-LL-A6FR" "27" "Negative" "Positive" "2+"       "Equivocal" "Positive" "Indeterminate"
[28,] "TCGA-BH-A1EW" "28" "Negative" "Negative" ""         "Negative"  ""         ""             
[29,] "TCGA-D8-A1JF" "29" "Negative" "Negative" "1+"       "Negative"  ""         "Negative"     
[30,] "TCGA-A2-A1G1" "30" "Negative" "Negative" "2+"       "Equivocal" "Positive" "Indeterminate"





colnames(x_origin924)[var.selected1]
  [1] "ENSG00000005187" "ENSG00000007372" "ENSG00000007933" "ENSG00000027869" "ENSG00000064787" "ENSG00000065485"
  [7] "ENSG00000065675" "ENSG00000072041" "ENSG00000073464" "ENSG00000074410" "ENSG00000082196" "ENSG00000082781"
 [13] "ENSG00000083454" "ENSG00000086991" "ENSG00000088386" "ENSG00000090520" "ENSG00000090621" "ENSG00000091831"
 [19] "ENSG00000092758" "ENSG00000099812" "ENSG00000100121" "ENSG00000100312" "ENSG00000100453" "ENSG00000100599"
 [25] "ENSG00000101224" "ENSG00000101425" "ENSG00000101812" "ENSG00000102034" "ENSG00000102387" "ENSG00000102837"
 [31] "ENSG00000103310" "ENSG00000103546" "ENSG00000103740" "ENSG00000105675" "ENSG00000106541" "ENSG00000106819"
 [37] "ENSG00000107485" "ENSG00000107807" "ENSG00000109576" "ENSG00000109861" "ENSG00000109917" "ENSG00000109927"
 [43] "ENSG00000109956" "ENSG00000110693" "ENSG00000111725" "ENSG00000112874" "ENSG00000115361" "ENSG00000115523"
 [49] "ENSG00000115590" "ENSG00000115596" "ENSG00000115648" "ENSG00000118513" "ENSG00000119396" "ENSG00000119403"
 [55] "ENSG00000119614" "ENSG00000120262" "ENSG00000120440" "ENSG00000123545" "ENSG00000123569" "ENSG00000124107"
 [61] "ENSG00000124134" "ENSG00000124243" "ENSG00000126067" "ENSG00000130234" "ENSG00000130943" "ENSG00000131748"
 [67] "ENSG00000132164" "ENSG00000132688" "ENSG00000134202" "ENSG00000134215" "ENSG00000134453" "ENSG00000134538"
 [73] "ENSG00000134802" "ENSG00000134830" "ENSG00000134873" "ENSG00000135898" "ENSG00000135912" "ENSG00000136425"
 [79] "ENSG00000136636" "ENSG00000136842" "ENSG00000136944" "ENSG00000137673" "ENSG00000138138" "ENSG00000139515"
 [85] "ENSG00000140067" "ENSG00000140451" "ENSG00000140464" "ENSG00000140511" "ENSG00000140545" "ENSG00000140905"
 [91] "ENSG00000141448" "ENSG00000141505" "ENSG00000141736" "ENSG00000141741" "ENSG00000141854" "ENSG00000142319"
 [97] "ENSG00000142700" "ENSG00000143412" "ENSG00000143653" "ENSG00000144029" "ENSG00000144119" "ENSG00000144354"
[103] "ENSG00000144681" "ENSG00000144840" "ENSG00000145194" "ENSG00000145220" "ENSG00000145949" "ENSG00000146147"
[109] "ENSG00000146648" "ENSG00000147324" "ENSG00000148459" "ENSG00000148468" "ENSG00000148602" "ENSG00000148677"
[115] "ENSG00000148704" "ENSG00000148965" "ENSG00000149043" "ENSG00000150551" "ENSG00000150656" "ENSG00000150667"
[121] "ENSG00000151287" "ENSG00000151379" "ENSG00000152455" "ENSG00000152457" "ENSG00000152977" "ENSG00000153266"
[127] "ENSG00000153292" "ENSG00000153574" "ENSG00000153774" "ENSG00000154127" "ENSG00000154229" "ENSG00000154548"
[133] "ENSG00000154645" "ENSG00000155868" "ENSG00000158164" "ENSG00000158553" "ENSG00000159753" "ENSG00000160182"
[139] "ENSG00000160282" "ENSG00000160345" "ENSG00000161395" "ENSG00000162755" "ENSG00000163032" "ENSG00000163235"
[145] "ENSG00000163737" "ENSG00000164117" "ENSG00000164306" "ENSG00000164509" "ENSG00000164609" "ENSG00000164825"
[151] "ENSG00000165066" "ENSG00000165125" "ENSG00000165186" "ENSG00000165568" "ENSG00000165685" "ENSG00000165695"
[157] "ENSG00000166573" "ENSG00000166794" "ENSG00000167633" "ENSG00000167703" "ENSG00000167755" "ENSG00000168135"
[163] "ENSG00000168268" "ENSG00000168309" "ENSG00000168515" "ENSG00000168703" "ENSG00000169752" "ENSG00000169994"
[169] "ENSG00000170426" "ENSG00000170516" "ENSG00000170624" "ENSG00000171234" "ENSG00000171604" "ENSG00000171847"
[175] "ENSG00000172350" "ENSG00000172738" "ENSG00000172987" "ENSG00000173467" "ENSG00000173660" "ENSG00000173673"
[181] "ENSG00000173762" "ENSG00000174607" "ENSG00000174963" "ENSG00000175329" "ENSG00000175592" "ENSG00000176076"
[187] "ENSG00000176165" "ENSG00000176428" "ENSG00000176485" "ENSG00000176532" "ENSG00000176540" "ENSG00000176597"
[193] "ENSG00000176842" "ENSG00000180287" "ENSG00000180318" "ENSG00000180644" "ENSG00000180773" "ENSG00000181577"
[199] "ENSG00000181830" "ENSG00000182175" "ENSG00000182185" "ENSG00000182898" "ENSG00000183638" "ENSG00000183831"
[205] "ENSG00000183888" "ENSG00000184254" "ENSG00000184661" "ENSG00000185924" "ENSG00000186442" "ENSG00000186675"
[211] "ENSG00000186854" "ENSG00000186897" "ENSG00000186977" "ENSG00000187323" "ENSG00000187715" "ENSG00000187741"
[217] "ENSG00000188338" "ENSG00000188581" "ENSG00000189013" "ENSG00000189167" "ENSG00000196208" "ENSG00000197191"
[223] "ENSG00000197241" "ENSG00000197943" "ENSG00000198729" "ENSG00000198954" "ENSG00000203760" "ENSG00000204279"
[229] "ENSG00000204385" "ENSG00000204479" "ENSG00000205177" "ENSG00000205268" "ENSG00000205549" "ENSG00000205683"
[235] "ENSG00000206105" "ENSG00000206535" "ENSG00000212864" "ENSG00000212938" "ENSG00000214842" "ENSG00000224107"
[241] "ENSG00000224960" "ENSG00000228623" "ENSG00000241697" "ENSG00000258436" "ENSG00000260007" "ENSG00000262406"
[247] "ENSG00000273047" "ENSG00000278566"


colnames(x_origin924)[var.selected2]
 [1] "ENSG00000074410" "ENSG00000091831" "ENSG00000094755" "ENSG00000102243" "ENSG00000106541" "ENSG00000107485"
 [7] "ENSG00000109436" "ENSG00000110484" "ENSG00000115648" "ENSG00000124664" "ENSG00000124935" "ENSG00000129514"
[13] "ENSG00000151892" "ENSG00000159763" "ENSG00000160180" "ENSG00000160182" "ENSG00000171428" "ENSG00000173467"
[19] "ENSG00000175356" "ENSG00000181617" "ENSG00000186832" "ENSG00000204385"


colnames(x_origin924)[var.selected3]
1] "ENSG00000012223" "ENSG00000074410" "ENSG00000091831" "ENSG00000094755" "ENSG00000101210" "ENSG00000102243"
 [7] "ENSG00000104332" "ENSG00000106541" "ENSG00000106819" "ENSG00000110484" "ENSG00000124107" "ENSG00000124935"
[13] "ENSG00000128422" "ENSG00000129514" "ENSG00000131771" "ENSG00000135374" "ENSG00000137673" "ENSG00000141736"
[19] "ENSG00000141738" "ENSG00000160180" "ENSG00000160182" "ENSG00000160678" "ENSG00000161395" "ENSG00000167754"
[25] "ENSG00000167755" "ENSG00000172551" "ENSG00000173467" "ENSG00000178372" "ENSG00000181617" "ENSG00000186847"
[31] "age.vect"    

 a<-intersect(var.selected1,var.selected2)

 variables_selected<-colnames(x_origin924)[intersect(a,var.selected3)]

 "ENSG00000074410" "ENSG00000091831" "ENSG00000106541" "ENSG00000160182" "ENSG00000173467"


 "ENSG00000074410" :CA12: carbonic anhydrase 12
"ENSG00000091831" ESR1: estrogen receptor 1
 "ENSG00000106541"AGR2: anterior gradient 2, protein disulphide isomerase family member
"ENSG00000160182" TFF1: trefoil factor 1
"ENSG00000173467"  AGR3: anterior gradient 3, protein disulphide isomerase family member

 a12<-intersect(var.selected1,var.selected2)

 setdiff(colnames(x_origin924)[a12],variables_selected)
   "ENSG00000107485" "ENSG00000115648" "ENSG00000204385"

  "ENSG00000107485":GATA3: GATA binding protein 3
  ENSG00000115648:MLPH: melanophilin
ENSG00000204385:SLC44A4: solute carrier family 44 member 4

 a23<-intersect(var.selected3,var.selected2)

  setdiff(colnames(x_origin924)[a23],variables_selected)
   "ENSG00000094755" "ENSG00000102243" "ENSG00000110484" "ENSG00000124935" "ENSG00000129514" "ENSG00000160180" "ENSG00000181617"

ENSG00000094755:GABRP: gamma-aminobutyric acid type A receptor subunit pi
ENSG00000102243:VGLL1: vestigial like family member 1
ENSG00000110484:SCGB2A2: secretoglobin family 2A member 2
ENSG00000124935:SCGB1D2: secretoglobin family 1D member 2
ENSG00000129514:FOXA1: forkhead box A1
ENSG00000160180:TFF3: trefoil factor 3
ENSG00000181617:FDCSP: follicular dendritic cell secreted protein

 a13<-intersect(var.selected3,var.selected1)
   setdiff(colnames(x_origin924)[a13],variables_selected)
 [1] "ENSG00000106819" "ENSG00000124107" "ENSG00000137673" "ENSG00000141736" "ENSG00000161395" "ENSG00000167755"

ENSG00000106819:ASPN: asporin
ENSG00000124107:SLPI: secretory leukocyte peptidase inhibitor
ENSG00000137673:MMP7: matrix metallopeptidase 7
ENSG00000141736:ERBB2: erb-b2 receptor tyrosine kinase 2
ENSG00000161395:PGAP3: post-GPI attachment to proteins phospholipase 3
ENSG00000167755:KLK6: kallikrein related peptidase 6




Marta_outliers<-read.table("Marta_outliers.txt")
Marta_outliers<-as.data.frame(Marta_outliers)
Marta_outliers<-as.character(Marta_outliers$V1)
##  Coincident with 22 outliers detected by Lopes
sub_outliers_ensemble[which(sub_outliers_ensemble %in% Marta_outliers)]

[1] "TCGA-E9-A1ND" "TCGA-AR-A1AJ" "TCGA-E9-A22G" "TCGA-OL-A97C" "TCGA-AC-A62X" "TCGA-BH-A42U" "TCGA-OL-A5S0" "TCGA-E2-A1II"
 [9] "TCGA-A2-A0EQ" "TCGA-A2-A1G6" "TCGA-A2-A0YJ" "TCGA-AR-A1AH" "TCGA-C8-A3M7" "TCGA-A2-A3Y0" "TCGA-LL-A5YP" "TCGA-D8-A1JF"


 