### Load relevant libraries ###
library(foreign)
library(Hmisc)
library(plyr)
library(mice)
library(VIM)
library(lattice)
library(ggplot2)

### Set working directory ###
setwd("final\\0. Imputation") 

### Load the data in ###
data <- read.csv(file="imputation_totalModel.csv")
data <- data.frame(data)

### Examine missing data patterns ###
md.pattern(data)

### Factorize categorical variables ###
# get the number of unique values for each variable to determine which variables to factorize (for me: categorical vars are < 12 levels)
count <- rep(NA,dim(data)[2])
for (i in 1:dim(data)[2]) {
  count[i] <- length(table(data[,i]))
  names(count) <- colnames(data)
}

factorize = function (x, columns = names(data)) {   # for better imputation -> categorical variables +1 (so 1 and 2 instead of 0 and 1)
  data[columns] = lapply(data[columns], function(x) as.factor(x));data}

data <- factorize(data, which(count<12))

### Imputation ###
set.seed(65486)  
data.imp <- aregImpute(as.formula(paste("~", paste(names(data), collapse='+'))),
                       n.impute=1, data=data, x=TRUE, nk = 0)  # use all variables in data frame for imputation
data.imp$nna # number of missing data per variable
data<-data.frame(data.imp$x)  # imputed data set

### Make factors numerical again ###
num = function (x, columns = names(data)) {  
  data[columns] = lapply(data[columns], function(x) as.numeric(x));data}

data <- num(data, which(count<12))

# reverse earlier conversion (i.e., now return to 0 and 1 instead of 1 and 2)
min1 =   function (x, columns = names(data)) { 
  data[columns] = lapply(data[columns], function(x){x-1}); data 
} 
data <- min1(data,which(count<12)); rm(min1);rm(num)

### Save the data! ###
write.csv(data, "imputed_totalModel.csv", row.names=FALSE)








