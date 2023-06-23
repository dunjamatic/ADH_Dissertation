# Clear environment
rm(list=ls())

# Load libraries
library(foreign);library(Hmisc); library(car); library(rms); library(cmprsk);
library(survival)
library(ggplot2)
library(plotly)
library(prodlim)
library(riskRegression)
library(pec)
library(ggsurvfit)

# Set working directory
setwd("final\\1. Objective 1") 

# Load data
data<-read.csv("imputed_totalModel.csv")

# re-name height, weight, crp columns
colnames(data)[colnames(data) == "height..cm."] ="height"
colnames(data)[colnames(data) == "weight..kg."] ="weight"
colnames(data)[colnames(data) == "crp.value..mg.L."] ="crp"

# Remove variables that are not needed like ethnicity and townsend 
myvars<-c("e_patid", "baseline_age", "followup_age", "followup_time", "gender_encoded",
          "height", "weight", "alcohol_status", "smoking_status", "diabetes_status",
          "antiplatelets_use", "crp", "total_outcomes")
data<-data[myvars]
rm(myvars)

################################ INTERACTIONS ##################################
# Establish the interactions

# establish interaction with smoking
data$smoking.1.age<-ifelse(data$smoking_status==1,data$followup_time,0) 
data$smoking.2.age<-ifelse(data$smoking_status==2,data$followup_time,0)
data$age.cent<-data$followup_time-mean(data$followup_time)
data$smoking.1.age<-ifelse(data$smoking_status==1,data$age.cent,0)
data$smoking.2.age<-ifelse(data$smoking_status==2,data$age.cent,0)

# establish interactions with CRP
data$logcrp.age<-(log(data$crp))*data$followup_time
data$crp.age<-data$crp*data$followup_time
data$crp.2.age<-(data$crp^2)*data$followup_time

# establish interactions with weight
data$weight.age<-data$weight*data$followup_time
data$weight.sq.age<-(data$weight^2)*data$followup_time
################################################################################


############################## PREP FOR FINE-GRAY ##############################
# useFine-Gray to check PH, transformations, get shrinkage
# Create competing risk endpoint
#  0: Alive or LFU (censored)
#  1: total cancer
#  2: Death by other cause

# Specify units
units(data$followup_time) <- "Year" 

data$comp.event <- data$total_outcomes
describe(data$comp.event)
data$comp.event2 <- factor(data$comp.event, 0:2, labels=c("censor", "cancer", "death")) # set competing events, 0 - censor, 1 - cancer, 2 - death  
data$follow.up<-data$followup_time
data$follow.up.age<-data$followup_age

# Make factors of all categorical data
factorize <-  function (x, columns = names(data), izer=as.factor) { data[columns] = lapply(data[columns], izer);   data } 
data<-factorize(data,c('gender_encoded','diabetes_status','smoking_status','alcohol_status', 'antiplatelets_use'))  
str(data); rm(factorize)
data$gender_encoded <- as.factor(data$gender_encoded)
data$alcohol_status <- as.factor(data$alcohol_status)
data$smoking_status <- as.factor(data$smoking_status)
data$diabetes_status <- as.factor(data$diabetes_status)
data$antiplatelets_use <- as.factor(data$antiplatelets_use)
################################################################################


############################### TRANSFORMATIONS ################################
# Variables in the model:
# Continuous: weight, height, CRP (age is the time axis)
# Categorical: smoking, diabetes, alcohol_status, antiplatelets, gender

# transformations help account for non-linearity of continuous variables
# non-linearity is checked with spline plots and p-value of spline plot which shows whether the spline is necessary for the model
# transformation options include linear, quadr, log, and spline 
# note that lower AIC is better

# Create variable to store AIC which we will use to determine the best transformation
AIC <- matrix(nrow=8,ncol=4)
colnames(AIC) <- c("linear","log","quadr", "spline")
rownames(AIC) <- c("weight1","weight2",
                   "height1","height2",
                   "crp1","crp2",
                   "age1", "age2")

# CANCER MODEL
# Create weighted datasets for cancer (timescale is age at entry until age at exit)
dataweight1 <- finegray(Surv(followup_time,comp.event2) ~ .,id=e_patid, data=data, etype="cancer")
dd <- datadist(dataweight1);options(datadist="dd")

# WEIGHT
# try linear
mod1 <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+weight+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=dataweight1,weight=fgwt,x=T,y=T)
# try log
mod2 <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+log(weight+1)+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=dataweight1,weight=fgwt,x=T,y=T)
# try quadr
mod3 <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+pol(weight,2)+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=dataweight1,weight=fgwt,x=T,y=T)
# try spline
spline <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+rcs(weight,4)+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=dataweight1,weight=fgwt,x=T,y=T)

spline <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+rcs(weight,4)+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=dataweight1,weight=fgwt,x=T,y=T)

# store AIC
AIC[1,1] <- extractAIC(mod1)[2]
AIC[1,2] <- extractAIC(mod2)[2]
AIC[1,3] <- extractAIC(mod3)[2]
AIC[1,4] <- AIC(spline)

plot(Predict(spline,weight, fun=exp), conf.int=T, main="Spline of Weight", xlab="Weight (kg)", ylab="Hazard ratio")

# anova on spline with significance in the spline variable indicates that there is
# a non-linear relationship and provides evidence for why the spline term should be included.
# especially given that F-G has a linearity assumption. 
anova(spline)
# Wald Statistics          Response: Surv(fgstart, fgstop, fgstatus) 
# Factor            Chi-Square d.f. P     
# baseline_age      1412.98     1   <.0001
# weight             396.44     3   <.0001
# Nonlinear         130.20     2   <.0001
# gender_encoded    1259.09     1   <.0001
# height             430.35     1   <.0001
# crp                808.29     1   <.0001
# smoking_status     368.85     2   <.0001
# diabetes_status    421.84     1   <.0001
# alcohol_status     129.55     2   <.0001
# antiplatelets_use  326.21     1   <.0001
# TOTAL             9351.47    13   <.0001


# HEIGHT
# try linear
mod1 <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+height+gender_encoded+weight+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=dataweight1,weight=fgwt,x=T,y=T)
# try log
mod2 <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+log(height+1)+gender_encoded+weight+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=dataweight1,weight=fgwt,x=T,y=T)
# try quadr
mod3 <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+pol(height,2)+gender_encoded+weight+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=dataweight1,weight=fgwt,x=T,y=T)
# try spline
spline <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+rcs(height,4)+gender_encoded+weight+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=dataweight1,weight=fgwt,x=T,y=T)

# store AIC
AIC[3,1] <- extractAIC(mod1)[2]
AIC[3,2] <- extractAIC(mod2)[2]
AIC[3,3] <- extractAIC(mod3)[2]
AIC[3,4] <- AIC(spline)

plot(Predict(spline,height, fun=exp), conf.int=T, main="Spline of Height", xlab="Height (cm)", ylab="Hazard ratio")

# anova on spline with significance in the spline variable indicates that there is
# a non-linear relationship and provides evidence for why the spline term should be included.
# especially given that F-G has a linearity assumption. 
anova(spline) #
# Wald Statistics          Response: Surv(fgstart, fgstop, fgstatus) 
# 
# Factor            Chi-Square d.f. P     
# baseline_age      1410.99     1   <.0001
# height             606.67     3   <.0001
# Nonlinear         125.24     2   <.0001
# gender_encoded     837.59     1   <.0001
# weight             246.47     1   <.0001
# crp                809.72     1   <.0001
# smoking_status     371.11     2   <.0001
# diabetes_status    416.06     1   <.0001
# alcohol_status     129.83     2   <.0001
# antiplatelets_use  331.69     1   <.0001
# TOTAL             9381.74    13   <.0001


# CRP
# try linear
mod1 <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+crp+gender_encoded+height+weight+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=dataweight1,weight=fgwt,x=T,y=T)
# try log
mod2 <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+log(crp+1)+gender_encoded+height+weight+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=dataweight1,weight=fgwt,x=T,y=T)
# try quadr
mod3 <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+pol(crp,2)+gender_encoded+height+weight+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=dataweight1,weight=fgwt,x=T,y=T)
# try spline
spline <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+rcs(crp,4)+gender_encoded+height+weight+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=dataweight1,weight=fgwt,x=T,y=T)

# store AIC
AIC[5,1] <- extractAIC(mod1)[2]
AIC[5,2] <- extractAIC(mod2)[2]
AIC[5,3] <- extractAIC(mod3)[2]
AIC[5,4] <- AIC(spline)

plot(Predict(spline,crp, fun=exp), conf.int=T, main="Spline of CRP", xlab="CRP (mg/L)", ylab="Hazard ratio")

# anova on spline with significance in the spline variable indicates that there is
# a non-linear relationship and provides evidence for why the spline term should be included.
# especially given that F-G has a linearity assumption. 
anova(spline) # 
# Wald Statistics          Response: Surv(fgstart, fgstop, fgstatus) 
# 
# Factor            Chi-Square d.f. P     
# baseline_age       1546.81    1   <.0001
# crp                1949.79    3   <.0001
# Nonlinear         1078.92    2   <.0001
# gender_encoded     1392.99    1   <.0001
# height              441.98    1   <.0001
# weight              290.64    1   <.0001
# smoking_status      373.94    2   <.0001
# diabetes_status     399.75    1   <.0001
# alcohol_status      132.04    2   <.0001
# antiplatelets_use   322.66    1   <.0001
# TOTAL             10383.30   13   <.0001

# AGE
# try linear
mod1 <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+crp+gender_encoded+height+weight+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=dataweight1,weight=fgwt,x=T,y=T)
# try log
mod2 <- cph(Surv(fgstart,fgstop,fgstatus)~log(baseline_age+1)+crp+gender_encoded+height+weight+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=dataweight1,weight=fgwt,x=T,y=T)
# try quadr
mod3 <- cph(Surv(fgstart,fgstop,fgstatus)~pol(baseline_age,2)+crp+gender_encoded+height+weight+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=dataweight1,weight=fgwt,x=T,y=T)
# try spline
spline <- cph(Surv(fgstart,fgstop,fgstatus)~rcs(baseline_age,4)+crp+gender_encoded+height+weight+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=dataweight1,weight=fgwt,x=T,y=T)

# store AIC
AIC[7,1] <- extractAIC(mod1)[2]
AIC[7,2] <- extractAIC(mod2)[2]
AIC[7,3] <- extractAIC(mod3)[2]
AIC[7,4] <- AIC(spline)

plot(Predict(spline,baseline_age, fun=exp), conf.int=T, main="Spline of Age", xlab="Age (years)", ylab="Hazard ratio")

# anova on spline with significance in the spline variable indicates that there is
# a non-linear relationship and provides evidence for why the spline term should be included.
# especially given that F-G has a linearity assumption. 
anova(spline)
# Factor            Chi-Square d.f. P     
# baseline_age       3062.79    3   <.0001
# Nonlinear         1777.21    2   <.0001
# crp                 797.57    1   <.0001
# gender_encoded     1285.49    1   <.0001
# height              474.97    1   <.0001
# weight              272.28    1   <.0001
# smoking_status      355.84    2   <.0001
# diabetes_status     455.79    1   <.0001
# alcohol_status      126.79    2   <.0001
# antiplatelets_use   393.04    1   <.0001
# TOTAL             11009.22   13   <.0001

# lower AIC is better so pick transformation with lowest AIC
print(AIC)
#         linear     log   quadr  spline
# weight1 4340600 4340583 4340489 4340462 quadr -- spline is lowest                         
# height1 4340600 4340624 4340569 4340505 quadr -- spline is lowest 
# crp1    4340600 4339444 4340051 4339210 log   -- spline is lowest 
                        

rm(AIC)
################################################################################


################################# CHECK PH #####################################
dataweight1 <- finegray(Surv(followup_time,comp.event2) ~ .,id=e_patid, data=data, etype="cancer")
dd <- datadist(dataweight1);options(datadist="dd")

# Create model with optimal transformations to check PH assumption on 
model1 <- coxph(Surv(fgstart,fgstop,fgstatus)~baseline_age+gender_encoded+smoking_status+weight+pol(crp,2)+            
                  diabetes_status+alcohol_status+antiplatelets_use+height, data=dataweight1,weight=fgwt,x=T,y=T)
cox.zph(model1) 

#                   chisq  df      p
# baseline_age      32947.4  1 < 2e-16
# gender_encoded      668.2  1 < 2e-16
# smoking_status       50.9  2 8.7e-12
# weight              192.1  1 < 2e-16
# pol(crp, 2)        1514.9  2 < 2e-16
# diabetes_status    2777.4  1 < 2e-16
# alcohol_status     1281.8  2 < 2e-16
# antiplatelets_use  7998.8  1 < 2e-16
# height              566.0  1 < 2e-16
# GLOBAL            37912.6 12 < 2e-16

# SMOKING
AIC(model1) # 4339942
model1b <- coxph(Surv(fgstart,fgstop,fgstatus)~baseline_age+gender_encoded+smoking_status+pol(weight,2)+pol(crp,2)+smoking.1.age+smoking.2.age+
                   height+diabetes_status+alcohol_status+antiplatelets_use, data=dataweight1,weight=fgwt,x=T,y=T) # add interactions w smoking
AIC(model1b) # 4335239 --> AIC is lower when we add smoking interactions 

# CRP 
model1c <- coxph(Surv(fgstart,fgstop,fgstatus)~baseline_age+gender_encoded+smoking_status+pol(weight,2)+crp+smoking.1.age+smoking.2.age+
                   height+diabetes_status+alcohol_status+antiplatelets_use, data=dataweight1,weight=fgwt,x=T,y=T) # try crp instead of crp^2
AIC(model1c) # 4335792 --> need quadratic crp (compare with AIC from model1) 

# WEIGHT
model1d <- coxph(Surv(fgstart,fgstop,fgstatus)~baseline_age+gender_encoded+smoking_status+weight+pol(crp,2)+smoking.1.age+smoking.2.age+
                   height+diabetes_status+alcohol_status+antiplatelets_use, data=dataweight1,weight=fgwt,x=T,y=T) # try weight instead og weight^2
AIC(model1d) # 4335328 --> weight is better than weight^2 (compare with AIC from model1) 
cox.zph(model1d) # significance is still there... not sure what that means 

model1e <- coxph(Surv(fgstart,fgstop,fgstatus)~baseline_age+gender_encoded+smoking_status+weight+pol(crp,2)+smoking.1.age+smoking.2.age+crp.age+crp.2.age+
                   diabetes_status+alcohol_status+height+antiplatelets_use, data=dataweight1,weight=fgwt,x=T,y=T)
AIC(model1e) # 4334368 try adding interactions between crp and age! 
# lower AIC when interactions are added

# Check Cancer model's PH with plots 
somePDFPath = "proportionalHazards_cancerModel.pdf"
pdf(file=somePDFPath)  
plot(cox.zph(model1), resid = F, se = T, lwd = 3)
dev.off() 
rm(somePDFPath)
################################################################################


# to be able to get shrinkage factor!
############################# FINAL CANCER MODEL ###############################
# Create new weighted dataset
dataweight1 <- finegray(Surv(followup_time,comp.event2) ~ .,id=e_patid, data=data, etype="cancer")  
dd <- datadist(dataweight1); options(datadist="dd")

# Model 1 (cancer): w appropriate transformations/interactions/etc.
model1 <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+gender_encoded+smoking_status+weight+pol(crp,2)+smoking.1.age+smoking.2.age+
                height+diabetes_status+alcohol_status+antiplatelets_use+crp.age+crp.2.age,data=dataweight1,weight=fgwt,x=T,y=T)
coef1 <- model1$coefficients;coef1

saveRDS(model1, "cancerModel.rds")
#model1 <- readRDS("cancerModel.rds")

########################## SHRINKAGE FOR CANCER MODEL ##########################
n.boot <- 50
b.coefs <- NULL
set.seed(3)

for (i in 1:n.boot) {
  print(i)
  b.sample <- sample(1:dim(dataweight1)[1],dim(dataweight1)[1],replace=TRUE)
  b.data <- dataweight1[b.sample,]
  b.model1 <- cph(Surv(fgstart,fgstop,fgstatus)~gender_encoded+smoking_status+weight+pol(crp,2)+smoking.1.age+smoking.2.age+crp.age+crp.2.age+
                    height+diabetes_status+alcohol_status+antiplatelets_use,data=b.data,weight=fgwt,x=T,y=T)
  dataweight1$b.lp <- predict(b.model1, newdata = dataweight1, type="lp")
  b.coefs[i] <- cph(Surv(fgstart,fgstop,fgstatus)~b.lp, data=dataweight1, weight=fgwt, x=T, y=T)$coefficients[1]
}

shrink <- mean(b.coefs); shrink 
model1.shrink <- model1
model1.shrink$coefficients <- model1.shrink$coefficients*shrink
coef1 <- model1.shrink$coefficients # to save to use for validations and predictions

saveRDS(model1.shrink, "cancerModel_shrunk.rds")

#################################### FGR WAY ###################################
############################# FINAL CANCER MODEL ###############################
data$crp_quadr <- pol(data$crp, 2)

# using follow up time as the underlying timescale because FGR does not support left truncation (as per its documentation)
finegray_final <- FGR(Hist(followup_time,total_outcomes)~baseline_age+gender_encoded+ 
                                 smoking_status+smoking.1.age+smoking.2.age+
                                 crp_quadr+crp.age+crp.2.age+
                                 height+
                                 diabetes_status+
                                 alcohol_status+antiplatelets_use+
                                 weight,data=data, cause = 1, maxiter = 500)

################################## SHRINKAGE ###################################
# Shrinkage --> shrinkage factor obtained from the finegray() method
shrink <- 0.9901548
finegray_final.shrink <- finegray_final
finegray_final.shrink$coefficients <- finegray_final.shrink$coefficients*shrink
coef1 <- finegray_final.shrink$crrFit$coef # to save to use for validations and predictions

saveRDS(finegray_final.shrink, "finegray_final_shrink.rds") 

# get hazard ratios and coefs of the final model (shrunk model)
coef2<-round((coef1),6)
HR2<-round(cbind(exp(finegray_final.shrink$crrFit$coef)),6)
table2 <- cbind (coef1,HR2); colnames(table2) <- c("coef","HR")
table2

units(data$followup_time) <- "Year" 

################################# VALIDATION ###################################
# Make the predictions
risk <- predictRisk(finegray_final.shrink,times=10, newdata=data)

# define relevant variables 
data$plotrisk <- risk
n.groups <- 10
year.risk <- 10

# C-statistic
# surv.obj <- with(data, Surv(followup_age-baseline_age,total_outcomes)) # type="mstate" is used for competing risks but requires expansion to as.numeric which results in non-equal length of vars in rcorr
surv.obj <- with(data, Surv(followup_time,total_outcomes))
rcorr1 <- rcorr.cens(x=(1-data$plotrisk)*10,S=surv.obj)
se.1 <- rcorr1["S.D."]/2
Low95.1 <- rcorr1["C Index"] - 1.96*se.1
Upper95.1 <- rcorr1["C Index"] + 1.96*se.1
c.stat<-cbind(rcorr1["C Index"], Low95.1, Upper95.1); c.stat

# get observed and predicted risks in deciles
data$groups1 <- as.numeric(cut2(data$plotrisk, g=n.groups, levels.mean=TRUE)) # divide data into equal number of groups based on predicted CV risk
time = data$followup_age - data$baseline_age
CIF1 <- cuminc(ftime=time, fstatus=data$total_outcomes, cencode="0", group=data$groups1)
error1 <- as.numeric(timepoints(CIF1[c(1:n.groups)], year.risk)$var)  # error = variance of cum. inc function of the number of years at risk per group
observeda1 <- as.numeric(timepoints(CIF1[c(1:n.groups)], year.risk)$est)  # observed = point estimate of cum. inc function of the number of years at risk per group
predicteda1<- tapply(data$plotrisk, data$groups1, mean)  # average of the predicted risk per group

# Calibration plot
plot(predicteda1, observeda1, xlab=" ", type = "b", ylab=" ", pch = 15, xaxt = "n",yaxt = "n",xlim=c(0,0.6),ylim=c(0,0.6), cex.axis=1, cex.lab=1)
title("Calibration: Fine-Gray Model", cex.main=1.2,line=1)
errbar(predicteda1,observeda1, (exp(log(observeda1+1.96*sqrt(error1)))), (exp(log(observeda1-1.96*sqrt(error1)))),add=T, lwd = 0.5)
abline(a=0,b=1, lty=4)
lines(predicteda1,observeda1)
mtext("Predicted 10-year risk",cex=1.2,side=1,adj=0.5,padj=+3.5)
mtext("Observed 10-year risk",cex=1.2,side=2,adj=0.5,padj=-3.5)
axis(1, at=c(0:10/10), labels=paste(c(0:10*10),"%",sep=""), lwd = 1, cex.axis=1.2)
axis(2, at=c(0:10/10), labels=paste(c(0:10*10),"%",sep=""), lwd = 1, cex.axis=1.2)
legend("topleft", paste("C-statistic: 0.671 (95%CI 0.670 - 0.672)"))


############################## OTHER CALCULATIONS ##############################
# Mean predicted risk and range 
risk_percent <- risk*100
mean(risk_percent) 
min(risk_percent) 
max(risk_percent) 

# Distribution of Predicted risk
h = hist(risk) # or hist(x,plot=FALSE) to avoid the plot of the histogram
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, ylim=c(0,50), main = "Distribution of Predicted Risk: Fine-Gray Model", xlab = "Predicted Risk", ylab = "Percent")
################################################################################










