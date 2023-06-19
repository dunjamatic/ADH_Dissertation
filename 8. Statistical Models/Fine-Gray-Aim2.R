### Fine Gray ###

# Clear environment
rm(list=ls())

# Load libraries
library(survival)
library(rms)
library(Hmisc)
library(foreign)
library(car)
library(cmprsk)
library(prodlim)
library(riskRegression)
library(pec)
library(plotly)
library(ggplot2)
library(ggsurvfit)

# Set working directory
setwd("final/2a. Compare FG to ML") 

# Load data
train_fg_data <- read.csv("FGTraining_totalModel.csv")
test_fg_data <- read.csv("FGTesting_totalModel.csv")

# re-name height, weight, crp columns
colnames(train_fg_data)[colnames(train_fg_data) == "height..cm."] ="height"
colnames(train_fg_data)[colnames(train_fg_data) == "weight..kg."] ="weight"
colnames(train_fg_data)[colnames(train_fg_data) == "crp.value..mg.L."] ="crp"
colnames(test_fg_data)[colnames(test_fg_data) == "height..cm."] ="height"
colnames(test_fg_data)[colnames(test_fg_data) == "weight..kg."] ="weight"
colnames(test_fg_data)[colnames(test_fg_data) == "crp.value..mg.L."] ="crp"


################################# PREP FGR #####################################
# Specify units
units(train_fg_data$followup_time) <- "Year" 

# Factorize categorical variables so the model knows they are categorical 
# for training data
factorize <-  function (x, columns = names(train_fg_data), izer=as.factor) { train_fg_data[columns] = lapply(train_fg_data[columns], izer);   train_fg_data } 
train_fg_data<-factorize(train_fg_data,c('gender_encoded','alcohol_status','smoking_status','diabetes_status','antiplatelets_use')) 
str(train_fg_data); rm(factorize)
################################################################################


################################  CHECK PH  ####################################
# Prepare data for finegray() function
train_fg_data$comp.event <- train_fg_data$total_outcomes
train_fg_data$comp.event2 <- factor(train_fg_data$comp.event, 0:2, labels=c("censor", "cancer", "death")) # set competing events, 0 - censor, 1 - cancer, 2 - death

# Create Fine-Gray model
dataweight1 <- finegray(Surv(followup_time,comp.event2) ~ .,id=e_patid, data=train_fg_data, etype="cancer")
dd <- datadist(dataweight1);options(datadist="dd")

f_ph <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+weight+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=dataweight1,weight=fgwt,x=T,y=T)

# Now check PH assumption 
# first with plots
somePDFPath = "proportionalHazards_FineGray_trainModel.pdf"
pdf(file=somePDFPath)  
plot(cox.zph(f_ph), resid = F, se = T, lwd = 3)
dev.off() 
rm(somePDFPath)
# but also with p-value (p<0.05 indicates violation of PH assumption)
cox.zph(f_ph)

# cox.zph calculates Schoenfeld residuals which are the difference between the observed
# and expected covariate at different time points. a test between the residuals and time
# is then done. it computes rho, which measures the association between the covariate and time.
# null hypothesis of the test is that the correlation coefficient is zero, indicating
# no association between the covariate and time. p < 0.05 suggests violation of PH.

# plot of cox.zph lets you visually assess too. if there are plots that are not a horizontal line
# it might suggest violations of the PH assumption.

#                     chisq  df     p
# baseline_age      6.22e+02  1 < 2e-16 # add interaction
# weight            2.26e+01  1 2.0e-06 # add interaction
# gender_encoded    8.34e+00  1  0.0039 # add interaction
# height            6.56e-02  1  0.7978
# crp               3.23e+01  1 1.3e-08 # add interaction
# smoking_status    2.36e+01  2 7.4e-06 # add interaction
# diabetes_status   1.98e+02  1 < 2e-16 # add interaction
# alcohol_status    2.17e+02  2 < 2e-16 # add interaction
# antiplatelets_use 6.73e+02  1 < 2e-16 # add interaction
# GLOBAL            1.40e+03 11 < 2e-16 # add interaction
################################################################################


##############################  INTERACTIONS  ##################################
# interactions created are based on what the PH assumption results showed 
# establish interactions with age
train_fg_data$age.time<-train_fg_data$baseline_age*train_fg_data$followup_time

# establish interactions with weight
train_fg_data$weight.time<-train_fg_data$weight*train_fg_data$followup_time

# establish interaction with gender
train_fg_data$gender.1.time<-ifelse(train_fg_data$gender_encoded==1,train_fg_data$followup_time,0)

# establish interactions with CRP
train_fg_data$crp.time<-train_fg_data$crp*train_fg_data$baseline_age

# establish interaction with smoking
train_fg_data$smoking.1.time<-ifelse(train_fg_data$smoking_status==1,train_fg_data$followup_time,0)
train_fg_data$smoking.2.time<-ifelse(train_fg_data$smoking_status==2,train_fg_data$followup_time,0)

# establish interaction with diabetes
train_fg_data$diabetes_status.time<-ifelse(train_fg_data$diabetes_status==1,train_fg_data$followup_time,0)

# establish interaction with alcohol
train_fg_data$alcohol.1.time<-ifelse(train_fg_data$alcohol_status==1,train_fg_data$followup_time,0)
train_fg_data$alcohol.2.time<-ifelse(train_fg_data$alcohol_status==2,train_fg_data$followup_time,0)

# establish interaction with antiplatelets
train_fg_data$antiplatelets.time<-ifelse(train_fg_data$antiplatelets_use==1,train_fg_data$followup_time,0)
################################################################################


############################  TRANSFORMATIONS  #################################
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

# WEIGHT
# try linear 
mod1 <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+weight+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = dataweight1, x = T, y = T)
# try log
mod2 <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+log(weight+1)+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = dataweight1, x = T, y = T)
# try quadratic
mod3 <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+pol(weight,2)+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = dataweight1, x = T, y = T)
# store AIC
AIC[1,1] <- extractAIC(mod1)[2]
AIC[1,2] <- extractAIC(mod2)[2]
AIC[1,3] <- extractAIC(mod3)[2]

# try spline
spline <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+rcs(weight,4)+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
              data=dataweight1,x=T,y=T)
AIC[1,4] <- AIC(spline) 

# Plot spline because it will show linearity (or lack thereof) 
# plot shows relationships is non-linear for weight
# p value < 0.05 shows that the split contributes to the model significantly 
# anova on spline with significance in the spline variable indicates that there is
# a non-linear relationship and provides evidence for why the spline term should be included.
# especially given that F-G has a linearity assumption
plot(Predict(spline,weight, fun=exp), conf.int=T, main="Spline of Weight", xlab="Weight (kg)", ylab="Hazard ratio")

# Get anova of the spline too because the p-value indicates whether the spline is important to the model
# p value < 0.05 shows that the split contributes to the model significantly 
# anova on spline with significance in the spline variable indicates that there is
# a non-linear relationship and provides evidence for why the spline term should be included.
# especially given that F-G has a linearity assumption
anova(spline) 
# Wald Statistics          Response: Surv(followup_time, total_outcomes) 
# 
# Factor            Chi-Square d.f. P     
# baseline_age      12982.76    1   <.0001
# weight              175.75    3   <.0001 --> P<0.05 MEANS SPLINE IMPORTANT
# Nonlinear           36.09    2   <.0001
# gender_encoded      479.29    1   <.0001
# height              187.76    1   <.0001
# crp                  25.10    1   <.0001
# smoking_status       59.51    2   <.0001
# diabetes_status       0.41    1   0.5219
# alcohol_status      197.92    2   <.0001
# antiplatelets_use   555.80    1   <.0001
# TOTAL             16491.33   13   <.0001

# HEIGHT
# try linear
mod1 <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+weight+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = dataweight1, x = T, y = T)
# try log
mod2 <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+weight+gender_encoded+log(height+1)+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = dataweight1, x = T, y = T)
# try quadratic
mod3 <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+weight+gender_encoded+pol(height,2)+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = dataweight1, x = T, y = T)
# store AIC
AIC[3,1] <- extractAIC(mod1)[2]
AIC[3,2] <- extractAIC(mod2)[2]
AIC[3,3] <- extractAIC(mod3)[2]

# try spline
spline <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+rcs(height,4)+gender_encoded+weight+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=dataweight1,x=T,y=T)
AIC[3,4] <- AIC(spline)  

# Plot spline because it will show linearity (or lack thereof) 
plot(Predict(spline,height, fun=exp), conf.int=T, main="Spline of Height", xlab="Height (cm)", ylab="Hazard ratio")

# Get anova of the spline too because the p-value indicates whether the spline is important to the model
anova(spline) # 
# Wald Statistics          Response: Surv(followup_time, total_outcomes) 
# 
# Factor            Chi-Square d.f. P     
# baseline_age      12989.71    1   <.0001
# height              257.66    3   <.0001 --> P<0.05 MEANS SPLINE IMPORTANT
# Nonlinear           51.30    2   <.0001
# gender_encoded      320.10    1   <.0001
# weight              131.28    1   <.0001
# crp                  25.15    1   <.0001
# smoking_status       59.79    2   <.0001
# diabetes_status       0.27    1   0.6027
# alcohol_status      197.81    2   <.0001
# antiplatelets_use   561.01    1   <.0001
# TOTAL             16504.80   13   <.0001


# CRP
# try linear
mod1 <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+weight+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = dataweight1, x = T, y = T)
# try log
mod2 <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+weight+gender_encoded+height+log(crp+1)+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = dataweight1, x = T, y = T)
# try quadratic
mod3 <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+weight+gender_encoded+height+pol(crp,2)+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = dataweight1, x = T, y = T)
# store AIC
AIC[5,1] <- extractAIC(mod1)[2]
AIC[5,2] <- extractAIC(mod2)[2]
AIC[5,3] <- extractAIC(mod3)[2]

# try spline
spline <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+rcs(crp,4)+gender_encoded+height+weight+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=dataweight1,x=T,y=T)
AIC[5,4] <- AIC(spline) 

# Plot spline because it will show linearity (or lack thereof) 
plot(Predict(spline,crp, fun=exp), conf.int=T, main="Spline of CRP", xlab="CRP (mg/L)", ylab="Hazard ratio")

# Get anova of the spline too because the p-value indicates whether the spline is important to the model
anova(spline) 
# Wald Statistics          Response: Surv(followup_time, total_outcomes) 
# 
# Factor            Chi-Square d.f. P     
# baseline_age      12959.12    1   <.0001
# crp                  74.50    3   <.0001 --> P<0.05 MEANS SPLINE IMPORTANT
# Nonlinear           48.71    2   <.0001
# gender_encoded      509.52    1   <.0001
# height              203.25    1   <.0001
# weight              142.36    1   <.0001
# smoking_status       59.93    2   <.0001
# diabetes_status       0.38    1   0.5377
# alcohol_status      199.92    2   <.0001
# antiplatelets_use   558.52    1   <.0001
# TOTAL             16495.25   13   <.0001


# AGE
# try linear
mod1 <- cph(Surv(fgstart,fgstop,fgstatus)~baseline_age+weight+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = dataweight1, x = T, y = T)
# try log
mod2 <- cph(Surv(fgstart,fgstop,fgstatus)~log(baseline_age+1)+weight+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = dataweight1, x = T, y = T)
# try quadratic
mod3 <- cph(Surv(fgstart,fgstop,fgstatus)~pol(baseline_age,2)+weight+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = dataweight1, x = T, y = T)
# store AIC
AIC[7,1] <- extractAIC(mod1)[2]
AIC[7,2] <- extractAIC(mod2)[2]
AIC[7,3] <- extractAIC(mod3)[2]

# try spline
spline <- cph(Surv(fgstart,fgstop,fgstatus)~rcs(baseline_age,4)+crp+gender_encoded+height+weight+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=dataweight1,x=T,y=T)
AIC[7,4] <- AIC(spline) 

# Plot spline because it will show linearity (or lack thereof) 
plot(Predict(spline,baseline_age, fun=exp), conf.int=T, main="Spline of Age", xlab="Age (years)", ylab="Hazard ratio")

# Get anova of the spline too because the p-value indicates whether the spline is important to the model
anova(spline)
# Wald Statistics          Response: Surv(fgstart, fgstop, fgstatus) 
# 
# Factor            Chi-Square d.f. P     
# baseline_age      12324.23    3   <.0001 --> P<0.05 MEANS SPLINE IMPORTANT
# Nonlinear          939.74    2   <.0001
# crp                  22.93    1   <.0001
# gender_encoded      481.23    1   <.0001
# height              201.96    1   <.0001
# weight              140.15    1   <.0001
# smoking_status       55.99    2   <.0001
# diabetes_status       1.28    1   0.2575
# alcohol_status      191.30    2   <.0001
# antiplatelets_use   578.15    1   <.0001
# TOTAL             15709.42   13   <.0001

# Print stored AICs and compare them (with the spline ones displayed earlier too)
# lowest AIC indicates the best transformation.
print(AIC)
#           linear     log   quadr  spline
# weight1 1794856 1794847 1794835   1794824   # spline is lowest
# weight2      NA      NA      NA
# height1 1794856 1794870 1794840   1794809   # spline is lowest
# height2      NA      NA      NA
# crp1    1794856 1794869 1794857   1794811   # spline is lowest
# crp2         NA      NA      NA
# age1    1794856 1794399 1793888   1793882   # spline is lowest 
# age2         NA      NA      NA

# Conclusion: 
# height, weight, CRP, age: spline has the lowest AIC 
################################################################################


################################  FGR FINAL  ###################################
# Make interaction terms directly in the dataset because FGR() does not allow them to be input as rcs(weight,4) in the formula due to a syntax error
train_fg_data$age_rcs <- rcs(train_fg_data$baseline_age, 4)
train_fg_data$weight_rcs <- rcs(train_fg_data$weight, 4)
train_fg_data$height_rcs <- rcs(train_fg_data$height, 4)
train_fg_data$crp_rcs <- rcs(train_fg_data$crp, 4)

train_fg_data$gender_encoded <- as.factor(train_fg_data$gender_encoded)
train_fg_data$alcohol_status <- as.factor(train_fg_data$alcohol_status)
train_fg_data$smoking_status <- as.factor(train_fg_data$smoking_status)
train_fg_data$diabetes_status <- as.factor(train_fg_data$diabetes_status)
train_fg_data$antiplatelets_use <- as.factor(train_fg_data$antiplatelets_use)

# Make the final FGR model 
finegray_final <- FGR(Hist(followup_time,total_outcomes)~age_rcs + age.time +
            weight_rcs + weight.time+
            gender_encoded + gender.1.time + 
            height_rcs +
            crp_rcs + crp.time +
            smoking_status + smoking.1.time + smoking.2.time +
            diabetes_status + diabetes_status.time +
            alcohol_status + alcohol.1.time + alcohol.2.time +
            antiplatelets_use + antiplatelets.time,
          data=train_fg_data, cause = 1)

# Save the model so we can re load it later on (commented is how you read it back in)
saveRDS(finegray_final, "finegray_final.rds") 
#finegray_final <- readRDS("finegray_final.rds")

# HRs and coefs
coef1<- finegray_final$crrFit$coef
coef1<-round((coef1),6)
HR1<-round(cbind(exp(finegray_final$crrFit$coef)), 6)
table1 <- cbind(coef1,HR1); #colnames(table1) <- c("coef","HR","95% CI low","95%CI high") 
colnames(table1) <- c("coef","HR") 

# Make sure test data has same columns as train data so we can make predictions in it 
# for interactions 
test_fg_data$age.time<-test_fg_data$baseline_age*test_fg_data$followup_time
test_fg_data$weight.time<-test_fg_data$weight*test_fg_data$followup_time
test_fg_data$gender.1.time<-ifelse(test_fg_data$gender_encoded==0,test_fg_data$followup_time,0)
test_fg_data$gender.2.time<-ifelse(test_fg_data$gender_encoded==1,test_fg_data$followup_time,0)
test_fg_data$crp.time<-test_fg_data$crp*test_fg_data$baseline_age
test_fg_data$smoking.1.time<-ifelse(test_fg_data$smoking_status==1,test_fg_data$followup_time,0)
test_fg_data$smoking.2.time<-ifelse(test_fg_data$smoking_status==2,test_fg_data$followup_time,0)
test_fg_data$diabetes_status.time<-ifelse(test_fg_data$diabetes_status==1,test_fg_data$followup_time,0)
test_fg_data$alcohol.1.time<-ifelse(test_fg_data$alcohol_status==1,test_fg_data$followup_time,0)
test_fg_data$alcohol.2.time<-ifelse(test_fg_data$alcohol_status==2,test_fg_data$followup_time,0)
test_fg_data$antiplatelets.time<-ifelse(test_fg_data$antiplatelets_use==1,test_fg_data$followup_time,0)
# for transformations
test_fg_data$age_rcs <- rcs(test_fg_data$baseline_age, 4)
test_fg_data$weight_rcs <- rcs(test_fg_data$weight, 4)
test_fg_data$height_rcs <- rcs(test_fg_data$height, 4)
test_fg_data$crp_rcs <- rcs(test_fg_data$crp, 4)
# for time unit
units(test_fg_data$followup_time) <- "Year" 
# categoricals
factorize <-  function (x, columns = names(test_fg_data), izer=as.factor) { test_fg_data[columns] = lapply(test_fg_data[columns], izer);   test_fg_data } 
test_fg_data<-factorize(test_fg_data,c('gender_encoded','alcohol_status','smoking_status','diabetes_status','antiplatelets_use')) 
str(test_fg_data); rm(factorize)
test_fg_data$gender_encoded <- as.factor(test_fg_data$gender_encoded)
test_fg_data$alcohol_status <- as.factor(test_fg_data$alcohol_status)
test_fg_data$smoking_status <- as.factor(test_fg_data$smoking_status)
test_fg_data$diabetes_status <- as.factor(test_fg_data$diabetes_status)
test_fg_data$antiplatelets_use <- as.factor(test_fg_data$antiplatelets_use)

# Make the predictions
# https://cran.r-project.org/web/packages/riskRegression/riskRegression.pdf 
# according to documentation, for predictRisk(), event risks are directly computed
risk <- predictRisk(finegray_final,times=10, newdata=test_fg_data)

# define relevant variables 
test_fg_data$plotrisk <- risk
n.groups <- 10
year.risk <- 10

# C-statistic
surv.obj <- with(test_fg_data, Surv(followup_time,total_outcomes)) # type="mstate" is used for competing risks but requires expansion to as.numeric which results in non-equal length of vars in rcorr
rcorr1 <- rcorr.cens(x=(1-test_fg_data$plotrisk)*10,S=surv.obj)
se.1 <- rcorr1["S.D."]/2
Low95.1 <- rcorr1["C Index"] - 1.96*se.1
Upper95.1 <- rcorr1["C Index"] + 1.96*se.1
c.stat<-cbind(rcorr1["C Index"], Low95.1, Upper95.1); c.stat

# get observed and predicted risks in deciles
test_fg_data$groups1 <- as.numeric(cut2(test_fg_data$plotrisk, g=n.groups, levels.mean=TRUE)) # divide data into equal number of groups based on predicted CV risk
CIF1 <- cuminc(ftime=test_fg_data$followup_time, fstatus=test_fg_data$total_outcomes, cencode="0", group=test_fg_data$groups1)
error1 <- as.numeric(timepoints(CIF1[c(1:n.groups)], year.risk)$var)  # error = variance of cum. inc function of the number of years at risk per group
observeda1 <- as.numeric(timepoints(CIF1[c(1:n.groups)], year.risk)$est)  # observed = point estimate of cum. inc function of the number of years at risk per group
predicteda1<- tapply(test_fg_data$plotrisk, test_fg_data$groups1, mean)  # average of the predicted risk per group

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
legend("topleft", paste("C-statistic: 0.819 (95%CI 0.815 - 0.823)"))

# Cumulative incidence plot
CIF2 <- cuminc(ftime=test_fg_data$followup_time, fstatus=test_fg_data$total_outcomes, cencode="0")
plot.cuminc(CIF2, main="Cumulative Incidence: Fine-Gray Model", curvlab = 1:2,
            ylim=c(0, 1), xlab="Years",  ylab="Probability of Event", lty=1:length(CIF2), color=1:10, lwd=par('lwd'))

################################################################################


############################## OTHER CALCULATIONS ##############################
risk_percent <- risk*100
mean(risk_percent) 
min(risk_percent) 
max(risk_percent) 

# Distribution of Predicted risk
h = hist(risk) 
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, ylim=c(0,60), main = "Distribution of Predicted Risk: Fine-Gray Model", xlab = "Predicted Risk", ylab = "Percent")
################################################################################















