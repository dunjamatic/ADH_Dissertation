### Cox ###

# Clear environment
rm(list=ls())

# Load libraries
library(survival)
library(rms)
library(pec)
library(cmprsk)

# Set working directory
setwd("final/2b. Compare Cox, ML")

# Load in the data
train_cox_data <- read.csv("CoxMLTraining_totalModel.csv")
test_cox_data <- read.csv("CoxMLTesting_totalModel.csv")

# re-name height, weight, crp columns (for convenience)
colnames(train_cox_data)[colnames(train_cox_data) == "height..cm."] ="height"
colnames(train_cox_data)[colnames(train_cox_data) == "weight..kg."] ="weight"
colnames(train_cox_data)[colnames(train_cox_data) == "crp.value..mg.L."] ="crp"
colnames(test_cox_data)[colnames(test_cox_data) == "height..cm."] ="height"
colnames(test_cox_data)[colnames(test_cox_data) == "weight..kg."] ="weight"
colnames(test_cox_data)[colnames(test_cox_data) == "crp.value..mg.L."] ="crp"

################################  COX SET UP  ##################################
# Necessary commands to run to create Cox model
units(train_cox_data$followup_time) <- "Year" 

dd <- datadist(train_cox_data)
options(datadist = "dd")

# Factorize categorical variables so the model knows they are categorical 
# for training data
factorize <-  function (x, columns = names(train_cox_data), izer=as.factor) { train_cox_data[columns] = lapply(train_cox_data[columns], izer);   train_cox_data } 
train_cox_data<-factorize(train_cox_data,c('gender_encoded','alcohol_status','smoking_status','diabetes_status','antiplatelets_use')) 
str(train_cox_data); rm(factorize)
################################################################################


########################### PROPORTIONAL HAZARDS  ##############################
# Checking PH assumption with baseline Cox model to determine if interactions 
# are needed

# Create Cox model
cox_model <- coxph(
  Surv(followup_time, total_outcomes) ~ baseline_age + gender_encoded + height +
    weight + alcohol_status + smoking_status + diabetes_status + antiplatelets_use + crp,
  data = train_cox_data
)

# Now check PH assumption 
# first with plots 
somePDFPath1 = "proportionalHazards_Cox_trainModel.pdf"
pdf(file=somePDFPath1)
plot(cox.zph(cox_model), resid = F, se = T, lwd = 3)
dev.off()
rm(somePDFPath1)
# now with P values 
cox.zph(cox_model)

# cox.zph calculates Schoenfeld residuals which are the difference between the observed
# and expected covariate at different time points. a test between the residuals and time
# is then done. it computes rho, which measures the association between the covariate and time.
# null hypothesis of the test is that the correlation coefficient is zero, indicating
# no association between the covariate and time. p < 0.05 suggests violation of PH.

# plot of cox.zph lets you visually assess too. if there are plots that are not a horizontal line
# it might suggest violations of the PH assumption.

#                   chisq df       p
# baseline_age      223.53  1 < 2e-16 ### need to add interaction 
# gender_encoded      6.73  1  0.0095 ### need to add interaction 
# height              1.24  1  0.2662
# weight             36.37  1 1.6e-09 ### need to add interaction 
# alcohol_status    144.71  2 < 2e-16 ### need to add interaction 
# smoking_status     41.34  2 1.1e-09 ### need to add interaction 
# diabetes_status    81.60  1 < 2e-16 ### need to add interaction 
# antiplatelets_use 419.00  1 < 2e-16 ### need to add interaction 
# crp                 7.43  1  0.0064 ### need to add interaction 
# GLOBAL            765.38 11 < 2e-16
################################################################################


##############################  INTERACTIONS  ##################################
# interactions created are based on what the PH assumption results showed 
# establish interactions with age
train_cox_data$age.time<-train_cox_data$baseline_age*train_cox_data$followup_time

# establish interactions with gender
train_cox_data$gender.1.time<-ifelse(train_cox_data$gender_encoded==1,train_cox_data$followup_time,0)

# establish interactions with weight
train_cox_data$weight.time<-train_cox_data$weight*train_cox_data$followup_time

# establish interaction with alcohol
train_cox_data$alcohol.1.time<-ifelse(train_cox_data$alcohol_status==1,train_cox_data$followup_time,0)
train_cox_data$alcohol.2.time<-ifelse(train_cox_data$alcohol_status==2,train_cox_data$followup_time,0)

# establish interaction with smoking
train_cox_data$smoking.1.time<-ifelse(train_cox_data$smoking_status==1,train_cox_data$followup_time,0)
train_cox_data$smoking.2.time<-ifelse(train_cox_data$smoking_status==2,train_cox_data$followup_time,0)

# establish interaction with diabetes
train_cox_data$diabetes.time<-ifelse(train_cox_data$diabetes_status==1,train_cox_data$followup_time,0)

# establish interaction with antiplatelets
train_cox_data$antiplatelets.time<-ifelse(train_cox_data$antiplatelets_use==1,train_cox_data$followup_time,0)

# establish interactions with crp
train_cox_data$crp.time<-train_cox_data$crp*train_cox_data$followup_time
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
mod1 <- cph(Surv(followup_time, total_outcomes)~weight+baseline_age+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = train_cox_data, x = T, y = T)
# try log
mod2 <- cph(Surv(followup_time, total_outcomes)~log(weight+1)+baseline_age+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = train_cox_data, x = T, y = T)
# try quadr
mod3 <- cph(Surv(followup_time, total_outcomes)~pol(weight,2)+baseline_age+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = train_cox_data, x = T, y = T)
# store AIC
AIC[1,1] <- extractAIC(mod1)[2]
AIC[1,2] <- extractAIC(mod2)[2]
AIC[1,3] <- extractAIC(mod3)[2]

# try spline
spline <- cph(Surv(followup_time, total_outcomes)~rcs(weight,4)+baseline_age+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
              data=train_cox_data,x=T,y=T)
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
# weight              149.30    3   <.0001 --> p < 0.05 means SPLINE IMPORTANT IN MODEL!
# Nonlinear            25.54    2   <.0001
# baseline_age      15694.43    1   <.0001
# gender_encoded      585.43    1   <.0001
# height              182.75    1   <.0001
# crp                   2.54    1   0.1107
# smoking_status       35.67    2   <.0001
# diabetes_status      34.01    1   <.0001
# alcohol_status      276.36    2   <.0001
# antiplatelets_use   844.40    1   <.0001
# TOTAL             20114.99   13   <.000

# HEIGHT
# try linear
mod1 <- cph(Surv(followup_time, total_outcomes)~weight+baseline_age+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = train_cox_data, x = T, y = T)
# try log
mod2 <- cph(Surv(followup_time, total_outcomes)~weight+baseline_age+gender_encoded+log(height+1)+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = train_cox_data, x = T, y = T)
# try quadr
mod3 <- cph(Surv(followup_time, total_outcomes)~weight+baseline_age+gender_encoded+pol(height,2)+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = train_cox_data, x = T, y = T)
# store AIC
AIC[3,1] <- extractAIC(mod1)[2]
AIC[3,2] <- extractAIC(mod2)[2]
AIC[3,3] <- extractAIC(mod3)[2]

# try spline
spline <- cph(Surv(followup_time, total_outcomes)~rcs(height,4)+baseline_age+gender_encoded+weight+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=train_cox_data,x=T,y=T)
AIC[3,4] <- AIC(spline) 

# Plot spline because it will show linearity (or lack thereof) 
plot(Predict(spline,height, fun=exp), conf.int=T, main="Spline of Height", xlab="Height (cm)", ylab="Hazard ratio")

# Get anova of the spline too because the p-value indicates whether the spline is important to the model
anova(spline) # 
# Wald Statistics          Response: Surv(followup_time, total_outcomes) 
# 
# Factor            Chi-Square d.f. P     
# height              247.33    3   <.0001 --> p < 0.05 means SPLINE IMPORTANT IN MODEL!
# Nonlinear           49.29    2   <.0001
# baseline_age      15704.18    1   <.0001
# gender_encoded      402.39    1   <.0001
# weight              115.39    1   <.0001
# crp                   2.55    1   0.1106
# smoking_status       35.72    2   <.0001
# diabetes_status      35.70    1   <.0001
# alcohol_status      276.70    2   <.0001
# antiplatelets_use   850.43    1   <.0001
# TOTAL             20135.56   13   <.0001

# CRP
# try linear
mod1 <- cph(Surv(followup_time, total_outcomes)~weight+baseline_age+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = train_cox_data, x = T, y = T)
# try log
mod2 <- cph(Surv(followup_time, total_outcomes)~weight+baseline_age+gender_encoded+height+log(crp+1)+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = train_cox_data, x = T, y = T)
# try quadr
mod3 <- cph(Surv(followup_time, total_outcomes)~weight+baseline_age+gender_encoded+height+pol(crp,2)+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = train_cox_data, x = T, y = T)
# store AIC
AIC[5,1] <- extractAIC(mod1)[2]
AIC[5,2] <- extractAIC(mod2)[2]
AIC[5,3] <- extractAIC(mod3)[2]

# try spline
spline <- cph(Surv(followup_time, total_outcomes)~rcs(crp,4)+baseline_age+gender_encoded+height+weight+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=train_cox_data,x=T,y=T)
AIC[5,4] <- AIC (spline) 

# Plot spline because it will show linearity (or lack thereof) 
plot(Predict(spline,crp, fun=exp), conf.int=T, main="Spline of CRP", xlab="CRP (mg/L)", ylab="Hazard ratio")

# Get anova of the spline too because the p-value indicates whether the spline is important to the model
anova(spline) 
# Wald Statistics          Response: Surv(followup_time, total_outcomes) 
# 
# Factor            Chi-Square d.f. P   
# crp                  52.13    3   <.0001 --> p < 0.05 means SPLINE IMPORTANT IN MODEL!
# Nonlinear           49.45    2   <.0001
# baseline_age      15619.75    1   <.0001
# gender_encoded      610.62    1   <.0001
# height              198.34    1   <.0001
# weight              124.13    1   <.0001
# smoking_status       36.05    2   <.0001
# diabetes_status      33.69    1   <.0001
# alcohol_status      279.20    2   <.0001
# antiplatelets_use   853.16    1   <.0001
# TOTAL             20130.34   13   <.0001

# AGE
# try linear
mod1 <- cph(Surv(followup_time, total_outcomes)~weight+baseline_age+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = train_cox_data, x = T, y = T)
# try log
mod2 <- cph(Surv(followup_time, total_outcomes)~weight+log(baseline_age+1)+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = train_cox_data, x = T, y = T)
# try quadr
mod3 <- cph(Surv(followup_time, total_outcomes)~weight+pol(baseline_age,2)+gender_encoded+height+crp+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,
            data = train_cox_data, x = T, y = T)
# store AIC
AIC[7,1] <- extractAIC(mod1)[2]
AIC[7,2] <- extractAIC(mod2)[2]
AIC[7,3] <- extractAIC(mod3)[2]

# try spline
spline <- cph(Surv(followup_time, total_outcomes)~rcs(baseline_age,4)+crp+gender_encoded+height+weight+smoking_status+diabetes_status+alcohol_status+antiplatelets_use,data=train_cox_data,x=T,y=T)
AIC[7,4] <- AIC (spline)  

# Plot spline because it will show linearity (or lack thereof) 
plot(Predict(spline,baseline_age, fun=exp), conf.int=T, main="Spline of Age", xlab="Age (years)", ylab="Hazard ratio")

# Get anova of the spline too because the p-value indicates whether the spline is important to the model
anova(spline) 
# Wald Statistics          Response: Surv(followup_time, total_outcomes) 
# 
# Factor            Chi-Square d.f. P 
# baseline_age      14733.45    3   <.0001 --> p < 0.05 means SPLINE IMPORTANT IN MODEL!
# Nonlinear          543.74    2   <.0001
# crp                   2.17    1   0.1405
# gender_encoded      584.91    1   <.0001
# height              195.12    1   <.0001
# weight              124.99    1   <.0001
# smoking_status       35.57    2   <.0001
# diabetes_status      28.91    1   <.0001
# alcohol_status      270.25    2   <.0001
# antiplatelets_use   858.77    1   <.0001
# TOTAL             19012.58   13   <.0001

# Print stored AICs and compare them (with the spline ones displayed earlier too)
# lowest AIC indicates the best transformation.
print(AIC)
#           linear     log   quadr  spline
# weight1 1783874 1783870 1783863 1783852 # spline is lowest
# weight2      NA      NA      NA      NA
# height1 1783874 1783888 1783857 1783829 # spline is lowest
# height2      NA      NA      NA      NA
# crp1    1783874 1783874 1783873 1783828 # spline is lowest
# crp2         NA      NA      NA      NA
# age1    1783874 1783517 1783314 1783315 # quadr is lowest
# age2         NA      NA      NA      NA

# Conclusion: 
# height, weight, CRP: spline has the lowest AIC 
# age: quadr has the lowest AIC 
################################################################################


################################  CPH FINAL  ###################################
# Make interaction terms directly in the dataset to be on same ground for Fine-Gray 
train_cox_data$age_quadr <- pol(train_cox_data$baseline_age, 2)
train_cox_data$weight_rcs <- rcs(train_cox_data$weight, 4)
train_cox_data$height_rcs <- rcs(train_cox_data$height, 4)
train_cox_data$crp_rcs <- rcs(train_cox_data$crp, 4)

# Make the final Cox model 
# need to run coxph to be able to run predict()
cox_model_final <- coxph(
  Surv(followup_time, total_outcomes) ~ age_quadr + age.time + 
    gender_encoded + gender.1.time + 
    height_rcs + 
    weight_rcs + weight.time + 
    alcohol_status + alcohol.1.time + alcohol.2.time +
    smoking_status + smoking.1.time + smoking.2.time +
    diabetes_status + diabetes.time +
    antiplatelets_use + antiplatelets.time +
    crp_rcs + crp.time,
  data = train_cox_data
)

# Save the model so we can re load it later on (commented is how you read it back in)
saveRDS(cox_model_final, "cox_model_final.rds") 
cox_model_final <- readRDS("cox_model_final.rds")

# HRs and coefs
coef1<- cox_model_final$coefficients
coef1<-round((coef1),6)
HR1<-round(cbind(exp(cox_model_final$coefficients),exp(confint(cox_model_final))) ,6)
table1 <- cbind(coef1,HR1); colnames(table1) <- c("coef","HR","95% CI low","95%CI high") 
table1

# Make categorical variables in test data
# and for testing data
factorize <-  function (x, columns = names(test_cox_data), izer=as.factor) { test_cox_data[columns] = lapply(test_cox_data[columns], izer);   test_cox_data } 
test_cox_data<-factorize(test_cox_data,c('gender_encoded','alcohol_status','smoking_status','diabetes_status','antiplatelets_use')) 
str(test_cox_data); rm(factorize)

# Make sure test data has same columns as train data so we can make predictions in it 
# for interactions
test_cox_data$age.time<-test_cox_data$baseline_age*test_cox_data$followup_time
test_cox_data$age2.time<-(test_cox_data$baseline_age^2)*test_cox_data$baseline_age
test_cox_data$gender.1.time<-ifelse(test_cox_data$gender_encoded==1,test_cox_data$followup_time,0)
test_cox_data$weight.time<-test_cox_data$weight*test_cox_data$followup_time
test_cox_data$smoking.1.time<-ifelse(test_cox_data$smoking_status==1,test_cox_data$followup_time,0)
test_cox_data$smoking.2.time<-ifelse(test_cox_data$smoking_status==2,test_cox_data$followup_time,0)
test_cox_data$alcohol.1.time<-ifelse(test_cox_data$alcohol_status==1,test_cox_data$followup_time,0)
test_cox_data$alcohol.2.time<-ifelse(test_cox_data$alcohol_status==2,test_cox_data$followup_time,0)
test_cox_data$diabetes.time<-ifelse(test_cox_data$diabetes_status==1,test_cox_data$followup_time,0)
test_cox_data$antiplatelets.time<-ifelse(test_cox_data$antiplatelets_use==1,test_cox_data$followup_time,0)
test_cox_data$crp.time<-test_cox_data$crp*test_cox_data$followup_time
# for transformations
test_cox_data$age_quadr <- pol(test_cox_data$baseline_age, 2)
test_cox_data$weight_rcs <- rcs(test_cox_data$weight, 4)
test_cox_data$height_rcs <- rcs(test_cox_data$height, 4)
test_cox_data$crp_rcs <- rcs(test_cox_data$crp, 4)
# also make sure test units are in Years like train
units(test_cox_data$followup_time) <- "Year"

# Make the predictions
# get predicted survivals and turn it into risk
exp_event <- predict(cox_model_final,newdata=test_cox_data, type="expected")
surv_prob <- exp(-exp_event)
risk <- (1-surv_prob)*100 

# define relevant variables
test_cox_data$plotrisk <- risk/100
n.groups <- 10
year.risk <- 10

# C-statistic
surv.obj <- with(test_cox_data,Surv(followup_time,total_outcomes))
rcorr1 <- rcorr.cens(x=(1-test_cox_data$plotrisk)*10,S=surv.obj)
se.1 <- rcorr1["S.D."]/2
Low95.1 <- rcorr1["C Index"] - 1.96*se.1
Upper95.1 <- rcorr1["C Index"] + 1.96*se.1
c.stat<-cbind(rcorr1["C Index"], Low95.1, Upper95.1); c.stat

# get predicted and observed risks in deciles 
test_cox_data$groups1 <- as.numeric(cut2(test_cox_data$plotrisk, g=n.groups, levels.mean=TRUE)) # divide data into equal number of groups based on predicted CV risk
CIF1 <- cuminc(ftime=test_cox_data$followup_time, fstatus=test_cox_data$total_outcomes, cencode="0", group=test_cox_data$groups1)
error1 <- as.numeric(timepoints(CIF1[c(1:n.groups)], year.risk)$var)  # error = variance of cum. inc function of the number of years at risk per group
observeda1 <- as.numeric(timepoints(CIF1[c(1:n.groups)], year.risk)$est)  # observeda = point estimate of cum. inc function of the number of years at risk per group
predicteda1<- tapply(test_cox_data$plotrisk, test_cox_data$groups1, mean)  # average of the predicted risk per group

# Calibration plot
plot(predicteda1, observeda1, xlab=" ", type = "b", ylab=" ", pch = 15, xaxt = "n",yaxt = "n",xlim=c(0,0.4),ylim=c(0,0.4), cex.axis=1, cex.lab=1)
title("Calibration: Cox Model", cex.main=1.2,line=1)
errbar(predicteda1,observeda1, (exp(log(observeda1+1.96*sqrt(error1)))), (exp(log(observeda1-1.96*sqrt(error1)))),add=T, lwd = 0.5)
abline(a=0,b=1, lty=4)
lines(predicteda1,observeda1)
mtext("Predicted 10-year risk",cex=1.2,side=1,adj=0.5,padj=+3.5)
mtext("Observed 10-year risk",cex=1.2,side=2,adj=0.5,padj=-3.5)
axis(1, at=c(0:10/10), labels=paste(c(0:10*10),"%",sep=""), lwd = 1, cex.axis=1.2)
axis(2, at=c(0:10/10), labels=paste(c(0:10*10),"%",sep=""), lwd = 1, cex.axis=1.2)
legend("topleft", paste("C-statistic: 0.920 (95%CI 0.917 - 0.923)")) 

# Cumulative incidence plot
CIF2 <- cuminc(ftime=test_cox_data$followup_time, fstatus=test_cox_data$total_outcomes, cencode="0")
plot.cuminc(CIF2, main="Cumulative Incidence: Cox Model", curvlab = 1,
            ylim=c(0, 1), xlab="Years",  ylab="Probability of Cancer Diagnosis", lty=1:length(CIF2), color=1, lwd=par('lwd'))

dev.off()
################################################################################


############################## OTHER CALCULATIONS ##############################
# Mean predicted risk and range 
risk_percent <- risk
mean(risk_percent) 
min(risk_percent) 
max(risk_percent) 

# Distribution of Predicted risk
h = hist(risk) # or hist(x,plot=FALSE) to avoid the plot of the histogram
h$density = h$counts/sum(h$counts)*100
plot(h, freq=FALSE, ylim=c(0,70), main = "Distribution of Predicted Risk: Cox Model", xlab = "Predicted Risk", ylab = "Percent")
################################################################################



