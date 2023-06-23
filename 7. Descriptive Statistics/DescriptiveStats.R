rm(list=ls())
library(tableone)
library(survival)
library(ggplot2)

############################## PRE-IMPUTATION ##################################
# Create Boxplots for weight, height, CRP to Examine Outliers
# get data from pre-clipping
setwd("final\\0. Imputation") 
dat2 <- read.csv('imputation_totalModel.csv')

# re-name height, weight, crp columns
colnames(dat2)[colnames(dat2) == "height..cm."] ="height"
colnames(dat2)[colnames(dat2) == "weight..kg."] ="weight"
colnames(dat2)[colnames(dat2) == "crp.value..mg.L."] ="crp"

# Height 
boxplot(dat2$height,
        ylab = "Height (cm)",
        main = "Boxplot of Height Prior to Clipping"
)

ggplot(dat2, aes(height)) +
  ggtitle("Histogram of Height Prior to Clipping") +
  geom_histogram(bins = 30L, fill = "#0c4c8a") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) + 
  labs(
    x = "Height (cm)",
    y = "Frequency"
  )

# Weight
boxplot(dat2$weight,
        ylab = "Weight (kg)",
        main = "Boxplot of Weight Prior to Clipping"
)

ggplot(dat2) +
  aes(x = weight) +
  geom_histogram(bins = 30L, fill = "#0c4c8a") +
  ggtitle("Histogram of Weight Prior to Clipping") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "Weight (kg)",
    y = "Frequency"
  )

# CRP
boxplot(dat2$crp,
        ylab = "CRP (mg/L)",
        main = "Boxplot of CRP Prior to Clipping"
)

ggplot(dat2) +
  aes(x = crp) +
  geom_histogram(bins = 30L, fill = "#0c4c8a") + 
  ggtitle("Histogram of CRP Prior to Clipping") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "CRP (mg/L)",
    y = "Frequency"
  )

dev.off()
################################################################################


############################### POST IMPUTATION ################################

################################### OBJECTIVE 1 ######################################
######################### Create Table One - OBJECTIVE 1 #############################
# Import data
dat <- read.csv("imputed_totalModel.csv")

# re-name height, weight, crp columns
colnames(dat)[colnames(dat) == "height..cm."] ="height (cm)"
colnames(dat)[colnames(dat) == "weight..kg."] ="weight (kg)"
colnames(dat)[colnames(dat) == "crp.value..mg.L."] ="CRP (mg/mL)"

# Indicate variables to include in the table
vars <- c("baseline_age", "gender_encoded", "ethnicity_encoded", "townsend", "height (cm)", "weight (kg)", "alcohol_status", "smoking_status", "diabetes_status", "antiplatelets_use", "CRP (mg/mL)", "total_outcomes") 

# Indicate categorical variables 
catVars <- c("ethnicity_encoded", "gender_encoded", "alcohol_status", "smoking_status", "diabetes_status", "antiplatelets_use", "total_outcomes")

# Create table
tab <- CreateTableOne(vars = vars, data = dat, factorVars = catVars)

# Print table 
out <- print(tab, showAllLevels = TRUE, formatOptions = list(big.mark = ","), smd = TRUE)


# Save tableone output to the working directory 
write.csv(out, file = "final/0. Imputation/TableOnes/imputedTotalModel_tableone.csv")

####### Distributions for Continuous Variables (post-imputation) - OBJECTIVE 1 #######
setwd("final/0. Imputation") 
dat <- read.csv("imputed_totalModel.csv")

# re-name height, weight, crp columns
colnames(dat)[colnames(dat) == "height..cm."] ="height"
colnames(dat)[colnames(dat) == "weight..kg."] ="weight"
colnames(dat)[colnames(dat) == "crp.value..mg.L."] ="CRP"

# height
ggplot(dat, aes(height)) +                               
  ggtitle("Distribution of Height") +
  geom_histogram(color = "black", fill = "white") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "Height (cm)",
    y = "Frequency"
  )

# weight
ggplot(dat, aes(weight)) +                               
  ggtitle("Distribution of Weight") +
  geom_histogram(color = "black", fill = "white") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "Weight (kg)",
    y = "Frequency"
  )

# CRP
ggplot(dat, aes(CRP)) +                               
  ggtitle("Distribution of CRP") +
  geom_histogram(color = "black", fill = "white") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "CRP (mg/L)",
    y = "Frequency"
  ) 

# Townsend
ggplot(dat, aes(townsend)) +                               
  ggtitle("Distribution of Towsend Deprivation Score") +
  geom_histogram(color = "black", fill = "white") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "Towsend Deprivation Score",
    y = "Frequency"
  ) 

dev.off()

############################# Median follow-up time ############################
median(dat$followup_time) 

###################### Event rates per 100 person years ########################
dat$death <- dat$total_outcomes
dat$death[dat$death == 1] <- 0
dat$death[dat$death == 2] <- 1

dat$cancer <- dat$total_outcomes
dat$cancer[dat$cancer == 2] <- 0

# death rates
deathrates <- pyears(Surv(followup_time, death) ~ 1, scale =100, data=dat)
round(deathrates$event/ deathrates$pyears,2)  

# cancer rates
cancerrates <- pyears(Surv(followup_time, cancer) ~ 1, scale=100, data=dat)
round(cancerrates$event/ cancerrates$pyears,2) 


################################ OBJECTIVES 2 AND 3 ##################################
################### Create table one - OBJECTIVES TWO AND THREE ######################

# Import data
setwd("final/2a. Compare FG to ML") 
dat_p1 <- read.csv("FGTraining_totalModel.csv")
dat_p2 <- read.csv("FGTesting_totalModel.csv")

# re-name height, weight, crp columns
colnames(dat_p1)[colnames(dat_p1) == "height..cm."] ="height (cm)"
colnames(dat_p1)[colnames(dat_p1) == "weight..kg."] ="weight (kg)"
colnames(dat_p1)[colnames(dat_p1) == "crp.value..mg.L."] ="CRP (mg/mL)"
colnames(dat_p2)[colnames(dat_p2) == "height..cm."] ="height (cm)"
colnames(dat_p2)[colnames(dat_p2) == "weight..kg."] ="weight (kg)"
colnames(dat_p2)[colnames(dat_p2) == "crp.value..mg.L."] ="CRP (mg/mL)"

# Indicate variables to include in the table
vars <- c("baseline_age", "gender_encoded", "ethnicity_encoded", "townsend", "height (cm)", "weight (kg)", "alcohol_status", "smoking_status", "diabetes_status", "antiplatelets_use", "CRP (mg/mL)", "total_outcomes") 

# Indicate categorical variables 
catVars <- c("ethnicity_encoded", "gender_encoded", "alcohol_status", "smoking_status", "diabetes_status", "antiplatelets_use", "total_outcomes")

# Create table
tab_p1 <- CreateTableOne(vars = vars, data = dat_p1, factorVars = catVars)
tab_p2 <- CreateTableOne(vars = vars, data = dat_p2, factorVars = catVars)

# Print table 
out_p1 <- print(tab_p1, showAllLevels = TRUE, formatOptions = list(big.mark = ","), smd = TRUE)
out_p2 <- print(tab_p2, showAllLevels = TRUE, formatOptions = list(big.mark = ","), smd = TRUE)

# Save tableone output to the working directory 
write.csv(out_p1, file = "final/0. Imputation/TableOnes/objective2_3_tableone_train.csv")
write.csv(out_p2, file = "final/0. Imputation/TableOnes/objective2_3_tableone_test.csv")

#### Distributions for Continuous Variables (post-imputation) - OBJECTIVE 2 and 3 ####
dat_p1 <- read.csv("FGTraining_totalModel.csv")
dat_p2 <- read.csv("FGTesting_totalModel.csv")

# re-name height, weight, crp columns
colnames(dat_p1)[colnames(dat_p1) == "height..cm."] ="height"
colnames(dat_p1)[colnames(dat_p1) == "weight..kg."] ="weight"
colnames(dat_p1)[colnames(dat_p1) == "crp.value..mg.L."] ="CRP"
colnames(dat_p2)[colnames(dat_p2) == "height..cm."] ="height"
colnames(dat_p2)[colnames(dat_p2) == "weight..kg."] ="weight"
colnames(dat_p2)[colnames(dat_p2) == "crp.value..mg.L."] ="CRP"

# height p1
ggplot(dat_p1, aes(height)) +                               
  ggtitle("Distribution of Height - Training Data") +
  geom_histogram(color = "black", fill = "white") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "Height (cm)",
    y = "Frequency"
  )
# height p2
ggplot(dat_p2, aes(height)) +                               
  ggtitle("Distribution of Height - Testing Data") +
  geom_histogram(color = "black", fill = "white") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "Height (cm)",
    y = "Frequency"
  )

# weight p1
ggplot(dat_p1, aes(weight)) +                               
  ggtitle("Distribution of Weight - Training Data") +
  geom_histogram(color = "black", fill = "white") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "Weight (kg)",
    y = "Frequency"
  )
# weight p2
ggplot(dat_p2, aes(weight)) +                               
  ggtitle("Distribution of Weight - Testing Data") +
  geom_histogram(color = "black", fill = "white") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "Weight (kg)",
    y = "Frequency"
  )

# CRP p1
ggplot(dat_p1, aes(CRP)) +                               
  ggtitle("Distribution of CRP - Training Data") +
  geom_histogram(color = "black", fill = "white") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "CRP (mg/L)",
    y = "Frequency"
  ) 
# CRP p2
ggplot(dat_p2, aes(CRP)) +                               
  ggtitle("Distribution of CRP - Testing Data") +
  geom_histogram(color = "black", fill = "white") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "CRP (mg/L)",
    y = "Frequency"
  ) 

# Townsend p1
ggplot(dat_p1, aes(townsend)) +                               
  ggtitle("Distribution of Towsend Deprivation Score - Training Data") +
  geom_histogram(color = "black", fill = "white") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "Towsend Deprivation Score",
    y = "Frequency"
  ) 
# Townsend p2
ggplot(dat_p2, aes(townsend)) +                               
  ggtitle("Distribution of Towsend Deprivation Score - Testing Data") +
  geom_histogram(color = "black", fill = "white") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "Towsend Deprivation Score",
    y = "Frequency"
  ) 

dev.off()

############################# Median follow-up time ############################
# p1
median(dat_p1$followup_time) 
# p2
median(dat_p2$followup_time) 

###################### Event rates per 100 person years ########################
# p1
dat_p1$death <- dat_p1$total_outcomes
dat_p1$death[dat_p1$death == 1] <- 0
dat_p1$death[dat_p1$death == 2] <- 1

dat_p1$cancer <- dat_p1$total_outcomes
dat_p1$cancer[dat_p1$cancer == 2] <- 0

# death rates
deathrates <- pyears(Surv(followup_time, death) ~ 1, scale =100, data=dat_p1)
round(deathrates$event/ deathrates$pyears,2)  

# cancer rates
cancerrates <- pyears(Surv(followup_time, cancer) ~ 1, scale=100, data=dat_p1)
round(cancerrates$event/ cancerrates$pyears,2) 


# p2
dat_p2$death <- dat_p2$total_outcomes
dat_p2$death[dat_p2$death == 1] <- 0
dat_p2$death[dat_p2$death == 2] <- 1

dat_p2$cancer <- dat_p2$total_outcomes
dat_p2$cancer[dat_p2$cancer == 2] <- 0

# death rates
deathrates <- pyears(Surv(followup_time, death) ~ 1, scale =100, data=dat_p2)
round(deathrates$event/ deathrates$pyears,2)  

# cancer rates
cancerrates <- pyears(Surv(followup_time, cancer) ~ 1, scale=100, data=dat_p2)
round(cancerrates$event/ cancerrates$pyears,2)



