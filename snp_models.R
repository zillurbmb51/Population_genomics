library(readr)
library(plyr)

# Call each of my interested individuals files:

# To delete scientific notation in data (To activate use 0 instead of 999)(scipen=999) [issues with ancestry data sciebtific notation]
options(scipen = 999)

# Instructions to merge small tables in bigger table called Data
Ancestry=read.csv(file="Ancestry.csv", header=T)
head(Ancestry)

Phenotypes=read.csv(file="Phenotypes 2.csv", header=T)
head(Phenotypes)

Final_Id=read.csv(file="Final_List_May_2018.csv", header=T)
head(Final_Id)

Scholarity=read.csv(file="Scholarity_1.csv", header=T)
head(Scholarity)

Pathologies=read.csv(file="Pathologies 2.csv", header=T)
head(Pathologies)


#To join my individual files in one general file

Data=join(Final_Id,Phenotypes, by="Sample_ID", type="left")
head(Data)

Data=join(Data, Ancestry, by="Sample_ID", type="left")
head(Data)

Data=join(Data, Pathologies, by="Sample_ID", type="left")
head(Data)

Data=join(Data, Scholarity, by="Sample_ID", type="left")
head(Data)

summary(Data)

#NOT STRATIFIED:

#rs6693065

#DOM
Data$DX=as.factor(Data$DX)
univariate=glm(relevel(DX, "CON") ~ relevel(rs6693065_D,"AA"), family = binomial, data = Data)
summary(univariate)
exp(cbind(OR = coef(univariate), confint(univariate)))

#ADD
Data$DX=as.factor(Data$DX)
univariate=glm(relevel(DX, "CON") ~ relevel(rs6693065_A,"AA"), family = binomial, data = Data)
summary(univariate)
exp(cbind(OR = coef(univariate), confint(univariate)))


#REC
Data$DX=as.factor(Data$DX)
univariate=glm(relevel(DX, "CON") ~ relevel(rs6693065_R,"OTHER"), family = binomial, data = Data)
summary(univariate)
exp(cbind(OR = coef(univariate), confint(univariate)))


#ALL variables alleles
Data$DX=as.factor(Data$DX)
univariate=glm(relevel(DX, "CON") ~ relevel(rs6693065_D,"AA") + Age + Scholarity + ERT + Breast_Fed.Dur, family = binomial, data = Data)
summary(univariate)
exp(cbind(OR = coef(univariate), confint(univariate)))


# without few

Data$DX=as.factor(Data$DX)
univariate=glm(relevel(DX, "CON") ~ relevel(rs6693065_D, "AA") + Age + Scholarity, family = binomial, data = Data)
summary(univariate)
exp(cbind(OR = coef(univariate), confint(univariate)))


#MODELO
Data$DX=as.factor(Data$DX)
FITALL=glm(relevel(DX, "CON") ~ relevel(rs6693065_D,"AA") + Age + Scholarity + Num_Child + ERT + EUR + NAT + Meno_status + Geographic_Region_code + Anticonc + Civil_status_code + Numb_Preg + Breast_Fed.Dur + BMI, family = binomial, data = Data)
summary(FITALL)
exp(cbind(OR = coef(FITALL), confint(FITALL)))
step(FITALL, direction = "backward")


#ALL variables alleles not best model yet
Data$DX=as.factor(Data$DX)
univariate=glm(relevel(DX, "CON") ~ relevel(rs6693065_D,"AA") + Age + EUR + Scholarity + ERT + Breast_Fed.Dur + BMI + Numb_Preg, family = binomial, data = Data)
summary(univariate)
exp(cbind(OR = coef(univariate), confint(univariate)))

# without few not best model yet

#ALL variables alleles new model
Data$DX=as.factor(Data$DX)
univariate=glm(relevel(DX, "CON") ~ relevel(rs6693065_D,"AA") + Age + Scholarity + ERT + Breast_Fed.Dur, family = binomial, data = Data)
summary(univariate)
exp(cbind(OR = coef(univariate), confint(univariate)))



# Regression Stratified

#SNPs Dominant model analysis# Invasive vs Controls

#1. rs6693065

Data$DX_Invasive[Data$Type == "Invasive"]<-"Invasive"
Data$DX_Invasive[Data$Type == "In situ"]<-NA
Data$DX_Invasive[Data$DX == "CON"]<-"CON"

Data$DX_Invasive=as.factor(Data$DX_Invasive)
univariate=glm(relevel(DX_Invasive, "CON") ~ relevel(rs6693065_D, "AA"), family = binomial, data = Data)
summary(univariate)
exp(cbind(OR = coef(univariate), confint(univariate)))

Data$DX_Invasive=as.factor(Data$DX_Invasive)
univariate=glm(relevel(DX_Invasive, "CON") ~ relevel(rs6693065_D,"AA") + Age + Scholarity, family = binomial, data = Data)
summary(univariate)
exp(cbind(OR = coef(univariate), confint(univariate)))

Data$DX_Invasive=as.factor(Data$DX_Invasive)
univariate=glm(relevel(DX_Invasive, "CON") ~ relevel(rs6693065_D,"AA") + Age + Scholarity + ERT + Breast_Fed.Dur, family = binomial, data = Data)
summary(univariate)
exp(cbind(OR = coef(univariate), confint(univariate)))


Data$DX_Invasive=as.factor(Data$DX_Invasive)
modall=glm(relevel(DX_Invasive, "CON") ~ relevel(rs6693065_D, "AA") + Age + EUR + Scholarity + ERT + Num_Child + Numb_Preg + Breast_Fed.Dur + Civil_status_code + Occupation_code, family = binomial, data = Data)
summary(modall)
exp(cbind(OR = coef(modall), confint(modall)))

Data$DX_Invasive=as.factor(Data$DX_Invasive)
univariate=glm(relevel(DX_Invasive, "CON") ~ relevel(rs6693065_D,"AA") + Age + Scholarity, family = binomial, data = Data)
summary(univariate)
exp(cbind(OR = coef(univariate), confint(univariate)))

Data$DX_Invasive=as.factor(Data$DX_Invasive)
univariate=glm(relevel(DX_Invasive, "CON") ~ relevel(rs6693065_D,"AA") + Age + Scholarity + ERT + Breast_Fed.Dur, family = binomial, data = Data)
summary(univariate)
exp(cbind(OR = coef(univariate), confint(univariate)))


###ADD Nuevo Not strat

#NOT STRATIFIED CRUDE
#rs6693065

Data$DX=as.factor(Data$DX)
univariate=glm(relevel(DX, "CON") ~ rs6693065_A, family = binomial, data = Data)
summary(univariate)
exp(cbind(OR = coef(univariate), confint(univariate)))

#rs280500

Data$DX=as.factor(Data$DX)
univariate=glm(relevel(DX, "CON") ~ rs280500_A, family = binomial, data = Data)
summary(univariate)
exp(cbind(OR = coef(univariate), confint(univariate)))

#For all snps together
proc_glm <- function(snps) {
  univariate <- glm(relevel(Data$DX, "CON") ~ relevel(snps, "AA"), family = binomial)
  
  return(exp(cbind(OR = coef(univariate), confint(univariate))))
}

glm_list2 <- sapply(Data[3:426], function(col)  
  tryCatch(proc_glm(col), error = function(e) NA))
glm_list2
