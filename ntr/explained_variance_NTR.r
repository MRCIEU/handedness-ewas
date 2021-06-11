# Explained Variance (R2) calculation for left-handedness (binomial trait)
# V.Odintsova, J.van Dongen, C.Dolan
# Based on Lee et al.
# Genetic Epidemiology 36 : 214–224 (2012)
# A Better Coefficient of Determination for Genetic Profile Analysis

########################### ############## ##############
# Load data
########################### ############## ##############
# polygenic risk scores (PGS)
# methylation scores (MS)
# EWAS covariates
# PGS covariates

########################### ############## ##############
# Select variables for the prediction models
########################### ############## ##############
# FISNumber
# FamilyNumber
# sex
# age (age at DNA bloodsampling, scaled)
# hand01 (left-handedness variable: right-handed = 0, left-handed = 1)

# PGS and GWAS covariates
# PGS
# PC1 - PC10 (principal components)
# platform2 - platform6 (platform dummy variables)

# Methylation scores and EWAS covariates from adjusted model
# MS1 (methylation score calculated with weights at p-value thresholds 10^-1)
# MS3 (methylation score calculated with weights at p-value thresholds 10^-3)
# MS5 (methylation score calculated with weights at p-value thresholds 10^-5)
# smoking
# BMI
# Neut_Perc (neutrophils percentage)
# Mono_Perc (monocytes percentage)
# Eos_Perc (eosinophils percentage)
# Array_rownum (array rownumber)
# Sample_Plate (sample plate)
# Sample_Plate1 - Sample_PlateN (dummy variables for sample plate)

########################### ############## ##############
# Logistic regression to calculate R2 for PGS and MS
############################ ############## ##############
# function for logistic regression model with PGS and covariates prints PGS effect size, se, p-value, R2
PGSlogit_function = function(NTRdata) {
                  dat = NTRdata
                  hand01=as.numeric(dat$hand01)
                  sex=as.numeric(dat$sex)
                  PGS=as.numeric(dat$PGS)
                  PC1 = as.numeric(dat$PC1)
                  PC2 = as.numeric(dat$PC2)
                  PC3 = as.numeric(dat$PC3)
                  PC4 = as.numeric(dat$PC4)
                  PC5 = as.numeric(dat$PC5)
                  PC6 = as.numeric(dat$PC6)
                  PC7 = as.numeric(dat$PC7)
                  PC8 = as.numeric(dat$PC8)
                  PC9 = as.numeric(dat$PC9)
                  PC10 = as.numeric(dat$PC10)
                  platform2=as.numeric(dat$platform2)
                  platform3=as.numeric(dat$platform3)
                  platform5=as.numeric(dat$platform5)
                  platform6=as.numeric(dat$platform6)
                  # MODEL 1: handedness ~ PGS + GWAS covariates
                  vare1 = pi^2 / 3     # residual (homoskedastic) variance
                  r1=glm(hand01~PGS+sex+ platform2+platform3+platform5+platform6+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                         family=binomial(link='logit'))  # logistic regression
                  vp1=var( r1$coefficients["PGS"]*PGS+
                          r1$coefficients["sex"]*sex+
                          r1$coefficients["platform2"]*platform2+
                          r1$coefficients["platform3"]*platform3+
                          r1$coefficients["platform5"]*platform5+
                          r1$coefficients["platform6"]*platform6+
                          r1$coefficients["PC1"]*PC1+
                          r1$coefficients["PC2"]*PC2+
                          r1$coefficients["PC3"]*PC3+
                          r1$coefficients["PC4"]*PC4+
                          r1$coefficients["PC5"]*PC5+
                          r1$coefficients["PC6"]*PC6+
                          r1$coefficients["PC7"]*PC7+
                          r1$coefficients["PC8"]*PC8+
                          r1$coefficients["PC9"]*PC9+
                          r1$coefficients["PC10"]*PC10)    # variance explained
                  R2logist1 = vp1 / (vp1+vare1)   # explained var  / total var
                  # MODEL 2: handedness ~ GWAS covariates
                  r2=glm(hand01~sex+platform2+platform3+platform5+platform6+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                         family=binomial(link='logit'))
                  vp2=var(
                            r2$coefficients["sex"]*sex+
                            r2$coefficients["platform2"]*platform2+
                            r2$coefficients["platform3"]*platform3+
                            r2$coefficients["platform5"]*platform5+
                            r2$coefficients["platform6"]*platform6+
                            r1$coefficients["PC1"]*PC1+
                            r1$coefficients["PC2"]*PC2+
                            r1$coefficients["PC3"]*PC3+
                            r1$coefficients["PC4"]*PC4+
                            r1$coefficients["PC5"]*PC5+
                            r1$coefficients["PC6"]*PC6+
                            r1$coefficients["PC7"]*PC7+
                            r1$coefficients["PC8"]*PC8+
                            r1$coefficients["PC9"]*PC9+
                            r1$coefficients["PC10"]*PC10)    # variance explained
                  R2logist2 = vp2 / (vp2+vare1)   # explained var  / total var
                  t = matrix(0, 1, 6)
                  t[1,1] = summary(r1)$coefficients[2,1]
                  t[1,2] = summary(r1)$coefficients[2,2]
                  t[1,3] = summary(r1)$coefficients[2,4]
                  t[1,4] = R2logist1 - R2logist2
                  print(t)
}


# function for logistic regression model with MS, PGS and covariates prints MS effect size, se, p-value, R2
MSlogit_function = function(NTRdata, MScore) {
                  dat = NTRdata
                  MS = as.numeric(MScore)
                    hand01=as.numeric(dat$hand01)
                    age=as.numeric(dat$age)
                    sex=as.numeric(dat$sex)
                    PGS=as.numeric(dat$PGS)
                    PC1 = as.numeric(dat$PC1)
                    PC2 = as.numeric(dat$PC2)
                    PC3 = as.numeric(dat$PC3)
                    PC4 = as.numeric(dat$PC4)
                    PC5 = as.numeric(dat$PC5)
                    PC6 = as.numeric(dat$PC6)
                    PC7 = as.numeric(dat$PC7)
                    PC8 = as.numeric(dat$PC8)
                    PC9 = as.numeric(dat$PC9)
                    PC10 = as.numeric(dat$PC10)
                    platform1=as.numeric(dat$platform1) # not used in models
                    platform2=as.numeric(dat$platform2)
                    platform3=as.numeric(dat$platform3)
                    platform5=as.numeric(dat$platform5)
                    platform6=as.numeric(dat$platform6)
                    smoking = as.numeric(dat$smoking)
                    BMI = as.numeric(dat$BMI)
                    Neut_Perc = as.numeric(dat$Neut_Perc)
                    Mono_Perc = as.numeric(dat$Mono_Perc)
                    Eos_Perc = as.numeric(dat$Eos_Perc)
                    Sample_Plate = as.factor(dat$Sample_Plate)
                    Sample_Plate1 = ifelse(Sample_Plate== 1, 1, 0)
                    Sample_Plate2 = ifelse(Sample_Plate== 2, 1, 0)
                    Sample_Plate3 = ifelse(Sample_Plate== 3, 1, 0)
                    Sample_Plate4 = ifelse(Sample_Plate== 4, 1, 0)
                    Sample_Plate5 = ifelse(Sample_Plate== 5, 1, 0)
                    Sample_Plate6 = ifelse(Sample_Plate== 6, 1, 0)
                    Sample_Plate7 = ifelse(Sample_Plate== 7, 1, 0)
                    Sample_Plate8 = ifelse(Sample_Plate== 8, 1, 0)
                    Sample_Plate9 = ifelse(Sample_Plate== 9, 1, 0)
                    Sample_Plate10 = ifelse(Sample_Plate== 10, 1, 0)
                    Sample_Plate11 = ifelse(Sample_Plate== 11, 1, 0)
                    Sample_Plate12 = ifelse(Sample_Plate== 12, 1, 0)
                    Sample_Plate13 = ifelse(Sample_Plate== 13, 1, 0)
                    Sample_Plate14 = ifelse(Sample_Plate== 14, 1, 0)
                    Sample_Plate15 = ifelse(Sample_Plate== 15, 1, 0)
                    Sample_Plate16 = ifelse(Sample_Plate== 16, 1, 0)
                    Sample_Plate17 = ifelse(Sample_Plate== 17, 1, 0)
                    Sample_Plate18 = ifelse(Sample_Plate== 18, 1, 0)
                    Sample_Plate19 = ifelse(Sample_Plate== 19, 1, 0)
                    Sample_Plate20 = ifelse(Sample_Plate== 20, 1, 0)
                    Sample_Plate21 = ifelse(Sample_Plate== 21, 1, 0)
                    Sample_Plate22 = ifelse(Sample_Plate== 22, 1, 0)
                    Sample_Plate23 = ifelse(Sample_Plate== 23, 1, 0)
                    Sample_Plate24 = ifelse(Sample_Plate== 24, 1, 0)
                    Sample_Plate25 = ifelse(Sample_Plate== 25, 1, 0)
                    Sample_Plate26 = ifelse(Sample_Plate== 26, 1, 0)
                    Sample_Plate27 = ifelse(Sample_Plate== 27, 1, 0)
                    Sample_Plate28 = ifelse(Sample_Plate== 28, 1, 0)
                    Sample_Plate29 = ifelse(Sample_Plate== 29, 1, 0)
                    Sample_Plate30 = ifelse(Sample_Plate== 30, 1, 0)
                    Sample_Plate31 = ifelse(Sample_Plate== 31, 1, 0)
                    Sample_Plate32 = ifelse(Sample_Plate== 32, 1, 0)
                    Sample_Plate33 = ifelse(Sample_Plate== 33, 1, 0)
                    Sample_Plate34 = ifelse(Sample_Plate== 34, 1, 0)
                    # MODEL 1: handedness ~ MS+PGS + covariates
                    vare1 = pi^2 / 3     # residual (homoskedastic) variance
                    r1=glm(hand01~PGS+age+sex+
                             platform2+platform3+platform5+platform6+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                             smoking+BMI+Neut_Perc+Mono_Perc+Eos_Perc+Sample_Plate+MS,
                             family=binomial(link='logit'))  # logistic regression
                    vp1=var( r1$coefficients["PGS"]*PGS+
                            r1$coefficients["age"]*age+
                            r1$coefficients["sex"]*sex+
                            r1$coefficients["platform2"]*platform2+
                            r1$coefficients["platform3"]*platform3+
                            r1$coefficients["platform5"]*platform5+
                            r1$coefficients["platform6"]*platform6+
                            r1$coefficients["PC1"]*PC1+
                            r1$coefficients["PC2"]*PC2+
                            r1$coefficients["PC3"]*PC3+
                            r1$coefficients["PC4"]*PC4+
                            r1$coefficients["PC5"]*PC5+
                            r1$coefficients["PC6"]*PC6+
                            r1$coefficients["PC7"]*PC7+
                            r1$coefficients["PC8"]*PC8+
                            r1$coefficients["PC9"]*PC9+
                            r1$coefficients["PC10"]*PC10+
                            r1$coefficients["smoking"]*smoking+
                            r1$coefficients["BMI"]*BMI+
                            r1$coefficients["Neut_Perc"]*Neut_Perc+
                            r1$coefficients["Mono_Perc"]*Mono_Perc+
                            r1$coefficients["Eos_Perc"]*Eos_Perc+
                            r1$coefficients["Sample_Plate2"]*Sample_Plate2+
                            r1$coefficients["Sample_Plate3"]*Sample_Plate3+
                            r1$coefficients["Sample_Plate4"]*Sample_Plate4+
                            r1$coefficients["Sample_Plate5"]*Sample_Plate5+
                            r1$coefficients["Sample_Plate6"]*Sample_Plate6+
                            r1$coefficients["Sample_Plate7"]*Sample_Plate7+
                            r1$coefficients["Sample_Plate8"]*Sample_Plate8+
                            r1$coefficients["Sample_Plate9"]*Sample_Plate9+
                            r1$coefficients["Sample_Plate10"]*Sample_Plate10+
                            r1$coefficients["Sample_Plate11"]*Sample_Plate11+
                            r1$coefficients["Sample_Plate12"]*Sample_Plate12+
                            r1$coefficients["Sample_Plate13"]*Sample_Plate13+
                            r1$coefficients["Sample_Plate14"]*Sample_Plate14+
                            r1$coefficients["Sample_Plate15"]*Sample_Plate15+
                            r1$coefficients["Sample_Plate16"]*Sample_Plate16+
                            r1$coefficients["Sample_Plate17"]*Sample_Plate17+
                            r1$coefficients["Sample_Plate18"]*Sample_Plate18+
                            r1$coefficients["Sample_Plate19"]*Sample_Plate19+
                            r1$coefficients["Sample_Plate20"]*Sample_Plate20+
                            r1$coefficients["Sample_Plate21"]*Sample_Plate21+
                            r1$coefficients["Sample_Plate22"]*Sample_Plate22+
                            r1$coefficients["Sample_Plate23"]*Sample_Plate23+
                            r1$coefficients["Sample_Plate24"]*Sample_Plate24+
                            r1$coefficients["Sample_Plate25"]*Sample_Plate25+
                            r1$coefficients["Sample_Plate26"]*Sample_Plate26+
                            r1$coefficients["Sample_Plate27"]*Sample_Plate27+
                            r1$coefficients["Sample_Plate28"]*Sample_Plate28+
                            r1$coefficients["Sample_Plate29"]*Sample_Plate29+
                            r1$coefficients["Sample_Plate30"]*Sample_Plate30+
                            r1$coefficients["Sample_Plate31"]*Sample_Plate31+
                            r1$coefficients["Sample_Plate32"]*Sample_Plate32+
                            r1$coefficients["Sample_Plate33"]*Sample_Plate33+
                            r1$coefficients["Sample_Plate34"]*Sample_Plate34+
                            r1$coefficients["MS"]*MS)   # variance explained
                    R2logist1 = vp1 / (vp1+vare1)   # explained var  / total var
                    # MODEL 2: handedness ~ PGS +  covariates
                    r2=glm(hand01~PGS5+age+sex+platform2+platform3+platform5+platform6+
                                PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                                smoking+BMI+Neut_Perc+Mono_Perc+Eos_Perc+Sample_Plate,
                                family=binomial(link='logit'))
                    
                    vp2=var(            r2$coefficients["PGS"]*PGS+
                                        r2$coefficients["age"]*age+
                                        r2$coefficients["sex"]*sex+
                                        r2$coefficients["platform2"]*platform2+
                                        r2$coefficients["platform3"]*platform3+
                                        r2$coefficients["platform5"]*platform5+
                                        r2$coefficients["platform6"]*platform6+
                                        r2$coefficients["PC1"]*PC1+
                                        r2$coefficients["PC2"]*PC2+
                                        r2$coefficients["PC3"]*PC3+
                                        r2$coefficients["PC4"]*PC4+
                                        r2$coefficients["PC5"]*PC5+
                                        r2$coefficients["PC6"]*PC6+
                                        r2$coefficients["PC7"]*PC7+
                                        r2$coefficients["PC8"]*PC8+
                                        r2$coefficients["PC9"]*PC9+
                                        r2$coefficients["PC10"]*PC10+
                                        r2$coefficients["smoking"]*smoking+
                                        r2$coefficients["BMI"]*BMI+
                                        r2$coefficients["Neut_Perc"]*Neut_Perc+
                                        r2$coefficients["Mono_Perc"]*Mono_Perc+
                                        r2$coefficients["Eos_Perc"]*Eos_Perc+
                                        r2$coefficients["Sample_Plate2"]*Sample_Plate2+
                                        r2$coefficients["Sample_Plate3"]*Sample_Plate3+
                                        r2$coefficients["Sample_Plate4"]*Sample_Plate4+
                                        r2$coefficients["Sample_Plate5"]*Sample_Plate5+
                                        r2$coefficients["Sample_Plate6"]*Sample_Plate6+
                                        r2$coefficients["Sample_Plate7"]*Sample_Plate7+
                                        r2$coefficients["Sample_Plate8"]*Sample_Plate8+
                                        r2$coefficients["Sample_Plate9"]*Sample_Plate9+
                                        r2$coefficients["Sample_Plate10"]*Sample_Plate10+
                                        r2$coefficients["Sample_Plate11"]*Sample_Plate11+
                                        r2$coefficients["Sample_Plate12"]*Sample_Plate12+
                                        r2$coefficients["Sample_Plate13"]*Sample_Plate13+
                                        r2$coefficients["Sample_Plate14"]*Sample_Plate14+
                                        r2$coefficients["Sample_Plate15"]*Sample_Plate15+
                                        r2$coefficients["Sample_Plate16"]*Sample_Plate16+
                                        r2$coefficients["Sample_Plate17"]*Sample_Plate17+
                                        r2$coefficients["Sample_Plate18"]*Sample_Plate18+
                                        r2$coefficients["Sample_Plate19"]*Sample_Plate19+
                                        r2$coefficients["Sample_Plate20"]*Sample_Plate20+
                                        r2$coefficients["Sample_Plate21"]*Sample_Plate21+
                                        r2$coefficients["Sample_Plate22"]*Sample_Plate22+
                                        r2$coefficients["Sample_Plate23"]*Sample_Plate23+
                                        r2$coefficients["Sample_Plate24"]*Sample_Plate24+
                                        r2$coefficients["Sample_Plate25"]*Sample_Plate25+
                                        r2$coefficients["Sample_Plate26"]*Sample_Plate26+
                                        r2$coefficients["Sample_Plate27"]*Sample_Plate27+
                                        r2$coefficients["Sample_Plate28"]*Sample_Plate28+
                                        r2$coefficients["Sample_Plate29"]*Sample_Plate29+
                                        r2$coefficients["Sample_Plate30"]*Sample_Plate30+
                                        r2$coefficients["Sample_Plate31"]*Sample_Plate31+
                                        r2$coefficients["Sample_Plate32"]*Sample_Plate32+
                                        r2$coefficients["Sample_Plate33"]*Sample_Plate33+
                                        r2$coefficients["Sample_Plate34"]*Sample_Plate34
                    )    # variance explained
                    R2logist2 = vp2 / (vp2+vare1)   # explained var  / total var
                    t = matrix(0, 1, 6)
                    t[1,1] = summary(r1)$coefficients[57,1]
                    t[1,2] = summary(r1)$coefficients[57,2]
                    t[1,3] = summary(r1)$coefficients[57,4]
                    t[1,6] = R2logist1 - R2logist2
                    print(t)
}


