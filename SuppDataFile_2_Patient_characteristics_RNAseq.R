#Table S5 - Patient characteristics in RNAseq subcohort

#Clear R environment
rm(list = ls())

#Import the master dataset of all combined patients, PID to fixed row ID
combined <- read.csv(file.choose(), header=TRUE)
rownames(combined) = combined$pid
combined$pid = NULL

#Select covid patients with RNAseq data
covid <- subset(combined, pathogencode==1 & rna_seq==1)

#Table of clinical characteristics for RNAseq subcohort - TableOne package
library(tableone)
#Create a variable list which we want in Table
listVars <- c("sex", 
              "age",
              "illnessdurationenroll_imp",
              "phasedelta",
              "historyoffever", 
              "nightsweats", 
              "headache", 
              "cough", 
              "sorethroat",
              "runnynose",
              "sob",
              "tempmax_imp", 
              "heartrate3_imp", 
              "resprate3_imp", 
              "sbp3_imp",
              "o2sat3", 
              "kpsadmit16plus_imp",
              "lactate_imp",
              "hgb_imp",
              "wbc_imp",
              "platelets_imp",
              "malariardtresult",
              "hivrdtresult", 
              "artprior",
              "priortuberculosis",
              "heartdisease",
              "hypertension",
              "diabetes",
              "whoseverity_imp",
              "supplementaloxygen",
              "supplementaloxygenamount",
              "abx",
              "steroids",
              "hospdeathtransf") 

#Define categorical variables
catVars <- c("sex", 
             "phasedelta",
             "historyoffever", 
             "nightsweats", 
             "headache", 
             "cough", 
             "sorethroat",
             "runnynose",
             "sob",
             "malariardtresult",
             "hivrdtresult", 
             "artprior",
             "priortuberculosis",
             "heartdisease",
             "hypertension",
             "diabetes",
             "whoseverity_imp",
             "supplementaloxygen",
             "abx",
             "steroids",
             "hospdeathtransf") 

tablernaseq <- CreateTableOne(vars = listVars, 
                                        data = covid, 
                                        factorVars = catVars,
                                        includeNA = TRUE,
                                        addOverall = TRUE)
tablernaseq
summary(tablernaseq)

#Now specify we need medians(IQR) for continuous variables and Fisher exact for small cell counts
tablernaseq <- print(tablernaseq, nonnormal = c("age",
                                                                    "tempmax_imp",
                                                                    "illnessdurationenroll_imp",
                                                                    "heartrate3_imp", 
                                                                    "resprate3_imp", 
                                                                    "sbp3_imp",
                                                                    "o2sat3",
                                                                    "lactate_imp",
                                                                    "hgb_imp",
                                                                    "wbc_imp",
                                                                    "platelets_imp",
                                                                    "kpsadmit16plus_imp",
                                                                    "supplementaloxygenamount"))

#Now export table to HTML
library(tableHTML)
write_tableHTML(tableHTML(tablernaseq), file = 'tablecovidrnaseq.html')