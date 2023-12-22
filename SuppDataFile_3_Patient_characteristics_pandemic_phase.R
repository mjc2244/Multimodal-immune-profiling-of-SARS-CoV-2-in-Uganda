#Table S8 - Patient characteristics stratified by SARS-CoV-2-driven pandemic phase

#Clear R environment
rm(list = ls())

#Import the master dataset of all combined patients, row ID to fixed
combined <- read.csv(file.choose(), header=TRUE)
rownames(combined) = combined$pid
combined$pid = NULL

#Select covid
covid <- subset(combined, pathogencode==1 )

#Table comparing clinical characteristics across illness phase - TableOne package
library(tableone)
#Create a variable list which we want in Table
listVars <- c("sex", 
              "age",
              "illnessdurationenroll_imp",
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
              "supplementaloxygen",
              "supplementaloxygenamount",
              "abx",
              "steroids",
              "hospdeathtransf") 

#Define categorical variables
catVars <- c("sex", 
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
             "heartdisease",
             "hypertension",
             "diabetes",
             "priortuberculosis",
             "whophase_imp",
             "supplementaloxygen",
             "abx",
             "steroids",
             "hospdeathtransf") 

tableclinicalphase <- CreateTableOne(vars = listVars, 
                                        data = covid, 
                                        factorVars = catVars,
                                        strata = "phase3",
                                        includeNA = TRUE,
                                        addOverall = TRUE)
tableclinicalphase
summary(tableclinicalphase)

#Now specify we need medians(IQR) for continuous variables and Fisher exact for small cell counts
tableclinicalphase <- print(tableclinicalphase, nonnormal = c("age",
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
write_tableHTML(tableHTML(tableclinicalphase), file = 'tablecovidclinicalphase.html')