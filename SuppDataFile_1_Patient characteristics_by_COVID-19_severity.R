#Table 1 - Patient characteristics stratified by COVID-19 severity 

#Clear R environment
rm(list = ls())

#Import the master dataset of all combined patients, PID to fixed row ID
combined <- read.csv(file.choose(), header=TRUE)
rownames(combined) = combined$pid
combined$pid = NULL

#Select covid patients
covid <- subset(combined, pathogencode==1 )

#Table comparing clinical characteristics by COVID severity - TableOne package
library(tableone)
#Create a variable list which we want in Table
listVars <- c("sex", 
              "age",
              "illnessdurationenroll_imp",
              "phase3",
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
              "hiv_suppressed",
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
             "phase3",
              "historyoffever", 
              "nightsweats", 
              "headache", 
              "cough", 
              "sorethroat",
              "runnynose",
              "sob",
              "malariardtresult",
              "hivrdtresult", 
              "hiv_suppressed",
              "artprior",
              "heartdisease",
              "hypertension",
              "diabetes",
              "priortuberculosis",
              "whoseverity_imp",
              "supplementaloxygen",
              "abx",
              "steroids",
              "hospdeathtransf") 

tableclinicalseverity <- CreateTableOne(vars = listVars, 
                                        data = covid, 
                                        factorVars = catVars,
                                        strata = "whoseverity_imp",
                                        includeNA = TRUE,
                                        addOverall = TRUE)
tableclinicalseverity
summary(tableclinicalseverity)

#Now specify we need medians(IQR) for continuous variables and Fisher exact for small cell counts
tableclinicalseverity <- print(tableclinicalseverity, nonnormal = c("age",
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
write_tableHTML(tableHTML(tableclinicalseverity), file = 'tablecovidclinicalseverity.html')