#Table S15 - Patient characteristics for non-COVID SARI cohort

#Clear R environment
rm(list = ls())

#Import the master dataset of all combined patients, row ID to fixed
combined <- read.csv(file.choose(), header=TRUE)
rownames(combined) = combined$pid
combined$pid = NULL

#Table for patient characteristics 
library(tableone)
#Create a variable list which we want in Table
listVars <- c("sex", 
              "age",
              "pathogen",
              "malariardtresult",
              "hivrdtresult", 
              "whoseverity_imp",
              "hospdeathtransf") 

#Define categorical variables
catVars <- c("sex", 
             "pathogen",
             "malariardtresult",
             "hivrdtresult", 
             "whoseverity_imp",
             "hospdeathtransf") 


tableclinicalsari <- CreateTableOne(vars = listVars, 
                                        data = subset(combined, covid==0),
                                        factorVars = catVars,
                                        includeNA = TRUE,
                                        addOverall = TRUE)
tableclinicalsari
summary(tableclinicalsari)

#Now specify we need medians(IQR) for continuous variables and Fisher exact for small cell counts
tableclinicalsari <- print(tableclinicalsari, nonnormal = c("age"))
                                                                    
#Now export table to HTML
library(tableHTML)
write_tableHTML(tableHTML(tableclinicalsari), file = 'tablecovidsariclinical.html')