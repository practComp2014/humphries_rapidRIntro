#============================================
# 
#  Practical Computing: Working with data.table.
#  John Eric Humphries
#  2014-04-17
#
#
#============================================

#========================
# Section 0: setup
#========================

#setwd("/mnt/ide0/home/johneric/sbox/projects/neighborhoodVision/")
rm(list=ls())           # Clear the workspace
#set.seed(907) 
library(data.table)
#library(plyr)
library(stargazer)
library(texreg)
library(xtable)
library(AER)

# This draws heavily  from the tutorial provided by Mathew Dowle here: http://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.pdf
# Talk about usefulness of dataframes in R compared to datasets in stata, etc. 

# A data frame
DF = data.frame(x=c("b","b","b","a","a"),v=rnorm(5))
DF

# A data table
DT = data.table(x=c("b","b","b","a","a"),v=rnorm(5))
DT

# A pre-loaded data set:
data(USConsump1993)
CONS =  data.frame(USConsump1993)
CONS <- data.table(CONS)
head(CONS)

# Tell me what tables I have in my environment (also shown in R-studio panel)
tables() # note that date.tables are held in ram!

# info on columns
sapply(DT,class)
sapply(CONS,class)

# Data tables can be manipulated like data frames for the most part
DF[2,]
DT[2,]
DT[DT$x=="a"]

#-------------------
# Keys in data.tables
#-------------------

# much like databases, each table can have a column set as the key. The key will be used for sorting (though it need not be unique).

setkey(DT,x)
DT  # note its now sorted by x.
haskey(DT)
key(DT)
tables()

# can now reference by key-value
DT["b",] 
DT["b"] 
DT["b", mult = "first"]
DT["b", mult = "last"]

# Now lets make a really big data.frame:
grpsize = ceiling(1e7/26^2)
system.time( DF <- data.frame( x=rep(LETTERS,each=26*grpsize),y=rep(letters,each=grpsize),v=runif(grpsize*26^2),w=rnorm(grpsize*26^2),stringsAsFactors=FALSE))
dim(DF)
head(DF,5)

system.time(
 ans1 <- DF[DF$x=="R" & DF$y == "h",]    
)
head(ans1,5)

# Now lets try it out as a data table
DT <- as.data.table(DF)
system.time( setkey(DT,x,y))

system.time(
    ans2 <-DT[data.table("R","h")]    # behold the power of binary search! (DF uses vector scan.)
)

system.time(
    ans2 <-DT[J("R","h")]    # J works as an alias to avoid wrting data.table all the time...
)
# 900 times faster on my work station.
head(ans2,5)

# can use data.frame syntax still (but will go back to our previous slow life... but even here not AS slow...)
system.time(
    ans1 <- DT[x=="R" & y=="h",]
)

#========================
# Behind the comma in data.tables!
#========================

# After the comma we can use expressions acting on col-names like they are variables:
DT[,sum(v)]
DT[,sum(v), by=x] # can also easily do by statements


# once again, this code is WAY faster tapply or doing things in dataframees...
system.time(
 tt <- tapply(DT$v,DT$x,sum)    
)

system.time(
 ss <- DT[,sum(v), by=x]    
)

head(tt)
head(ss)

# Grouping by two columns:
system.time(tt <- tapply(DT$v,list(DT$x,DT$y),sum))
system.time(ss <- DT[,sum(v),by="x,y"])

# can return serveral expressions
DT[,c(sum(v),min(v))]

# There are special things about data.tables that help them be really useful/fast for time series, i just don't know much about them. 


# The with command
with(DT,v + w)


# Data tables lookups return data tables, so we can have recursive calls:
#DT[where,select,by=...][order(...)][...][...]...
DT["A",sum(v),by=y]
DT[,"w"] <- NULL

# fread() is much better than read.table or read.csv

#=============================
# Working with some real data. 
#=============================

# Reading in 1988 CPS data
data(CPS1988)
CPS <- data.table(CPS1988)
names(CPS)
head(CPS)

CPS[,table(education)]
CPS[,length(education), by=education]

#deleting an entry
CPS[,smsa:=NULL]

# adding an entry
CPS[,newcol:=mean(wage),by=region]
CPS[,unique(newcol),by=region]


mincer1 <- lm(wage ~ education + experience + I(experience^2), data=CPS)
mincer2 <- lm(wage ~ education + experience + I(experience^2) + region, data=CPS)

xtable(mincer1) # see examples here: http://cran.r-project.org/web/packages/xtable/vignettes/xtableGallery.pdf
stargazer(mincer1, mincer2) # http://cran.r-project.org/web/packages/stargazer/vignettes/stargazer.pdf
texreg(list(mincer1, mincer2)) # http://www.jstatsoft.org/v55/i08/paper
# Many other ways to output latex: http://stackoverflow.com/questions/5465314/tools-for-making-latex-tables-in-r

wagesummary =CPS[,
    list(count=.N, mean_wage =mean(wage), sd_wage =sd(wage), max_wage = max(wage), min_wage = min(wage)),
    by=list(education, ethnicity)]
setkey(wagesummary,education,ethnicity)
wagesummary
xtable(wagesummary)
stargazer(CPS)

dt1 <- data.table(cbind(1:100,LETTERS))
setkey(dt1,V1)
dt2 <- data.table(cbind(1:100,letters))
setkey(dt2,V1)

dt1[dt2]