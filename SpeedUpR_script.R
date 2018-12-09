#=================================================
#
# How to speed up your R code?
#
# Josephine Daub
# R-Ladies meetup December 11, 2018
#
# Some examples taken from: 
# https://www.r-bloggers.com/strategies-to-speedup-r-code/ 
#
# See also: 
# http://adv-r.had.co.nz/Profiling.html
# https://www.alexejgossmann.com/benchmarking_r/
# https://www.r-bloggers.com/using-the-microbenchmark-package-to-compare-the-execution-time-of-r-expressions/
# https://csgillespie.github.io/efficientR/programming.html
#
#=================================================

#=================================================
# TO PREPARE
#=================================================

# Install packages and test if they can be loaded

# install.packages("reshape2")
# install.packages("microbenchmark")
# install.packages("RColorBrewer")
# install.packages("Rcpp")
# install.packages("data.table")

# Set working directory to location of scripts
#setwd(<location of SpeedUpR scripts>)


#=================================================
# Ex 1: Count co-occurrences in gene-sample matrix
#=================================================

source("SpeedUpR_functions.R")

# Create binary gene-sample data frame 
res<-CreateGeneSampleDF(ng=9, ns=5, prob1=0.4)
df<-res$df
df

# count per gene-pair the number of co-occurences
df.count<-GetCO_0(df)
df.count
df.orig<-df.count

#------------------------------------------------
# Timing
#------------------------------------------------
res<-CreateGeneSampleDF(ng=100, ns=25, prob1=0.4)
df<-res$df

# how long does my function take?
system.time({
  df.count<-GetCO_0(df)
})
# store in df.orig to test if other solutions
# give similar result

df.orig<-df.count

#------------------------------------------------
# Profiling
#------------------------------------------------
# start profiling to this file
Rprof("Rprof_COcount.txt")
df.count<-GetCO_0(df)
Rprof(NULL) # stop profiling
summaryRprof("Rprof_COcount.txt")

# Profiling using line profiling
Rprof("Rprof_COcount.txt", line.profiling = T) 
df.count<-GetCO_0(df)
Rprof(NULL) 
summaryRprof("Rprof_COcount.txt", lines="show")

# Profiling using sub functions
Rprof("Rprof_COcount.txt") 
df.count<-GetCO_0_func(df)
Rprof(NULL) 
summaryRprof("Rprof_COcount.txt")


#------------------------------------------------
# allocate memory: initiate df.count
#------------------------------------------------
system.time({
  df.count<-GetCO_1_alloc_mem(df)
})
# test if df.count=df.orig
identical(df.count,df.orig)

# Check with profiling the difference
Rprof("Rprof_COcount.txt", line.profiling = T) 
df.count<-GetCO_1_alloc_mem(df)
Rprof(NULL) 
summaryRprof("Rprof_COcount.txt", lines="show")


#------------------------------------------------
# use matrix instead of df
#------------------------------------------------
res<-CreateGeneSampleDF(ng=dim(df)[1], 
                        ns=dim(df)[2], prob1=0.4)
mtx<-res$mtx

mtx[1:5,1:4]
df[1:5,1:4]
str(mtx)
str(df)

system.time({
  df.count<-GetCO_0(mtx)
})
# test if df.count=df.orig
identical(df.count,df.orig)


#------------------------------------------------
# reserve memory: initiate df.count with mtx
#------------------------------------------------
system.time({
  df.count<-GetCO_1_alloc_mem(mtx)
})
# test if df.count=df.orig
identical(df.count,df.orig)


#------------------------------------------------
# use vectorization
#------------------------------------------------
g1<-c(1,1,0,1,0)
g1
5*g1
g2<-c(0,1,0,1,0)
g1
g2
g1*g2
sum(g1*g2)

system.time({
  df.count<-GetCO_2_vect(mtx)
})
# test if df.count=df.orig
identical(df.count,df.orig)


#------------------------------------------------
# use matrix manipulation
#------------------------------------------------
res<-CreateGeneSampleDF(ng=6, ns=3, prob1=0.4)
mtx.small<-res$mtx
mtx.small
df.count.small<-GetCO_2_vect(mtx.small)
df.count.small
mtx.small %*% t(mtx.small)

require(reshape2)
system.time({
  df.count<-GetCO_3_mtx(mtx)
})
identical(df.count,df.orig)
# sometimes faster: function crossprod/tcrossprod




#------------------------------------------------
# Benchmarking for different matrix sizes
#------------------------------------------------
source("SpeedUpR_functions.R")
library(microbenchmark)
library(reshape2)
library(RColorBrewer)

# vbenchmark = number of samples
# number of genes: number of samples * 4
vbenchmark<-c(5,10,15,20)
# the next one might take some time
#vbenchmark<-c(5,10,15,20,25,30,40,50,60,70)
res<-BenchMark_GetCO(vbenchmark, times=1, maxsec=50,
                     run=T, plot=T)
df.benchmark<-res$df.benchmark
nc.lims<-res$nc.lims




#=================================================
# Ex 2: Which rows have sum > 4 ?
# Examples taken from: 
# https://www.r-bloggers.com/strategies-to-speedup-r-code/ 
#=================================================

source("SpeedUpR_functions.R")

# create data frame
df<-CreateBigDF(n=10)
df

#------------------------------------------------
# use for loop
#------------------------------------------------

# naive solution
N <- nrow(df)
output<-character()
for (i in 1:N){
  if (sum(df[i,])>4) {
    output[i]<-"greater"
  } else {
    output[i]<-"less_or_equal"  
  }
}

output

#------------------------------------------------
# profiling
#------------------------------------------------

# create big data frame
df<-CreateBigDF(n=1e4)
dim(df)
system.time({
  output<-Row4_naive(df)
})

Rprof("Rprof_Row4.txt", line.profiling = T) 
output<-Row4_naive(df)
Rprof(NULL)
summaryRprof("Rprof_Row4.txt", lines="show")



# create bigger data frame
df<-CreateBigDF(n=1e5)
system.time({
  output<-Row4_naive(df)
})
output.orig<-output


#------------------------------------------------
# remove condition out of loop
#------------------------------------------------
df.S<-CreateBigDF(n=5)
df.S

# use vectorization
df.S$col1 + df.S$col2 + df.S$col3 + df.S$col4
(df.S$col1 + df.S$col2 + df.S$col3 + df.S$col4)>4

# or: use build-in functions
# rowSums, colSums, rowMeans, etc
rowSums(df.S)
rowSums(df.S)>4

system.time({
  N <- nrow(df)
  output<-character()
  precheck<-rowSums(df)>4
  for (i in 1:N){
    if (precheck[i]) {
      output[i]<-"greater"
    } else {
      output[i]<-"less_or_equal"  
    }
  }
})
identical(output,output.orig)


#------------------------------------------------
# use memory allocation and pre-check
#------------------------------------------------
system.time({
  N <- nrow(df)
  output <- character(N)
  precheck<-rowSums(df)>4
  for (i in 1:N){
    if (precheck[i]) {
      output[i]<-"greater"
    } else {
      output[i]<-"less_or_equal"  
    }
  }
})
identical(output,output.orig)

#------------------------------------------------
# use default value
#------------------------------------------------
system.time({
  N <- nrow(df)
  output <- rep("less_or_equal",N)
  precheck<-rowSums(df)>4
  for (i in 1:N){
    if (precheck[i]) {
      output[i]<-"greater"
    } 
  }
})
identical(output,output.orig)


#------------------------------------------------
# use sapply without precheck
#------------------------------------------------
system.time({
  N <- nrow(df)
  output<-sapply(1:N, FUN=function(x){
    if (sum(df[x,])>4) return("greater")
    else return("less_or_equal")
    })
})
identical(output,output.orig)


#------------------------------------------------
# use sapply with precheck
#------------------------------------------------
system.time({
  N <- nrow(df)
  precheck<-rowSums(df)>4
  output<-sapply(1:N, FUN=function(x){
    if (precheck[x]) return("greater")
    else return("less_or_equal")
  })
})
identical(output,output.orig)


#------------------------------------------------
# use apply
#------------------------------------------------
system.time({
  output<-apply(df, 1, FUN=function(x){
    if (sum(x)>4) return("greater")
    else return("less_or_equal")
  })
})
identical(output,output.orig)

# create bigger data frame
df<-CreateBigDF(n=1e7)


#------------------------------------------------
# use Vectorize
#------------------------------------------------
CheckFour <- function(v1,v2,v3,v4){
  if (sum(c(v1,v2,v3,v4))>4) {
    return("greater") 
  } else {
    return("less_or_equal")
  }
}
vCheckFour <-Vectorize(CheckFour, SIMPLIFY = TRUE)

system.time({
  output<-vCheckFour(df$col1,df$col2,df$col3,df$col4)
})
identical(output,output.orig)


#------------------------------------------------
# use ifelse and rowSums
#------------------------------------------------
system.time({
  output<-ifelse(rowSums(df)>4,"greater",
                 "less_or_equal")
})
identical(output,output.orig)


#------------------------------------------------
# use which
#------------------------------------------------
system.time({
  N <- nrow(df)
  output <- rep("less_or_equal",N)
  output[which(rowSums(df)>4)]<-"greater"
  #comparably fast:
  #output[rowSums(df)>4]<-"greater"
})
identical(output,output.orig)


#------------------------------------------------
# use Rcpp
#------------------------------------------------

library(Rcpp)
sourceCpp("SpeedUpR.cpp")
system.time({output <- cppCheckFour(df)})
identical(output,output.orig)




#------------------------------------------------
# benchmarking
#------------------------------------------------
source("SpeedUpR_functions.R")
require(microbenchmark)
require(Rcpp)
sourceCpp("SpeedUpR.cpp")

funcs<-c("naive", "mem","mem_def","sapply",
         "precheck","precheck_mem", 
         "precheck_mem_def", "sapply_precheck", 
         "vect", "apply", "ifelse", "which", "rcpp")
vbenchmark<-c(1e3,5e3,1e4)
# this one will take some time
#vbenchmark<-c(1e3,5e3,1e4,2e4,5e4,8e4,1e5,2e5,5e5,8e5,
#              1e6,2e6,5e6,8e6,1e7,2e7,3e7,5e7)

txt<-"test"
#txt<-"R33_windows"
#txt<-"R34_mac"
res<-BenchMark_Row4(vbenchmark, funcs, times=3, 
                  maxsec=10, run=F, plot=T, 
                  txt=txt)


#=================================================
# Ex 3: How to avoid aggregate
#=================================================

# how to run function per grouping?
# aggregate, tapply, data.table, home made function

source("SpeedUpR_functions.R")

# create big matrix
df<-CreateBigDF(n=1e07)
head(df)
# add groupings
df$group<-sample(1:4,nrow(df), replace = TRUE)
head(df)

#------------------------------------------------
# use aggregate function
#------------------------------------------------
system.time({
  df.agg<-aggregate(df$col1,
                     by=list(group=df$group),
                     FUN = mean)
  })
df.agg

#------------------------------------------------
# use tapply function
#------------------------------------------------
system.time({
  tab.agg<-tapply(df$col1,df$group,
                             FUN=mean)
})
tab.agg

#------------------------------------------------
# use home made function
#------------------------------------------------
system.time({
  df.agg.fast<-FastAggr_1col(df,col="group",
                             fcol="col1", 
                             func="mean", 
                             ncol="mean1")
  
})
df.agg.fast

#------------------------------------------------
# use data table
#------------------------------------------------
require(data.table)
dt<-as.data.table(df)
system.time({
  dt.agg<-dt[, .(mean = mean(col1)), 
             by = group]
})
dt.agg

#------------------------------------------------
#------------------------------------------------

# how do we add another summary statistic?
# max[col4]

# create big matrix
df<-CreateBigDF(n=2e07)
# add groupings
df$group<-sample(1:4,nrow(df), replace = TRUE)

#------------------------------------------------
# use aggregate function
#------------------------------------------------
system.time({
  df.agg1<-aggregate(df[,c("col1")],
                    by=list(group=df$group),
                    FUN = mean)
  
  df.agg2<-aggregate(df[,c("col4")],
                     by=list(group=df$group),
                     FUN = max)
  df.agg<-merge(df.agg1, df.agg2, by="group") 
  colnames(df.agg)<-c("group","mean1","max4")
})
df.agg

#------------------------------------------------
# use tapply
#------------------------------------------------
system.time({
  tab.agg1<-tapply(df$col1,df$group,
                  FUN=mean)
  tab.agg2<-tapply(df$col4,df$group,
                   FUN=max)
  tab.agg<-data.frame(group=names(tab.agg1),
                      mean1=tab.agg1,max4=tab.agg2)
})
tab.agg


#------------------------------------------------
# use home made function
#------------------------------------------------
system.time({
  df.agg.fast<-FastAggr_1col(df,col="group",
                             fcol=c("col1","col4"), 
                             func=c("mean","max"), 
                             ncol=c("mean1","max4"))
  
})
df.agg.fast

#------------------------------------------------
# use data.table
#------------------------------------------------
dt<-data.table(df)
system.time({
  dt.agg<-dt[, list(mean1=mean(col1), max4=max(col4)), 
             by=c("group")]
})
dt.agg
