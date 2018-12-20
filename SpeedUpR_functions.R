#=================================================
#
# How to speed up your R code? - FUNCTIONS
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
# Count co-occurrences in gene-sample matrix
#=================================================

CreateGeneSampleDF <- function (ng=10, ns=6, 
                                prob1=0.4, 
                                seed=1){
  set.seed(seed)
  mtx<-matrix(sample(c(0,1), ng*ns, replace = T, 
                     prob = c(1-prob1,prob1)), 
                     nrow = ng, ncol=ns, 
                     dimnames = 
                       list(paste0("gene",1:ng),
                            paste0("sample",1:ns)))
  df<-as.data.frame(mtx)
  return(list(mtx=mtx,df=df))
}

# naive function
GetCO_0 <- function(x){
  N <- nrow(x)
  M <- ncol(x)
  genes<-rownames(x)
  df.CO <- data.frame()
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      vg1<-x[i,]
      vg2<-x[j,]
      COcount<-0
      for (k in 1:M){
        if (vg1[k]==1 & vg2[k]==1){
          COcount <- COcount+1
        }
      }
      df.CO<-rbind(df.CO,
                   data.frame(g1=genes[i], 
                              g2=genes[j],
                              co=COcount, 
                              stringsAsFactors=F))
    }
  }
  return(df.CO)
}

#------------------------------------------------

# naive function with inner functions
GetCO_0_func <- function(x){
  N <- nrow(x)
  M <- ncol(x)
  genes<-rownames(x)
  df.CO <- data.frame()
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      vg1<-GetGeneRow(x,i)
      vg2<-GetGeneRow(x,j)
      COcount<-AddCoCount(M,vg1,vg2)
      df.CO<-AddGenePair(df.CO,g1=genes[i], 
                         g2=genes[j], co=COcount)
    }
  }
  return(df.CO)
}

GetGeneRow<-function(x, g){
  return(x[g,])
}

AddCoCount<-function(M,vg1,vg2){
  COcount<-0
  for (k in 1:M){
    if (vg1[k]==1 & vg2[k]==1){
      COcount <- COcount+1
    }
  }
  return(COcount)
}

AddGenePair<-function(df,g1,g2,co){
  return(rbind(df,data.frame(g1=g1,g2=g2,co=co,
                             stringsAsFactors=F)))
}

#------------------------------------------------

# pre-allocate memory: initiate df.CO in advance
GetCO_1_alloc_mem <- function(x){
  N <- nrow(x)
  M <- ncol(x)
  genes<-rownames(x)
  
  # create an empty data.frame in advance
  nb_pairs<-N*(N-1)/2
  df.CO <- data.frame(g1=rep("",nb_pairs),
                      g2=rep("",nb_pairs),
                      co=0,stringsAsFactors = F)
  p <- 0
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      vg1<-x[i,]
      vg2<-x[j,]
      COcount<-0
      for (k in 1:M){
        if (vg1[k]==1 & vg2[k]==1){
          COcount <- COcount+1
        }
      }
      p <- p+1
      df.CO[p,"g1"]<-genes[i]
      df.CO[p,"g2"]<-genes[j]
      df.CO[p,"co"]<-COcount
    }
  }
  return(df.CO)
}

#------------------------------------------------

# use vectorisation
GetCO_2_vect <- function(x){
  N <- nrow(x)
  M <- ncol(x)
  genes<-rownames(x)
  
  # create an empty data.frame in advance
  nb_pairs<-N*(N-1)/2
  df.CO <- data.frame(g1=rep("",nb_pairs),
                      g2=rep("",nb_pairs),
                      co=0,stringsAsFactors=F)
  p <- 0
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      vg1<-x[i,]
      vg2<-x[j,]
      # use vectorization
      COcount<-sum(vg1*vg2)
      p <- p+1
      df.CO[p,"g1"]<-genes[i]
      df.CO[p,"g2"]<-genes[j]
      df.CO[p,"co"]<-COcount
    }
  }
  return(df.CO)
}

#------------------------------------------------

# use matrix manipulation
GetCO_3_mtx <- function(x){
  # get product of mtx x and t(x)
  mtx.tmp<-x %*% t(x)
  # set the lower triangle of the matrix 
  # and diagonal to NA
  mtx.tmp[lower.tri(mtx.tmp, diag = TRUE)] <- NA 
  # melt the matrix (from reshape2 package)
  df.CO <- melt(mtx.tmp, na.rm = TRUE)
  colnames(df.CO) <- c("g1", "g2", "co")
  
  # next lines just to have same order, data types
  # and rownames as other GetCO functions
  # order on columns g1, g2
  df.CO<-df.CO[order(df.CO$g1,df.CO$g2),]
  # remove factors
  df.CO$g1<-as.character(df.CO$g1)
  df.CO$g2<-as.character(df.CO$g2)
  # remove rownames
  rownames(df.CO)<-NULL
  return(df.CO)
}

#------------------------------------------------

BenchMark_GetCO <- function(vbenchmark, times=1, 
                            maxsec=50, run=T, plot=T){
  library(microbenchmark)
  library(reshape2)

  funcs<-c("naive", "mtx", "allocmem", "vect", 
           "mtx_multipl")
  N<-length(vbenchmark)
  M<-length(funcs)
  nc.lims<- rep(max(vbenchmark*vbenchmark*4), M)
  
  
  if (run) {
    # create dataframe to store results of 
    # benchmark
    df.benchmark<-data.frame(ng=vbenchmark*4, 
                             ns=vbenchmark, 
                             nc=vbenchmark*
                               vbenchmark*4)
    for (f in funcs){
      df.benchmark[[f]]<-rep(NA,N)
    }
    
    # run functions for different matrix sizes
    for (i in 1:N){
      res<-CreateGeneSampleDF(ng=df.benchmark$ng[i], 
                              ns=df.benchmark$ns[i], 
                              prob1=0.4)
      df<-res$df
      mtx<-res$mtx
      nc<-df.benchmark$nc[i]
      
      bm<-microbenchmark(
        "naive"={df.count<-ifelse(nc>nc.lims[1],0,
                           GetCO_0(df))},
        "mtx"={df.count<-ifelse(nc>nc.lims[2],0,
                         GetCO_0(mtx))},
        "allocmem"={df.count<-ifelse(nc>nc.lims[3],0,
                              GetCO_1_alloc_mem(mtx))},
        "vect"={df.count<-ifelse(nc>nc.lims[4],0,
                          GetCO_2_vect(mtx))},
        "mtx_multipl"={df.count<-ifelse(nc>nc.lims[5],
                                 0,GetCO_3_mtx(mtx))},
        times=times, unit="s"
      )
      bm.summ<-summary(bm)
      
      df.benchmark[i,funcs] <- bm.summ$mean
      # if mean > maxsec sec, set future limit to nc+1
      nc.lims[which(bm.summ$mean>maxsec)]<-nc+1
      ix<-which(nc.lims<nc)
      df.benchmark[i,funcs[ix]]<-NA
    }
    
    save(df.benchmark, nc.lims, file="dfbenchmark.R")
  }
  
  if (plot){
    require(RColorBrewer)
    nfuncs<-length(funcs)
    if (!run) load(file="dfbenchmark.R")
    cols<-brewer.pal(nfuncs,"Set1")
    #cols<-c("black", "green", "blue", "orange", "red")
    plot(NULL, 
         xlab="size of matrix (nb genes x nb samples)", 
         ylab="time (sec)", 
         xlim=c(0,max(df.benchmark$nc)),
         ylim=c(0, max(df.benchmark[,-(1:3)], na.rm=T)))
    legend("topleft",
           legend = c("dataframe", "matrix",
                      "mem_alloc", "vectorize",
                      "matrix_multipl"),
           col=cols,lwd=2,
           lty=1,cex=1, bty="n")
    
    for (i in 1:nfuncs){
      ix<-which(!is.na(df.benchmark[[funcs[i]]]))
      lines(df.benchmark$nc[ix],
            df.benchmark[[funcs[i]]][ix], col=cols[i],
            lwd=2)
    }
  }
  return(list(df.benchmark=df.benchmark, 
              nc.lims=nc.lims))
}

#=================================================
# Which rows have sum > 4 ?
#=================================================

CreateBigDF <- function (n=100){
  set.seed(10)
  col1 <- runif (n, 0, 2)
  col2 <- rnorm (n, 0, 2)
  col3 <- rpois (n, 3)
  col4 <- rchisq (n, 2)
  return(data.frame (col1, col2, col3, col4))
}

#------------------------------------------------

Row4_naive <- function (df){
  N <- nrow(df)
  output<-character()
  for (i in 1:N){
    if (sum(df[i,])>4) {
      output[i]<-"greater"
    } else {
      output[i]<-"less_or_equal"  
    }
  }
  return(output)
}

#------------------------------------------------

Row4_mem <- function (df){
  N <- nrow(df)
  output <- character(N)
  for (i in 1:N){
    if (sum(df[i,])>4) {
      output[i]<-"greater"
    } else {
      output[i]<-"less_or_equal"  
    }
  }
  return(output)
}

#------------------------------------------------

Row4_mem_def <- function (df){
  N <- nrow(df)
  output <- rep("less_or_equal",N)
  for (i in 1:N){
    if (sum(df[i,])>4) {
      output[i]<-"greater"  
    }
  }
  return(output)
}

#------------------------------------------------

Row4_precheck <- function (df){
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
  return(output)
}

#------------------------------------------------

Row4_precheck_mem <- function (df){
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
  return(output)
}

#------------------------------------------------

Row4_precheck_mem_def <- function (df){
  N <- nrow(df)
  output <- rep("less_or_equal",N)
  precheck<-rowSums(df)>4
  for (i in 1:N){
    if (precheck[i]) {
      output[i]<-"greater"  
    }
  }
  return(output)
}

#------------------------------------------------

Row4_sapply <- function (df){
  N <- nrow(df)
  output<-sapply(1:N,FUN = function(x){
    if (sum(df[x,])>4) return("greater")
    else return("less_or_equal")
  })
  return(output)
}

#------------------------------------------------

Row4_sapply_precheck <- function (df){
  N <- nrow(df)
  precheck<-rowSums(df)>4
  output<-sapply(1:N,FUN = function(x){
    if (precheck[x]) return("greater")
    else return("less_or_equal")
  })
  return(output)
}

#------------------------------------------------

Row4_apply <- function (df){
  output<-apply(df,1,FUN = function(x){
    if (sum(x)>4) return("greater")
    else return("less_or_equal")
  })
  return(output)
}

#------------------------------------------------

CheckFour <- function(v1,v2,v3,v4){
  if (sum(c(v1,v2,v3,v4))>4) {
    return("greater") 
  } else {
    return("less_or_equal")
  }
}

vCheckFour <-Vectorize(CheckFour, SIMPLIFY = TRUE)

Row4_vect <- function (df){
  output<-vCheckFour(df$col1,df$col2, 
                     df$col3,df$col4)
  return(output)
}

#------------------------------------------------

Row4_ifelse <- function (df){
  output<-ifelse(rowSums(df)>4,"greater",
                 "less_or_equal")
  return(output)
}

#------------------------------------------------

Row4_which <- function (df){
  N <- nrow(df)
  output <- rep("less_or_equal",N)
  output[which(rowSums(df)>4)]<-"greater"
  return(output)
}

#------------------------------------------------

Row4_rcpp <- function (df){
  output <- cppCheckFour(df)
  return(output)
}

#------------------------------------------------

manualcolors<-c('black', 'red2', 'orange', 
                'cornflowerblue', 'seagreen', 
                'darkolivegreen4', 'indianred1', 
                'tan4', 'darkblue', 'mediumorchid1',
                'firebrick4', 'yellowgreen', 
                'lightsalmon', 'blue',
                'mediumvioletred', 
                'darkolivegreen1', 
                'tomato3','#7CE3D8'
                )

#------------------------------------------------

BenchMark_Row4 <- function(vbenchmark, funcs, times=1, 
                           maxsec=50, run=T, plot=T,
                           txt=""){

  library(microbenchmark)
  library(Rcpp)
  sourceCpp("SpeedUpR.cpp")

  N<-length(vbenchmark)
  M<-length(funcs)
  nr.lims<-rep(max(vbenchmark),M)+1
  
  if (run) {
    # create dataframe to store results of 
    # benchmark
    df.benchmark<-data.frame(nr=vbenchmark)
    for (f in funcs){
      df.benchmark[[f]]<-rep(NA,N)
    }
    
    # run functions for different df sizes
    for (i in 1:N){
      R<-df.benchmark$nr[i]
      df<-CreateBigDF(R)
      ls<-lapply(funcs, FUN=function(f){
                        if (R > nr.lims[which(funcs==f)]){
                          fname<-"nrow" 
                        } else {
                          fname<-paste0("Row4_",f) 
                        }
                        bquote(eval(call(.(fname),df)))
                        })
      names(ls)<-funcs
      bm<-microbenchmark(list=ls,times=times, unit="s")
      bm.summ<-summary(bm)
      
      df.benchmark[i,funcs] <- bm.summ$median
      # if mean > maxsec sec, set future limit to nc+1
      nr.lims[which(bm.summ$mean>maxsec)]<-R+1
      ix<-which(nr.lims<R)
      df.benchmark[i,funcs[ix]]<-NA
    }
    
    save(df.benchmark, nr.lims, 
         file=paste0("dfbenchmark_Row4",txt,".R"))
  }
  
  if (plot){
    load(file=paste0("dfbenchmark_Row4",txt,".R"))
    pdf(file=paste0("dfbenchmark_Row4",txt,".pdf"),
        width=9, height=7)
    
    # order funcs on speed
    # take nr.lims, order on nr.lims, if tie,
    # order on times at nr.lims-1
    lim.mtch<-match(nr.lims-1,df.benchmark$nr)
    timings<-rep(0,M)
    for (i in 1:M){
      timings[i]<-df.benchmark[lim.mtch[i],i+1]
    }
    ord_func<-order(nr.lims,-1*timings)
    funcs<-funcs[ord_func]

    cols<-manualcolors[1:M][ord_func]
    plot(NULL, main=txt,
         log="x",
         xlab="number of rows", 
         ylab="time (sec)", 
         xlim=c(min(df.benchmark$nr), max(df.benchmark$nr)),
         ylim=c(0, max(df.benchmark[,-1], na.rm=T)))
    legend("topleft",
           legend = funcs,
           col=cols,lwd=2,pch=ord_func,
           lty=1,cex=1, bty="n")
    
    for (i in 1:M){
      ix<-which(!is.na(df.benchmark[[funcs[i]]]))
      lines(df.benchmark$nr[ix],
            df.benchmark[[funcs[i]]][ix], col=cols[i],
            lwd=2, pch=ord_func[i], type="b")
    }
  }
  dev.off()
  return(list(df.benchmark=df.benchmark, 
              nr.lims=nr.lims))
}

#=====================================================
# Function: FastAggr_1col(df, col, fcols, funcs, ncol)
# Calculates fast aggregate on one column and returns 
# result dataframe 
# - df: dataframe
# - col: column name to aggregate on
# - fcol: column names to apply function on
# - func: functions to apply (same length as fcol)
# - ncol: vector of names of new columns
#=====================================================
FastAggr_1col<-function(df,col,fcol,func,ncol) {
  result <- lapply(split(seq(nrow(df)), df[[col]]), 
            function(.d)
            {
              c(
                df[[col]][.d][1], 
                sapply(1:length(fcol), FUN=function(x){
                  do.call(func[x],
                          list(df[[fcol[x]]][.d]))
                })
              )
            })
  mat <- do.call('rbind',result)   
  df.res<-as.data.frame(mat)
  colnames(df.res)<-c(col,ncol)
  if (!is.numeric(mat)){
    for (c in colnames(df.res)){
      if (is.numeric(df.res[[c]])){
        df.res[[c]]<-as.numeric(df.res[[c]])
      }
    }
  } 
  return(df.res)
}