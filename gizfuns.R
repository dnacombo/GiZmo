gizlm <- function (formula, basename = 'noname', nblocks = 1000, family=gaussian(), asfactors = NULL, ...){
  df = readdf(basename,asfactors)
  nobs = nrow(df)
  
  formula <- update(formula,'y ~ .')
  

  f <- function (y,formula,df,specs,at){
    df$y <- y
    res <- lm(formula, df)
    return(res)
  }
  
  delifexist(paste0(basename,'_coefs.dat'))
  delifexist(paste0(basename,'_resids.dat'))
  
  pb <- txtProgressBar(min=0
                       ,max=ceiling(file.info(paste0(basename,'.dat'))$size/4/nobs/nblocks)
                       ,style = 3)
  pbi <- 0
  fid <- file(description = paste0(basename,'.dat'),open="rb" )
  while (T) {
    setTxtProgressBar(pb,pbi)
    dat <- readBin(con = fid,
                   what="numeric",
                   n=nblocks * nobs,
                   size=4,
                   endian="little")
    if (length(dat) == 0) break
    
    dim(dat) <- c(nobs,length(dat)/nobs)
    
    res <- apply(dat,2,f,formula,df,specs,at)
    
    coefs = sapply(res,function(x)coefficients(x))
    resids = sapply(res,function(x)residuals(x))

    dat2bin(paste0(basename,'_coefs.dat'),coefs)
    dat2bin(paste0(basename,'_resids.dat'),resids)

    pbi <- pbi + 1
  }
  close(fid)
  write.table(x = rownames(coefs),file = paste0(basename,'_coefnames.txt'),row.names = F,col.names = F)
}

gizlm_withcontrasts <- function (formula, basename = 'noname', nblocks = 1000, family=gaussian(), asfactors = NULL, ...){
  library(lsmeans)
  df = readdf(basename,asfactors)
  nobs = nrow(df)
  
  formula <- update(formula,'y ~ .')
  
  specs <- attr(terms(formula),'term.labels') # extract RHS of the formula
  specs <- specs[! grepl(':',specs)]       # remove any interaction
  # create conditions for continuous predictors (for factors, it will be all levels)
  at <- list()
  for (s in specs){
    if (!is.factor(df[s])) { 
      at[[s]] <- c(0,1)
    }
  }
  
  
  
  f <- function (y,formula,df,specs,at){
    df$y <- y
    res <- lm(formula, df)
    
    lsm <- summary(lsmeans(res,specs,at=at,data=df))
    return(list(res,lsm))
  }
  
  delifexist(paste0(basename,'_coefs.dat'))
  delifexist(paste0(basename,'_resids.dat'))
  
  pb <- txtProgressBar(min=0
                       ,max=ceiling(file.info(paste0(basename,'.dat'))$size/4/nobs/nblocks)
                       ,style = 3)
  pbi <- 0
  fid <- file(description = paste0(basename,'.dat'),open="rb" )
  while (T) {
    setTxtProgressBar(pb,pbi)
    dat <- readBin(con = fid,
                   what="numeric",
                   n=nblocks * nobs,
                   size=4,
                   endian="little")
    if (length(dat) == 0) break
    
    dim(dat) <- c(nobs,length(dat)/nobs)
    
    res <- apply(dat,2,f,formula,df,specs,at)
    
    coefs = sapply(res,function(x)coefficients(x[[1]]))
    resids = sapply(res,function(x)residuals(x[[1]]))
    conditions = sapply(res,function(x)x[[2]]$lsmean)
    conditions.SE = sapply(res,function(x)x[[2]]$SE)
    
    dat2bin(paste0(basename,'_coefs.dat'),coefs)
    dat2bin(paste0(basename,'_resids.dat'),resids)
    dat2bin(paste0(basename,'_conditions.dat'),conditions)
    dat2bin(paste0(basename,'_conditions_SE.dat'),conditions.SE)
    
    pbi <- pbi + 1
  }
  close(fid)
  write.table(x = rownames(coefs),file = paste0(basename,'_coefnames.txt'),row.names = F,col.names = F)
}

gizlme <- function (formula, random, basename = 'noname', nblocks = 10, asfactors = NULL, ...){
  library(nlme)
  df = readdf(basename,asfactors)
  
  f <- function (y,formula,df,random){
    df$y = y
    res = lme(formula, df, random)
  }
  nobs = nrow(df)
  
  delifexist(paste0(basename,'_fixefs.dat'))
  delifexist(paste0(basename,'_ranefs.dat'))
  delifexist(paste0(basename,'_resids.dat'))
  delifexist(paste0(basename,'_contrastsF.dat'))
  delifexist(paste0(basename,'_contrastsp.dat'))
  
  pb <- txtProgressBar(min=0
                       ,max=ceiling(file.info(paste0(basename,'.dat'))$size/4/nobs/nblocks)
                       ,style = 3)
  pbi <- 0
  fid <- file(description = paste0(basename,'.dat'),open="rb" )
  res <- list()
  fixefs <- numeric()
  ranefs <- numeric()
  resids <- numeric()
  anovas <- list()
  ctrstF <- numeric()
  ctrstp <- numeric()
  while (T) {
    setTxtProgressBar(pb,pbi)
    dat <- readBin(con = fid,
                   what="numeric",
                   n=nblocks * nobs,
                   size=4,
                   endian="little")
    if (length(dat) == 0) break
    
    dim(dat) <- c(nobs,length(dat)/nobs)
    allres <- tryCatch({res <- apply(dat,2,f,formula,df,random)
                        #                         if (pbi==2) stop('test')
                        fixefs <- sapply(res,fixef)
                        ranefs <- as.numeric(unlist(sapply(res,ranef)))
                        resids <- sapply(res,resid)
                        anovas <- sapply(res,anova)
                        #     ctrstnumdf = sapply(res,function(x){anova(x)[,1]})
                        #     ctrstden = sapply(res,function(x){anova(x)[,2]})
                        ctrstF <- apply(anovas,2,function(x){x[[3]]})
                        ctrstp <- apply(anovas,2,function(x){x[[4]]})
    },
    error=function(e){
      fixefs <<- rep(NaN,length(fixefs))
      ranefs <<- rep(NaN,length(ranefs))
      resids <<- rep(NaN,length(resids))
      ctrstF <<- rep(NaN,length(ctrstF))
      ctrstp <<- rep(NaN,length(ctrstp))
    })
    dat2bin(paste0(basename,'_fixefs.dat'),fixefs)
    dat2bin(paste0(basename,'_ranefs.dat'),ranefs)
    dat2bin(paste0(basename,'_contrastsF.dat'),ctrstF)
    dat2bin(paste0(basename,'_contrastsp.dat'),ctrstp)
    dat2bin(paste0(basename,'_resids.dat'),resids)
    #         if (pbi > 3){
    #         break}
    pbi <- pbi + 1
  }
  close(fid)
  write.table(x = rownames(fixefs),file = paste0(basename,'_fixefsnames.txt'),row.names = T,col.names = F,sep = '\t')
  write.table(x = rownames(ranef(res[[1]])),file = paste0(basename,'_ranefsnames.txt'),row.names = F,col.names = F,sep = '\t')
  write.table(x = rownames(anova(res[[1]])),file = paste0(basename,'_contrastsnames.txt'),row.names = F,col.names = F,sep = '\t')
}

gizfastlm <- function (formula, basename = 'noname', nblocks = Inf, asfactors = NULL, ...){
  library(MASS)
  
  df = readdf(basename,asfactors)
  nobs = nrow(df)
  
  x <- model.matrix(formula, data = df)
 
  delifexist(paste0(basename,'_design.dat'))
  dat2bin(paste0(basename,'_design.dat'),x)
  
  gx <- ginv(x)
  
  
  f <- function (y){
    coefs <- gx %*% y
    resids <- y - (x %*% coefs)
    list(coefs=coefs,resids=resids)
  }
  delifexist(paste0(basename,'_coefs.dat'))
  delifexist(paste0(basename,'_resids.dat'))
  
  pb <- txtProgressBar(min=0
                       ,max=ceiling(file.info(paste0(basename,'.dat'))$size/4/nobs/nblocks)
                       ,style = 3)
  pbi <- 0
  fid <- file(description = paste0(basename,'.dat'),open="rb" )
  while (T) {
    setTxtProgressBar(pb,pbi)
    dat <- readBin(con = fid,
                   what="numeric",
                   n=nblocks * nobs,
                   size=4,
                   endian="little")
    if (length(dat) == 0) break
    
    dim(dat) <- c(nobs,length(dat)/nobs)
    res <- apply(dat,2,f)
    
    coefs = sapply(res,function(x){x$coefs})
    resids = sapply(res,function(x){x$resids})
    
    dat2bin(paste0(basename,'_coefs.dat'),coefs)
    dat2bin(paste0(basename,'_resids.dat'),resids)
    pbi <- pbi + 1
  }
  close(fid)
  write.table(x = attr(x,'dimnames')[[2]],file = paste0(basename,'_coefnames.txt'),row.names = F,col.names = F)
}

gizglm <- function (formula, basename = 'noname', nblocks = 1000, family=gaussian(), asfactors = NULL, ...){
  
  df = readdf(basename,asfactors)
  nobs = nrow(df)
  
  x <- as.numeric(as.matrix(model.frame(formula, data = df)))
  dim(x) <- c(nobs,length(x)/nobs)
  delifexist(paste0(basename,'_design.dat'))
  dat2bin(paste0(basename,'_design.dat'),x)
  
  f <- function (y){
    res = glm.fit(x,y,family=family)
  }
  delifexist(paste0(basename,'_coefs.dat'))
  delifexist(paste0(basename,'_resids.dat'))
  
  pb <- txtProgressBar(min=0
                       ,max=ceiling(file.info(paste0(basename,'.dat'))$size/4/nobs/nblocks)
                       ,style = 3)
  pbi <- 0
  fid <- file(description = paste0(basename,'.dat'),open="rb" )
  while (T) {
    setTxtProgressBar(pb,pbi)
    dat <- readBin(con = fid,
                   what="numeric",
                   n=nblocks * nobs,
                   size=4,
                   endian="little")
    if (length(dat) == 0) break
    
    dim(dat) <- c(nobs,length(dat)/nobs)
    res <- apply(dat,2,f)
    
    coefs = sapply(res,coef)
    resids = sapply(res,resid)
    
    dat2bin(paste0(basename,'_coefs.dat'),coefs)
    dat2bin(paste0(basename,'_resids.dat'),resids)
    pbi <- pbi + 1
  }
  close(fid)
  write.table(x = rownames(coefs),file = paste0(basename,'_coefnames.txt'),row.names = F,col.names = F)
}

gizlmBF <- function (formula, basename = 'noname', nblocks = 10, asfactors = NULL, whichRandom = NULL, ...){
  
  library(BayesFactor)
  options(BFprogress = F)
  df = read.table(paste0(basename,'.df'),header=T)
  if (!any(is.na(asfactors))) {
    for (i in 1:length(asfactors)){
      if (asfactors[i]){
        df[,i] = as.factor(df[,i])
      } else {
        df[,i] = as.numeric(df[,i])
      }
    }
  }
  nobs = nrow(df)
  f <- function (y,formula,df){
    df$y = y
    anovaBF(formula, df,whichRandom,multicore = T, noSample = F)
  }
  delifexist(paste0(basename,'_BF.dat'))
  
  fid <- file(description = paste0(basename,'.dat'),open="rb" )
  while (T) {
    dat <- readBin(con = fid,
                   what="numeric",
                   n=nblocks * nobs,
                   size=4,
                   endian="little")
    if (length(dat) == 0) break
    
    dim(dat) <- c(nobs,length(dat)/nobs)
    res <- apply(dat,2,f,formula,df)
    BF <- sapply(res,function(x){extractBF(x,onlybf=T)})
    dat2bin(paste0(basename,'_BF.dat'),BF)
  }
  close(fid)
  tmp = extractBF(res[[1]])
  write.table(x = rownames(tmp),file = paste0(basename,'_BFnames.txt'),row.names = F,col.names = F)
}


## Helpers

dat2bin <- function(fname,dat){
  fud <- file(description = fname, open="ab")
  ok <- writeBin(object = as.vector(dat), con=fud, size=4, endian="little")
  close(fud)
}

delifexist <- function(fname){
  if (file.exists(fname)) file.remove(fname)
}

readdf <- function(basename,asfactors){
  
  df = read.table(paste0(basename,'.df'),header=T)
  if (!any(is.na(asfactors))) {
    if (is.character(asfactors)){
      # named columns
      for (i in 1:length(asfactors)){
        df[,asfactors[i]] = as.factor(df[,asfactors[i]])
      }
    }
    else{
      # logicals
      for (i in 1:length(asfactors)){
        if (asfactors[i]){
          df[,i] = as.factor(df[,i])
        } else {
          df[,i] = as.numeric(df[,i])
        }
      }
    }
  }
  df
}


