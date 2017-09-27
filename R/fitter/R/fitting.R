INTERCEPT="Intercept"

read.sampled <- function(d, sampler.dir, ignore.empty.file){
  sampled.file = sub("YYYYMMDD", d, file.path( sampler.dir, "RAW.YYYYMMDD.rds"))
  if( !file.exists(sampled.file)){
    if( ignore.empty.file){
      return(NULL)
    }else{
      stop("Sampled data does not exists")
    }
  }
  sampled = readRDS(sampled.file)
}

fit.coef.OLS <- function(data, xy.spec, d = "19991231", grp.name = "00:00:00.000", pa.list = list()){
  x.vars = xy.spec$x.vars
  fit.y.vars = xy.spec$fit.y.vars
  wgt.var = xy.sepc$wgt.var
  
  if( !all(c(x.vars, fit.y.vars) %in% dimnames(data)[[2]])){
    stop(paste("fitting vars not complete", paste(setdiff(c(x.vars, fit.y.vars), dimnames(data)[[2]]),collapse = ",")))
  }
  
  intercept= !!as.integer(null.replace(pa.list$INTERCEPT,0))
  
  data[is.na(data) | is.infinite(data)] = 0
  
  ## fit as simple regression
  y = data[,fit.y.vars,drop=FALSE]
  X = data[,x.vars,drop=FALSE]
  x = cbind(rep(1,dim(x)[1]), x)
  colnames(x)[1] = INTERCEPT
  if ( is.null(wgt.var)){
    w = rep(1, dim(x)[1])
  }else{
    w = data[,wgt.var, drop=FALSE]
  }
  xx <- crossprod(x,x*w)
  xy <- crossprod(x,y*w)
  nonzero.x = which(abs(apply(xx, 2, sum, na.rm=TRUE)) > 1e-09)
  
  nz.xx = xx[nonzero.x,nonzero.x,drop=FALSE]
  nz.xy = xy[nonzero.x,,drop=FALSE]
  
  ridge = null.replace(pa.list$RIDGE, 1e-09)
  nz.coefs = crossprod(solve(nz.xx + diag(dim(nz.xx)[1])*ridge),nz.xy)
  
  coefs = xy*0
  coefs[nonzero.x,] = nz.coefs
  dim(coefs) = c(dim(coefs),1,1)
  dimnames(coefs) = list(c(INTERCEPT,x.vars),fit.y.vars,d,grp.name)
  coefs = aperm(coefs,c(3,4,1,2))
  names(dimnames(coefs)) = c("D","G", "X", "Y")
  coefs
}


fit.coef.LM <- function(data, xy.spec, d = "19991231", grp.name="00:00:00.000", pa.list=list()){
  x.vars = xy.spec$x.vars
  fit.y.vars = xy.spec$fit.y.vars
  wgt.var = xy.spec$wgt.var
  
  if( !all(c(x.vars, fit.y.vars) %in% dimnames(data)[[2]])){
    stop(paste("fitting vars not complete", paste(setdiff(c(x.vars, fit.y.vars), dimnames(data)[[2]]),collapse = ",")))
  }
  
  intercept = !!as.integer(null.replace(pa.list$INTERCEPT,0))
  
  data[is.na(data) | is.infinite(data)] = 0
  
  ##fit as simple regression
  y = data[,fit.y.vars,drop=FALSE]
  x = data[,x.vars,drop=FALSE]
  x = cbind(rep(1,dim(x)[1]),x)
  colnames(x)[1] = INTERCEPT
  if( is.null(wgt.var)){
    w = rep(1,dim(x)[1])
  }else{
    w = data[,wgt.var,drop=FALSE]
  }
  
  coefs = t(t(coef(lm(y~0+x,weights=w))))
  coefs[is.na(coefs)]=0
  dim(coefs) = c(dim(coefs),1,1)
  dimnames(coefs) = list(c(INTERCEPT, x.vars), fit.y.vars, d, grp.name)
  coefs = aperm(coefs, c(3,4,1,2))
  names(dimnames(coefs)) = c("D","G","X","Y")
  coefs
}

fit.coef.NNLS <- function(data,xy.spec, d = "19991231", grp.name = "00:00:00.000", pa.list=list()){
  x.vars = xy.spec$x.vars
  fit.y.vars = xy.spec$fit.y.vars
  wgt.var = xy.spec$wgt.var
  
  if( !all(c(x.vars, fit.y.vars) %in% dimnames(data)[[2]])){
    stop(paste("fitting vars not complete", paste(setdiff(c(x.vars,fit.y.vars),dimnames(data)[[2]]),collapse = ",")))
  }
  
  intercept= !!as.integer(null.replace(pa.list$INTERCEPT,0))
  
  data[is.na(data) | is.infinite(data)] = 0
  
  ## fit as simple regression
  y = data[, fit.y.vars, drop=FALSE]
  x = data[, x.vars, drop=FALSE]
  x = cbind(rep(1,dim(x)[1]),x)
  colnames(x)[1] = INTERCEPT
  if(is.null(wgt.var)){
    w = rep(1, dim(x)[1])
  }else{
    w = data[,wgt.var,drop=FALSE]
  }
  
  coefs = t(tcoef(nnls(x*w,y*w)))
  dim(coefs) = c(dim(coefs),1,1)
  dimnames(coefs) = list(c(INTERCEPT, x.vars), fit.y.vars, d, grp.name)
  coefs = aperm(coefs, c(3,4,1,2))
  names(dimnames(coefs)) = c("D","G","X","Y")
  coefs
}

fit.coef.cvnet <- function(data, xy.spec, d="19991231", grp.name="00:00:00.000", pa.list=list(MODE="FP",IDEAL_FR=1.02),cores=1){
  library(glmnet)
  library(doMC)
  registerDoMC(cores = cores)
  
  ideal.fr = as.numeric(null.replace(pa.list$IDEAL_FR,1.02))
  x.vars = xy.spec$x.vars
  fit.y.vars = xy.spec$fit.y.vars
  wgt.var = xy.spec$wgt.var
  
  if( !all(c(x.vars, fit.y.vars) %in% dimnames(data)[[2]])){
    stop(paste("fitting vars not complete",paste(setdiff(c(x.vars, fit.y.vars),dimnames(data)[[2]]),collapse = ",")))
  }
  
  intercept=!!as.integer(null.replace(pa.list$INTERCEPT,0))
  if(length(x.vars)==1){
    coefs = fit.coef.OLS(data,d=d,xy.spec=xy.spec,grp.name=grp.name,pa.list=list(INTERCEPT=intercept))
  }else{
    coefs = laply(fit.y.vars,function(y.var){
      x = data[,x.vars,drop=FALSE]
      y = data[,y.var, drop=FALSE]
      if( is.null(wgt.var)){
        w = t(t(rep(1,dim(x)[1])))
      }else{
        w = data[,wgt.var,drop=FALSE]
      }
      
      x[is.na(x)|is.infinite(x)]=0;
      y[is.na(y)|is.infinite(y)]=0;
      w[is.na(w)|is.infinite(w)]=0;
      
      cvfit = cv.glmnet(x=x,y=y,weights=w,parallel=cores>1, intercept = intercept)
      coefm = as.matrix(coef(cvfit,s=cvfit$lambda))
      dimnames(coefm)[[1]][1] = INTERCEPT
      a = data[,x.vars,drop=FALSE] %*% coefm[-1,]
      fr = cov(y,a,use="pair")/apply(a,2,var)
      fr.idx = which.min(abs(fr-ideal.fr))
      lmin.idx = which.min(abs(cvfit$lambda - cvfit$lambda.min))
      coefm[,min(fr.idx,lmin.idx),drop=FALSE]
    },.drop=FALSE)
    dimnames(coefs)[[1]] = fit.y.vars
    dim2 = dimnames(coefs)[[2]]
    dim(coefs) = c(dim(coefs),1)
    dimnames(coefs) = list(fit.y.vars, dim2, d,grp.name)
    coefs = aperm(coefs, c(3,4,2,1))
    names(dimnames(coefs))=c("D","G","X","Y")
  }
  coefs
}

fit.coef.POLY <- function(data,xy.spec, d="19991231", grp.name="00:00:00.000",pa.list=list()){
  
}

fitBeta <- function(stDate,edDate,sampler.dir,model.cfg=list(ALL="ALL_"),
                    ins.days =getTradingDayRange(20000101,20141231)[getTradingDayRange(20000101,20141231)%/%100%%2==1],
                    mode="ALL",
                    fit.y.vars = c("fwd.Ret.DAILY.1","fwd.Ret.DAILY.5"),wgt.var=NULL,
                    ignore.empty.file=TRUE,
                    fit.func = fit.coef.cvnet, fit.para = list(), use.cache=TRUE, cache.coef.ptn=NULL,
                    coef.max=Inf, coef.min=-Inf,coes=1,verbose=TRUE){
  
  if( is.character(fit.func)){
    fit.func = switch(fit.func, CVNET=fit.coef.cvnet, OLS = fit.coef.OLS, LM = fit.coef.LM, NNLS = fit.coef.NNLS, stop(paste("Unsupported fit method name",fut.func)) )
  }else if( !is.function(fit.func)){ stop("Invalid fit func type!")}
  
  library(doMC)
  doMC::registerDoMC(cores = cores)
  
  tradingDays = getTradingDayRange(stDate,edDate)
  ins.tradingDays = tradingDays[ tradingDays %in% ins.days]
  
  if(mode="ALL"){
    sampled=NULL
    cache.file = file.path(sampler.dir, paste("cache.RAW.",stDate,".to.",edDate,".N",length(ins.tradingDays),".rds",sep=""))
    if(use.cache){
      if(file.exists(paste(cache.file,".end",sep="")) && file.exists(cache.file)){
        sampled=readRDS(cache.file)
      }
    }
    
    if(is.null(sampled)){
      data.g = llply(ins.tradingDays, function(d){
        sampled = read.sampled(d, sampler.dir=sampler.dir, ignore.empty.file = ignore.empty.file)
      },.parallel=cores>1)
      which.null = unlist(lapply(data.g, is.null))
      data.g = data.g[!which.null]
      
      sampled = llply(names(data.g[[1]]), function(g){
        do.call(rbind, lapply(data.g, function(dd) dd[[g]]) )
      })
      names(sampled) = names(data.g[[1]])
      
      ## save cache
      if(use.cache){
        dir.create(dirname(cache.file),FALSE,TRUE)
        saveRDS(sampled, cache.file)
        file.create(paste(cache.file,".end",sep=""))
      }
    }
    coef.d = panel.combine(llply(names(model.cfg),function(model){
      coefs.g = abind(llply(names(sampled),function(g){
        fitdata = sampled[[g]]
        fitdata = fitdata[sapply(fitdata, is.numeric)]
        
        
        fitdata = as.matrix(fitdata)
        x.vars = model.cfg[[model]]
        print(x.vars)
        if(is.null(x.vars)){
          x.vars = model
        }else if(identical(x.vars,"ALL_")){
          x.vars = setdiff(dimnames(fitdata)[[2]],c(fit.y.vars,wgt.var) )
        }else if(!is.null(x.vars) && length(x.vars)==1 && length(grep("\\*",x.vars))>0){
          x.vars = setdiff(grep(x.vars,dimnames(fitdata)[[2]],value=TRUE),c(fit.y.vars,wgt.var))
        }
        
        ## handle high order condition
        order = null.replace(fit.para$order,1)
        if(length(x.vars)>1 & order>1) stop("high order fitting is only compatible with 1 independent variable")
        
        if(order > 1){
          new.fitdata.list <- lapply(2:order,function(o){
            fitdata[,x.vars,drop=FALSE]^o
          })
          new.fitdata <- do.call(cbind,new.fitdata.list)
          colnames(new.fitdata) <- paste0(x.vars,".",2:order)
          fitdata <- cbind(fitdata,new.fitdata)
          x.vars <- c(x.vars,paste0(x.vars,".",2:order))
        }
        
        xy.spec=list(x.vars = x.vars, fit.y.vars = fit.y.vars)
        coefs = fit.func(fitdata,xy.spec,d="19991231",grp.name=g,pa.list=fit.para)
      }),along=2)
      names(dimnames(coefs.g)) = c("D","G","X","Y")
      coefs.g = panel.add.dim(coefs.g,"M",model)
    }),default=0)
  }else if(mode = "ROLL"){
    cat("Fitting")
    cache.set = !is.null(cache.coef.ptn)
    use.cache = use.cache && cache.set
    
    coef.d = abind(llply(ins.tradingDays, function(d){
      cat(d)
      cat(" ")
      coef.out.file = sub("YYYYMMDD",d,cache.coef.ptn)
      
      if( use.cache && file.exists(coef.out.file)){
        return(panel.read(coef.out.file, verbose = verbose))
      }
      sampled = read.sampled(d,sampler.dir = sampler.dir, ignore.empty.file=ignore.empty.file)
      coefs.m = panel.combine(llply(names(model.cfg),function(model){
        coefs.g = abind(llply(names(sampled),function(g){
          fitdata = sampled[[g]]
          fitdata =fitdata[sapply(fitdata,is.numeric)]
          fitdata = as.matrix(fitdata)
          x.vars = model.cfg[[model]]
          if(is.null(x.vars) || identical(x.vars,"ALL_")){
            x.vars = setdiff( dimnames(fitdata)[[2]],c(fit.y.vars,wgt.var))
          }else if(!is.null(x.vars) && length(x.vars)==1 && length(grep("\\*",x.vars))>0) x.vars = setdiff(grep(x.vars,dimnames(fitdata)[[2]],value=TRUE),c(fit.y.vars,wgt.var))
          xy.spec=list(x.vars=x.vars, fit.y.vars = fit.y.vars)
          coefs = fit.func(fitdata,xy.spec,d=d,grp.name=g,pa.list=fit.para)
        }),along=2)
        names(dimnames(coefs.g))=c("D","G","X","Y")
        coefs.g = panel.add.dim(coefs.g,"M",model)
      }),default = 0)
      
      if(cache.set){
        dir.create(dirname(coef.out.file),FALSE,TRUE)
        panel.write(coefs.m, coef.out.file, verbose = verbose)
      }
      coefs.m
    }, .parallel = cores>1),along=1)
    cat("\n")
    
    if(any(is.na(coef.d))) stop("Coef error")
    
    if( !is.null(fit.para$WINDOW)){
      W = as.integer(fit.para$WINDOW)
      HL = as.integer(fit.para$HL)
      if(!is.null(HL) & length(HL)!=0){
        W = (0.5)^((1:W)/HL)
        w = w/sum(w)
      }else{
        w = rep(1/W,W)
      }
      print(paste0("smoothing ",W," HL=",HL))
      pad = rep(0,W-1)
      
      ## apply window average and set leading NA to zero
      coef.d.new = aperm(aaply(coef.d, 2:5, function(x){fiter(c(pad,x),w,sides=1)[-(1:(W-1))]},.drop=FALSE),c(5,1:4))
      dimnames(coef.d.new)[[1]]=dimnames(coef.d)[[1]]
      names(dimnames(coef.d.new)) = names(dimnames(coef.d))
      coef.d = coef.d.new
    }
    if(!is.null(fit.para$NDAY)){
      step = max(1,as.integer(fit.para$NDAY))
      print(paste0("stepping ",step))
      coef.d = coef.d[seq(1,dim(coef.d)[1],by=step),,,,,drop=FALSE]
    }
  }else{
    stop(paste("mode not supported",mode))
  }
  names(dimnames(coef.d)) = c("D","G","X","Y","M")
  coef.d[coef.d>coef.max] = coef.max
  coef.d[coef.d<coef.min] = coef.min
  return(coef.d)
}
  
## TODO: handle the group name to group mapping, for now assume same to all K,T
gen.alp.on.coef <- function(stDate,edDate,coef,term.path,
                            model,grp.name,alpha.path,alphaname=model,
                            lag=1,mkt="CHINA_STOCK",
                            include.days=NULL,
                            use.cache=TRUE,
                            verbose=FALSE,
                            alp.only=FALSE,
                            cores=1,
                            indExpoPath = "/mnt/analysis/ics/risk/CHINA_STOCK/DAILY/industry/YYYYMMDD.h5",
                            onecache.only=FALSE,fit.para = list()){
  
  library(doMC)
  doMC::registerDoMC(cores = cores)
  
  if(length(grp.name)!=1 || length(model)!=1)stop("Only support one group or model")
  if(alp.only && stDate != edDate){
    use.cache =FALSE ## do not use cache when alp.only ==TRUE
    stop("Only support one day call when alp.only=TRUE")
  }
  
  coef = coef[,,,,model,drop=FALSE]
  
  tradingDays = getTradingDayRange(stDate,edDate)
  if( !is.null(include.days)) tradingDays = tradingDays[tradingDays %in% include.days]
  
  alpha.path = sub("MODELNAME", model, alpha.path)
  alpha.path = sub("GROUPNAME", grp.name, alpha.path)
  
  alpha.names = sub("(fwd)\\.(Ret)\\.(.*)\\.(.*)", paste(alpname,".\\3.\\4",sep=""),dimnames(coef)[[4]] )
  
  readOneCache = FALSE
  if(onecache.only && file.exists(sub("YYYYMMDD", paste0("cache_from_",tradingDays[1],'_to_',tail(tradingDays,1)),term.path))){
    rawAlphaCache = panel.read(sub("YYYYMMDD",paste0('cache_from_',tradingDays[1],'_to_',tail(tradingDays,1)),term.path),verbose=verbose)
    readOneCache = TRUE
  }
  order = null.replace(fit.para$order,1)
  alpha.files = llply(tradingDays, function(d){
    alpha.file = sub("YYYYMMDD", d, alpha.path)
    if(!alp.only && use.cache && file.exists(alpha.file)) return(NULL);
    
    term.file = sub("YYYYMMDD", d, term.path)
    if(any(!file.exists(term.file))){
      print(paste("WARN: term file does not exists", term.file[!file.exists(term.file)]))
      if(alp.only) stop("Cannot load term to generate alp")
      return(NULL)
    }
    
    if(!readOneCache){
      terms = readicsdata(d,vers=NULL,paths = term.file,verbose = verbose)
    }else{
      terms = rawAlphaCache[,as.character(d),,,drop=FALSE]
    }
    
    terms[is.na(terms)]=0
    coefd = as.integer(dimnames(coef)[[1]])
    coefd = coefd[max(1,which(coefd<=d ) - lag)]
    if(gsub('\\d+',"",grp.name) == "G"){
      indExpo = panel.read(sub("YYYYMMDD",d,indExpoPath)) # need ind expo to determine which stock is in which group
      grp = icsUtil:getIndGroup()[[grp.name]]
      sg = llply(grp,function(g){
        k = which(indExpo==1, arr.ind =T)[,4] %in% g
        names(k) = dimnames(indExpo)$K
        k
      })
    }
    if(order > 1 && length(alpha.names)==1){
      new.terms.list = llply(2:order,function(o){
        tmp = terms^o
        dimnames(tmp)[[4]] = paste0(dimnames(terms)[[4]],'.',o)
        tmp
      })
      terms = panel.combine(append(list(terms),new.terms.list))
    }else if(order>1 & length(alpha.names)!=1){
      stop('high order mode is noly compatible with 1 dependent variable')
    }
    
    alp.v = aperm(aaply(terms,2:3,function(KV){
      if(is.null(dim(KV))){
        KV =t(t(KV))
        dimnames(KV)[[2]] = dimnames(terms)[[4]]
        dimnames(KV)[[1]] = dimnames(terms)[[1]]
      }
      
      if(grp.name == "EACH"){
        common.keys = intersect(dimnames(KV)[[1]], dimnames(coef)[["G"]])
        alp = abind( llply(dimnames(coef)[["Y"]], function(fit.y){
          alp.y = aaply(KV[common.keys, dimnames(coef)[['X']][-1] ]*adrop(coef[as.character(coefd),common.keys,-1,fit.y,model,drop=FALSE],c(1,4,5)),1,sum,na.rm=TRUE,drop=FALSE)
        }),along = 2)
        dimnames(alp)[[2]] = dimnames(coef)[['Y']]
        alp
      }else if(gsub('\\d+',"",grp.name) = 'G'){
        common.keys = intersect(dimnames(KV)[[1]],names(sg[[1]]))
        alp = abind(llply(dimnames(coef)$Y, function(fit.y){
          alp.y = ldply(dimnames(coef)$G, function(cg){
            KVtimeBeta = sweep(KV[common.keys,,drop=FALSE][sg[[cg]][common.keys],dimnames(coef)$X[-1],drop=FALSE],MARGIN = 2,coef[,cg,-1,,model],`*`)
            adply(KVtimeBeta,1,sum,.id=NULL)
          })
        }),along=2)
        dimnames(alp)[[2]] = dimnames(coef)$Y
        alp
      }else{
        beta = adrop( coef[ as.character(coefd),grp.name,-1,,model,drop=FALSE],c(1,2,5))
        alp = KV[,rownames(beta),drop=FALSE] %*% beta
      }
    },.parallel = coes>1,.drop=FALSE),c(3,1,2,4))
    dimnames(alp.v)[[4]] = alpha.names
    names(dimnames(alp.v)) = c("K","D","T","V")
    if(alp.only | onecache.only) return(alp.v)
    dir.create(dirname(alpha.file),FALSE,TRUE)
    panel.write(alp.v, alpha.file, verbose = verbose, overwrite = !use.cache)
    alpha.file
  }, .parallel=cores>1)
  
  if(alp.only){
    alpha.files[[1]]
  }else{
    if(onecache.only){
      alpha.combine = panel.combine(alpha.files)
      dir.create(dirname(alpha.path),FALSE,TRUE)
      panel.write(alpha.combine,file.path(dirname(alpha.path),paste0("cache_from_",tradingDays[1],"_to_",tail(tradingDays,1),'.h5')),verbose=verbose,overwrite=TRUE)
    }
    alpha.names
  }
}
  










