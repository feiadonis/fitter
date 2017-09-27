library(abind)
library(plyr)
library(ggplot2)
library(stringr)

library(icsUtil)
library(sampler)
library(fitter)
library(alpEval)

if("icsUtil" %in% loadedNamespaces()) detach("package:icsUtil",unload=TRUE);library(icsUtil)
if("sampler" %in% loadedNamespaces()) detach("package:sampler",unload=TRUE);library(sampler)
if("fitter" %in% loadedNamespaces()) detach("package:fitter",unload=TRUE);library(fitter)
if("alpEval" %in% loadedNamespaces()) detach("package:alpEval",unload=TRUE);library(alpEval)

cfg = list(
  univ = "samllUniv",
  mkt = "CHINA_FUTURE",
  Freq = "DAILY",
  root.dir = "~/proj/myFP/alphazf",
  univ.ver = "windzf",
  periods =  5,
  #grpType = "ALL",
  #grpName = "nonFinance"
  dates = list(sdate = 20090101, edate = 20141231, os.sdate = "20090101", os.edate = "20170615",
               omit.days = NULL, ins.grp=c("INS1"),oos.grp=c("INS1","INS2","OOS")),
  
  sampler.list = list(
    ALL=list(group.type="ALL", dir.list = c("/mnt/analysis/ics/alpha/CHINA_STOCK/DAILY/zfb_smallUniv/YYYYMMDD.h5"))
  ),
  
  fitter.list = list(
    OLS = list(mode="ALL",fit.func="OLS"),
    LM = list(mode="ALL",fit.func="LM"),
    CVNET=list(mode="ALL",fit.func="CVNET"),
    NNLS = list(mode="ALL",fit.func="NNLS"),
    OLSROLL63 = list(mode="ROLL",fit.func="OLS",WINDOW=63),
    CVNETROLL = list(mode = "ROLL", fit.func = "CVNET", WINDOW = 63),
    CVNETROLL100HL30N1=list(mode="ROLL",fit.func="CVNET",WINDOW=100,HL=30,NDAY=1)
    ),
  
  model.list = list(
    ALL = "ALL_",
    zfb="alphazf2b"
  ),
  
  fit.list = list(
    alpha = list(sampler="ALL",fitter="CVNET",model="zfb")
  )
)

if(TRUE){
  fit.process(cfg=cfg,rpt.name="test.pdf",ver="v1")
}

