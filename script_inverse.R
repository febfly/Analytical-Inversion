#Compute Bayesian inversion to optimize emission spatial distributions
#Computation is performed to achieve monthly emissions.



#Load processed tropomi data
load.tropomi=function(datestr,type='fullinv3'){
   if (type!='default'){
      path=paste0('/net/seasasfs02/srv/export/seasasfs02/share_root/yuzhongzhang/Permian/',type,'/')
   } else {
     path='/n/holyscratch01/jacob_lab/yuzhongzhang/Permian/invdata/'
   }
   file=paste0(path,'data_tropomi_',datestr,'.RData')
   if (!file.exists(file)) return(NULL)
   load(file)
   data=data.frame(data,datestr=datestr)
   return(data)
}

#Load processed geos-chem data
load.base.model=function(datestr,type='fullinv3'){
   if (type!='default'){
     path=paste0('/net/seasasfs02/srv/export/seasasfs02/share_root/yuzhongzhang/Permian/',type,'/')
   }else{
     path='/n/holyscratch01/jacob_lab/yuzhongzhang/Permian/invdata/'
   }
   file=paste0(path,'data_model_',datestr,'.RData')
   if (!file.exists(file)) return(NULL)
   load(file)
   d=as.vector(data[1,])
   return(d)
}

#Compute prior XCH4 (ya) if the prior emission differs from that used in GEOS-Chem
#Use emission-based jacobian (k) 
#ya'=ya+k%*%(xa'-xa)
calc.base.model=function(datestr,e,e0){
   path=paste0('/net/seasasfs02/srv/export/seasasfs02/share_root/yuzhongzhang/Permian/fullinv4/')
   file=paste0(path,'data_model_',datestr,'.RData')
   if (!file.exists(file)) return(NULL)
   load(file)
   d=as.vector(data[1,])
   ex=e
   ex$emisall.vec=e$emisall.vec-e0$emisall.vec 
   jac=load.jac(datestr,ex)
   d2=d+apply(jac,1,FUN=sum)
   return(d2)
}

#Load pre-computed emission-based jacobian (script_emis_based_jac.R)
#compute k%*%(xa'-xa)
load.jac=function(datestr,e){
  cmt='fullinv4'
  path='/net/seasasfs02/srv/export/seasasfs02/share_root/yuzhongzhang/Permian/jac_emis_based/'

  file=paste0(path,'data_emis_based_jac_',datestr,'_',cmt,'.RData')
  if (!file.exists(file)) return(NULL)
  load(file)
  id=which(jac<0)
  if (length(id)>0) jac[id]=0
  for (i in 1:nrow(jac)){
    jac[i,]=jac[i,]*e$emisall.vec[i]
  }
  jac=t(jac) 
  jac=cbind(jac,1)
  return(jac)
}

#Load archived prior emission inventory over the Permian Basin
load.prior.emis=function(type='fullinv'){
  if (type=='fullinv'|type=='fullinv4'){
     load('data_ei_mark.RData')
  } else if (type=='fullinv3'){
     load('data_ei_bram.RData')
  } else {
     load(paste0('../ogdistribution/data_ei_mark_',type,'.RData'))
  }
  return(e2)
}

#Compute diagonal elements for prior error matrix
#Assume SA is diagonal
#prior errors for emissions are relative error here
#Now use difference between Bram's and Mark's bottom-up emission
#to represent prior errors
get.sa=function(type,scl=1){
  load('data_ei_mark.RData')
  ea=e2$emisall.vec
  load('data_ei_bram.RData')
  eb=e2$emisall.vec

  if (type=='fullinv4') {
      unc=abs(ea-eb)/ea/sqrt(2)
  } else if (type=='fullinv3'){
      unc=abs(ea-eb)/eb/sqrt(2)
  } else {
      load(paste0('../ogdistribution/data_ei_mark_',type,'.RData'))
      unc=abs(ea-eb)/e2$emisall.vec/sqrt(2)
  }
   unc=unc*scl
   sa=unc^2

   #assume systematic bias over the whole domain as 10 ppbv
   sa=c(sa,10^2)
   
   #test uniform prior error (75%)
   if (scl==0) {
      sa=c(rep(0.75^2, length(unc)),10^2)
   }
   return(sa)
}

#Prepare diagonal elements for observation error matrix
#Assume SO is diagonal
get.so=function(inserr,modelerr,dtrop,l=40){
#   so=inserr^2+modelerr^2
   #so=matrix(0,nrow=nrow(dtrop),ncol=nrow(dtrop))
   #library(geosphere)
   #v=dtrop[c('tlon','tlat')]
   #dis=distm(v, v, fun = distHaversine)/1000. #km
   #so=exp(-dis/l)*modelerr*modelerr
   
   #sodiag=diag(so)
   #so=so*0.9
   #diag(so)=sodiag
   #diag(so)=diag(so)+inserr*inserr
   #diag(so)=modelerr*modelerr+inserr*inserr
   so=rep(modelerr*modelerr+inserr*inserr,nrow(dtrop))
   return(so)
}

#Alternative SO representation
#full So with correlation spatial length of 40 km
get.so.cor=function(inserr,modelerr,dtrop,l=40){
   nn=nrow(dtrop)
   require(Matrix)
   unidate=unique(dtrop$datestr)
   so=Matrix(0,nrow=nn,ncol=nn)
   library(geosphere)
   v=dtrop[c('tlon','tlat')]
   for (iday in unidate){
      id=which(dtrop$datestr==iday)
      dis=distm(v[id,], v[id,], fun = distHaversine)/1000. #km
      so[id,id]=exp(-dis/l)*modelerr*modelerr
      diag(so)[id]=diag(so)[id]+inserr*inserr
   }
   return(so)
}

#Compute distance between two points on a sphere
calc.distance.lon.lat=function(lon1, lat1, lon2, lat2){
  library(geosphere)
  distm(c(lon1, lat1), c(lon2, lat2), fun = distHaversine)/1000. #km
}

#Compute inverse of alternative So by block (day)
get.so.cor.inv=function(inserr,modelerr,dtrop,l=40){
   nn=nrow(dtrop)
   require(Matrix)
   unidate=unique(dtrop$datestr)
   soinv=Matrix(0,nrow=nn,ncol=nn)
   library(geosphere)
   v=dtrop[c('tlon','tlat')]
   for (iday in unidate){
      print(iday)
      id=which(dtrop$datestr==iday)
      dis=distm(v[id,], v[id,], fun = distHaversine)/1000. #km
      so.block=exp(-dis/l)*modelerr*modelerr
      diag(so.block)=diag(so.block)+inserr*inserr
      soinv[id,id]=solve(so.block)
   }
   return(soinv)
}


#Compute analytic inversion 
calc.inverse=function(dtrop,dbase,jac,type,lambda=1,bgprob=0.25, soerr=11, saerr=0.5, arrSo=F){
  m0=dbase
  #remove anomalous area for fitting regional bias
  #id=which((dtrop$lon>=(-106) & dtrop$lon<=(-104) &
  #         dtrop$lat>=32 & dtrop$lat<=35) |
  #         (dtrop$lon>=(-105) & dtrop$lon<=(-102) &
  #          dtrop$lat>=29 & dtrop$lat<=30))
  #if (length(id)>0) jac[id,113]=0

  #fit a mean bias
  #biasbg=quantile(dtrop$sat,probs=0.25) - quantile(m0, probs=0.25)
  #print(paste('bias:',biasbg))
  bias=0
  dif=dtrop$sat-m0-bias

  ####################################################
  #exclude extremely small values, which are spurious
  #id=which(dif<(-20))
  #if (length(id)>0){
  #  print('small')
  #  print(length(id))
  #  dtrop=dtrop[-id,]
  #  dmodel=dmodel[,-id]
  #}
  #m0=dmodel[1,]
  ####################################################

  nobs=nrow(dtrop)
  nstate=ncol(jac)

  sa=get.sa(type,scl=saerr)
  if (!arrSo){
     so=get.so(0, soerr, dtrop)
     so.lam.inv=diag(1./so*lambda)
  } else {
     sigm=4 # ppbv, model error, Cusworth et al.(2018) 
     sigi=0 # instrument error
     if (soerr>sigm) sigi=sqrt(soerr^2-sigm^2)
     so.lam.inv=get.so.cor.inv(sigi,sigm,dtrop) * lambda
  }

  jac.so=t(jac) %*% so.lam.inv 
  KK = jac.so %*% jac
  diag(KK)=diag(KK)+1/sa
  KKT=solve(KK) #Shat

  dif=matrix(dif,nrow=nobs,ncol=1)
  del=jac.so %*% dif
  xhat=KKT%*%del # (adjustment factor - 1) for emissions
                 # ppb CH4 for background bias correction
  #jac[,113]=1
  
  #Reconstruct posterior ya
  yrecon=m0+bias+jac%*%matrix(xhat,nrow=nstate,ncol=1)
  ydata=data.frame(dtrop,mod=m0, y0=m0+bias,yhat=yrecon)
  
  #Extract information for bias correction
  biasxhat=xhat[113]
  biasprec=sqrt(KKT[113,113])
  result=list(xhat=xhat,shat=KKT,nobs=nobs,biasbg25=biasbg,ydata=ydata,soerr=soerr,saerr=saerr,bias=biasxhat,bias.sd=biasprec)
  return(result)
}

#Compile useful data for a day
do.a.day=function(d,m,y,emispri, lonrange=c(-106,-100),latrange=c(29,35),lambda=0.01,type='fullinv',emisfi4=NULL){
  datestr=paste0(sprintf('%04d%02d%02d',y,m,d))
  if (type=='fullinv3' | type=='fullinv4'|type=='default'){
     dtrop=load.tropomi(datestr,type=type)
     dmodel=load.base.model(datestr,type=type)
  } else {
     dtrop=load.tropomi(datestr,type='fullinv4')
     dmodel=calc.base.model(datestr, emispri, emisfi4)
  }
  djac=load.jac(datestr, emispri)
  id=which(
           dtrop$lon>=lonrange[1]& dtrop$lon<=lonrange[2] &
           dtrop$lat>=latrange[1]& dtrop$lat<=latrange[2] )
  if (length(id)<1) return(NULL)
  ntrop=length(id)
  dtrop=dtrop[id,]
  dmodel=dmodel[id]
  djac=djac[id,]
  
  return(list(dtrop=dtrop,dmodel=dmodel,djac=djac))
}

#Compile data for a month and compute the inversion
do.a.month=function(m,y, emispri, lonrange=c(-106,-100),latrange=c(29,35),bgprob=0.25,lambda=0.01,type='fullinv',saerr=0.5,arrSo=F){
  dom=c(31,28,31,30,31,30,31,31,30,31,30,31)
  efi4=NULL
  if (type!='fullinv3' & type!='fullinv4' & type!='default'){
     efi4=load.prior.emis('fullinv4')
  }
  re=lapply(1:dom[m], do.a.day, m=m, y=y, emispri=emispri, 
                   lonrange=lonrange, latrange=latrange, 
                   lambda=lambda, type=type, emisfi4=efi4)
  dtrop=do.call(rbind, lapply(re, function(i) i$dtrop))

  dmodel=do.call(c, lapply(re, function(i) i$dmodel))
  djac=do.call(rbind, lapply(re,function(i) i$djac))
  print(quantile(dtrop$sat-dmodel,probs=seq(0,0.5,0.02)))
  soerr=mean(aggregate(dtrop$sat-dmodel,by=list(dtrop$i,dtrop$j),FUN=sd)$x,na.rm=T)
  print(mean(aggregate(dtrop$sat-dmodel,by=list(dtrop$i,dtrop$j),FUN=sd)$x,na.rm=T))
  re=calc.inverse(dtrop,dmodel,djac, type,lambda=1,bgprob=bgprob, soerr=soerr, saerr=saerr,arrSo=arrSo)
  return(re)
}


################## MAIN ##########################
# 1. base case
#      type='fullinv3'
#      lonrange=c(-106,-100)
#      latrange=c(29,35)
#      saerr=1
#        'fullinv3/data_inv_201805_b0.25_v2.RData'
# 2. change prior
#        'TYPE/data_inv_201805_b0.25_v2.RData'
# 3. change background region
#        'fullinv3/data_inv_201805_bA_v2.RData'
#        A=1: lonrange=c(-108,-98)
#             latrange=c(27,37)
#        A=11: lonrange=c(-110, -96)
#             latrange=c(25,39)
#        A=12: lonrange=c(-110, -100)
#             latrange=c(25,39)
# 4. change prior error
#        'fullinv3/data_inv_201805_bB_v2.RData'
#        B=0.5 saerr=0.5
#        B=1.5 saerr=1.5
#        B=2 saerr=2
#        B=0 saerr=0; use spatial uniform relative error 75%
# 5. change prior error correlation
#        'fullinv3/data_inv_201805_bC_v2.RData'
#        C=20: with off-diagnal elements 40km
# bgprob is now used only as a sensitivity identifier
# saerr is now used as a scaling factor for prior error
###################################################

require(methods)
type='fullinv3'

outpath=paste0('data_inv/',type,'/')

lonrange=c(-106,-100)
latrange=c(29,35)

#change background region
#lonrange=c(-110,-100)
#latrange=c(25,39)
########################
saerr=1
#change prior error
#saerr=0
########################
#use off-diagnol in So?
arrSo=F
#arrSo=T
########################
bgprob=12
lambda=1

args=commandArgs(trailingOnly=T)
months=c(5,6,7,8,9,10,11,12,1,2,3)
years=c(rep(2018,8),rep(2019,3))

if (length(args)==1){
   i=as.integer(args[1])
   m=months[i]
   y=years[i]
}

e=load.prior.emis(type)

print(paste(m,y))
mstr=sprintf('%04d%02d',y,m)
bgstr=as.character(bgprob)
result=do.a.month(m,y,e,lonrange=lonrange,latrange=latrange,bgprob=bgprob,lambda=lambda,type=type, saerr=saerr, arrSo=arrSo)
file=paste0(outpath,'data_inv_',mstr,'_b',bgstr,'_v2.RData')

if (arrSo){
  result$xhat=as.vector(result$xhat)
  result$shat=as.matrix(result$shat)
}
save(result,file=file)


