#Convert relative jacobian matrix (ppbv CH4/100% perturbation of emission)
#to absolute jacobian matrix (ppbv CH4/(kg/m2/s) change of emission)

#Load processed and archived GEOS-Chem results (from script_proc_simu.R)
load.model=function(datestr,type='fullinv3'){
   path='/n/holyscratch01/jacob_lab/yuzhongzhang/Permian/invdata/'
   file=paste0(path,'data_model_',datestr,'.RData')
   if (!file.exists(file)) return(NULL)
   load(file)
   return(data)
}


#Load archived prior emission inventory over the Permian Basin
load.prior.emis=function(type='fullinv'){
  if (type=='fullinv'|type=='fullinv4'){
  load('data_ei_mark.RData')
  } else if (type=='fullinv3'){
  load('data_ei_bram.RData')
  } else if (type=='default'){
  load('data_ei_default.RData')
  } else if (type=='default2'){
  load('data_ei_default2.RData')
  }
  return(e2)
}

#Compute emission based jacobian from model output
#d: matrix containing all GEOS-chem results
#e: list containing emission inventory used in the actual GEOS-Chem simulation
calc.emis.based.jac=function(d,e){
  dm=dim(d)
  jac=matrix(0,nrow=dm[1]-1,ncol=dm[2])
  for (i in 1:(dm[1]-1)){
      print(e$emisall.vec[i])
      if (e$emisall.vec[i]==0) {
         jac[i,]=0
      } else {
         jac[i,]=(d[i+1,]-d[1,])/e$emisall.vec[i]*2 #2->50% perturbation
      }
  }
  return(jac)
}

#Load model result and compute emission-based jacobian
do.calc=function(date,e,cmt,path){
  print(date)
  m=load.model(date,type=cmt)
  if (is.null(m)) return(NULL)
  jac=calc.emis.based.jac(m,e) #ppb/(kg/m2/s)
  file=paste0(path,'data_emis_based_jac_',date,'_',cmt,'.RData')
  save(jac, file=file)
}

###################################################
cmt='default2' 
#Inventory used in GEOS-Chem baseline simulation

if (cmt!='default'&cmt!='default2'){
  pathout='/net/seasasfs02/srv/export/seasasfs02/share_root/yuzhongzhang/Permian/jac_emis_based/'
}else {
  pathout='/n/holyscratch01/jacob_lab/yuzhongzhang/Permian/invdata/'
}

e=load.prior.emis(type=cmt)
dates=seq(from=as.Date('2019-03-28'),
            to=as.Date('2019-03-31'),
            by='day')

datestrs=sapply(dates,format, format='%Y%m%d')
x=lapply(datestrs, do.calc, e=e, cmt=cmt, path=pathout)
