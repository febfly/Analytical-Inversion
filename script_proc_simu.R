#Process GEOS-Chem output and match with TROPOMI observations

#Read TROPOMI data and extract useful variables
read.tropomi=function(file){
  #qa_value > 0.5
  #methane_mixing_ratio ppb
  #surface_altidue m
  #methane_profile_apriori mol m-2
  #dry_air_subcolumns mol m-2
  #surface_pressure Pa
  #pressure_interval Pa

  path='/n/holyscratch01/jacob_lab/yuzhongzhang/tropomi/data/'
  require(ncdf4)
  pathlist=c(rep('PRODUCT/',7),rep('DETAILED_RESULTS/',1),
             rep('INPUT_DATA/',6))
  varlist=c('latitude','longitude','qa_value','time_utc',
            'methane_mixing_ratio','methane_mixing_ratio_precision',
            'methane_mixing_ratio_bias_corrected',
            'column_averaging_kernel', 'surface_altitude',
            'methane_profile_apriori','surface_pressure',
            'pressure_interval','dry_air_subcolumns',
            'altitude_levels')
  fulllist=paste0(pathlist,varlist)
  fid=nc_open(paste0(path,file))
  re=lapply(fulllist,function(iv) ncvar_get(fid,iv))
  nc_close(fid)
  names(re)=varlist
  return(re)
}

#Choose TROPOMI data with QA_value>0.5
process.tropomi=function(file){
  a=read.tropomi(file)
  id=which(a$qa_value>0.5)
  idarr=which(a$qa_value>0.5,arr.ind=T)
  dimak=dim(a$column_averaging_kernel)
  redim1=c(dimak[1],dimak[2]*dimak[3])
  redim2=c(dimak[1]+1,dimak[2]*dimak[3])
  
  re=list(
   lat=a$latitude[id],
   lon=a$longitude[id],
   xch4=a$methane_mixing_ratio[id],
   xch4_bias_corrected=a$methane_mixing_ratio_bias_corrected[id],
   xch4_precision=a$methane_mixing_ratio_precision[id],
   ak=array(a$column_averaging_kernel,dim=redim1)[,id],
   dryair=array(a$dry_air_subcolumns,dim=redim1)[,id],
   prior=array(a$methane_profile_apriori,dim=redim1)[,id],
   altitude=array(a$altitude_levels,dim=redim2)[,id],
   surface_pressure=a$surface_pressure[id],
   pressure_interval=a$pressure_interval[id],
   surface_altitude=a$surface_altitude[id]
  )
  return(re)
}

#Read GEOS-FP met data to get GEOS-Chem model surface pressure
#IMIN, IMAX of region    :  34 128
#JMIN, JMAX of region    :  51 121
read.geosfp=function(datestr,irange=c(34,128), jrange=c(51,121)){
  path='/n/holylfs/EXTERNAL_REPOS/GEOS-CHEM/gcgrid/gcdata/ExtData/GEOS_0.25x0.3125_NA/GEOS_FP/'
  yrstr=substr(datestr,1,4)
  mnstr=substr(datestr,5,6)
  file=paste0(path,yrstr,'/',mnstr,'/GEOSFP.',datestr,'.I3.025x03125.NA.nc')
  require(ncdf4)
  fid=nc_open(file)
  ps=ncvar_get(fid,'PS')
  #hPa, 0, 3, 6, 9, 12, 15, 18, 21...
  nc_close(fid)
  ir=seq(irange[1],irange[2])
  jr=seq(jrange[1],jrange[2])
  psx=(ps[ir,jr,8]*2+ps[ir,jr,7])/3.
  return(psx)
}

#Read GEOS-Chem coordinates from its output
read.geoschem.lonlatlev=function(){
  path='/n/holyscratch01/jacob_lab/yuzhongzhang/Permian/Inv/Pert_0000/ts/nc_ts_satellite.20190101.nc'
  require(ncdf4)
  fid=nc_open(path)
  lon=ncvar_get(fid,'LON')
  lat=ncvar_get(fid,'LAT')
  etae=ncvar_get(fid,'ETAE')
  etac=ncvar_get(fid,'ETAC')
  nc_close(fid)
  return(list(lon=lon,lat=lat,etac=etac,etae=etae))
}

#Read CH4 dry mixing ratio fields from GEOS-Chem output
#For a specific date (datestr YYYYMMDD)
read.geoschem=function(datestr,path){
  file=paste0(path,'nc_ts_satellite.',datestr,'.nc')

  require(ncdf4)
  fid=nc_open(file)
  ch4=ncvar_get(fid,'IJ-AVG-S__CH4')
  nc_close(fid)
  return(ch4)
}

#Read CH4 dry mixing ratio fields from GEOS-Chem output
#for a specific date and 
#from baseline and all perturbed simulations
read.geoschem.allsimu=function(datestr,path,simids=seq(0,112)){
  re=lapply(simids, function(i){
      path=paste0(path,'Pert_',sprintf('%04d',i),'/ts/')
      read.geoschem(datestr,path)
  })
  return(re)
}

#Compute model XCH4 from simulated CH4 mixing ratio vertical profile
#First interpolate model profile to tropomi layers
#Second, account for satellite vertical sensitivities (ak) 
#and prior vertical profile (pri)
calc.model.cmr=function(mod, sat){
  modlvl=approx(mod$p, mod$ch4, sat$p, rule=2)$y
  nlvl=length(modlvl)
  modlay=(modlvl[1:(nlvl-1)]+modlvl[2:nlvl])/2.
  totair=sum(sat$air)
  satlay=sat$pri/sat$air*1e9 #ppb
  cmrpri=sum(sat$pri)/totair*1e9
  cmrcor=sum((modlay-satlay)*sat$air/totair*sat$ak)
  cmr=cmrpri+cmrcor
  return(cmr)
}

#Match one TROPOMI observation to GEOS-Chem simulation
#by lon, lat, and date
#Exclude if surface pressure differ by 50 hPa or more
#Compute simulated column methane mixing ratio (XCH4)
process.a.profile=function(ir,tp,gc,var='xch4_bias_corrected'){
  i=which.min(abs(tp$lon[ir]-gc$lon))
  j=which.min(abs(tp$lat[ir]-gc$lat))
  gcps=gc$ps[i,j]
  tpps=tp$surface_pressure[ir]/100 #Pa->hPa
  #disregard when simulation and retrieval surface are off by too much
  if (abs(gcps-tpps)>50) return(NULL) 
  tppint=tp$pressure_interval[ir]/100 #Pa->hPa
  tpak=tp$ak[,ir]
  satellite=list(
    xch4=tp[[var]][ir], #ppb
    xch4prec=tp$xch4_precision[ir],
    ak=tpak,
    pri=tp$prior[,ir], #mol/m2
    air=tp$dryair[,ir], #mol/m2
    p=rev(seq(tpps,by=tppint*(-1.),length.out=length(tpak)+1)) #hPa
  )
 
  modp=gc$etac*gcps
  modcmr=vector('numeric',length=length(gc$ch4))

  for (im in 1:length(gc$ch4)){ 
    model=list(
      p=modp,
      ch4=gc$ch4[[im]][i,j,]
    )
    modcmr[im]=calc.model.cmr(model,satellite)
  }
  result1=data.frame(
     i=i,j=j,lon=gc$lon[i],lat=gc$lat[j],
     tlon=tp$lon[ir], tlat=tp$lat[ir],
     pssat=tpps,psmod=gcps,
     sat=satellite$xch4, prec=satellite$xch4prec )
  result=list(tropomi=result1,model=modcmr)
  return(result)
}

#Process a TROPOMI data file
process.a.file=function(file,gc){
  lonmin=min(gc$lon)
  lonmax=max(gc$lon)
  latmin=min(gc$lat)
  latmax=max(gc$lat)
  tropomi=process.tropomi(file)
  id=which(tropomi$lon>lonmin & tropomi$lon<lonmax &
           tropomi$lat>latmin & tropomi$lat<latmax)
  if (length(id)==0) return(NULL)
  rex=lapply(id, process.a.profile, tp=tropomi, gc=gc)
  
  #postprocess output
  retropo=do.call(rbind,lapply(rex, function(i) i$tropomi))
  xmodel=lapply(rex, function(i) i$model)
  xmodel$along=2
  remodel=do.call(abind,xmodel)
  re=list(tropomi=retropo, model=remodel)
  return(re)
}

#Process all TROPOMI files within a day
#Save corresponding tropomi obs and geos-chem results
#after data filtering and processing
process.a.day=function(date,gclonlatlev,option){
  require(abind)
  datestr=format(date,'%Y%m%d')
  print(datestr)
  files=list.files(
           path='/n/holyscratch01/jacob_lab/yuzhongzhang/tropomi/data/',
           pattern=paste0('S5P_(RPRO|OFFL)_L2__CH4____',datestr,'.*.nc'))
  print(files)
  if (length(files)>0){
     geoschem=read.geoschem.allsimu(datestr,option$gcpath,option$gcsimu)
     ps=read.geosfp(datestr)
     gc=list(
        datestr=datestr,
        lon=gclonlatlev$lon,
        lat=gclonlatlev$lat,
        etac=gclonlatlev$etac,
        etae=gclonlatlev$etae,
        ps=ps,
        ch4=geoschem
    )
  } else {
    return(NULL)
  }
  re=lapply(files,process.a.file, gc=gc)
  
  if (option$save.tropomi){
     data=do.call(rbind,lapply(re,function(i) i$tropomi))
     if (!is.null(data)){
       file=paste0(option$outpath,'data_tropomi_',datestr,'.RData')
       save(data,file=file)
     }
  }
  if (option$save.model){
     model=lapply(re, function(i) i$model)
     model$along=2
     data=do.call(abind,model)
     if (!is.null(data)){
       file=paste0(option$outpath,'data_model_',datestr,'.RData')
       save(data,file=file)
     }
  }
}

###############################################
#Required argument: datestr YYYYMMDD
args=commandArgs(trailingOnly=T)
date=as.Date(args[1])

option=list(
  outpath='/n/holyscratch01/jacob_lab/yuzhongzhang/Permian/invdata/',
  gcpath='/n/holyscratch01/jacob_lab/yuzhongzhang/Permian/Inv/',
  gcsimu=seq(0,112),
  save.tropomi=T,
  save.model=T
)

gclonlatlev=read.geoschem.lonlatlev()
x=process.a.day(date, gclonlatlev, option)

