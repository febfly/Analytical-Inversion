#!/bin/bash

# ------------------------------------------------------------------
### Make run directories for Jacobian run
# ------------------------------------------------------------------

##=======================================================================
## 1. Set variables
START_I=1
END_I=281
pPERT="1.5"
BFname=submit_run_template
InputFile="input.geos.pert"
nameB="Inv"
mcwd=`pwd`
geos="gcclassic"
fclust="permian_clust_281.nc"

date1=20190101
date2=20200101
   
### Number of Clusters
start=$START_I
stop=$END_I
x=$start

### Define the base directory
baseDir="run_template"

##=======================================================================
## 2. Create run directories

while [ $x -le $stop ];do

   ### Positive or negative perturbation
   if [ $x -eq 0 ]; then
      PERT="1.0"
      xUSE=1
   elif [ $x -eq -1 ]; then
      PERT="1.0"
      xUSE=$x
   else
      PERT=$pPERT
      xUSE=$x
   fi

   ### Add zeros to string name
   if [ $x -lt 10 ]; then
      xstr="000${x}"
   elif [ $x -lt 100 ]; then
      xstr="00${x}"
   elif [ $x -lt 1000 ]; then
      xstr="0${x}"
   else
      xstr="${x}"
   fi

   ### Define the name
   name="${nameB}_${xstr}"

   ### Make the directory
   runDir="../${name}"
   mkdir -p ${runDir}
      
   ### Copy and point to the necessary data
   cp -r ${baseDir}/*   ${runDir}

   ### Create input.geos file from template

   cd $runDir
#   sed -e "s:11111111:${date1}:g" \
#       -e "s:22222222:${date2}:g" \
   sed -e "s:pertnumpertnum:${xUSE}:g" \
       -e "s:pertfacpertfac:${PERT}:g" \
       $InputFile > input.geos
   rm $InputFile

   ### Create run script from template
   sed -e "s:namename:${name}:g" \
       ${BFname} > ${name}.run
       chmod 755 ${runDir}/${name}.run
   rm ${BFname}

   ### Create symbolic links to executable and cluster file
   ln -s ${mcwd}/bin/${geos} ${runDir}/
   ln -s ${mcwd}/$fclust ${runDir}/$fclust
   ln -s ~/zhang/data_geoschem/GEOSChem_BC_tropomi BC
   cd $mcwd

   ### Increment
   x=$[$x+1]

   ### Print diagnostics
   echo "CREATED: ../${name}/"

done

exit 0
