#!/bin/bash

START_I=-1
END_I=2
nameB="InvOH"
mcwd=`pwd`

x=$START_I
stop=$END_I

while [ $x -le $stop ];do
   if [ $x -eq -1 ]; then
      xstr="truth"
   elif [ $x -lt 10 ]; then
      xstr="000${x}"
   elif [ $x -lt 100 ]; then
      xstr="00${x}"
   elif [ $x -lt 1000 ]; then
      xstr="0${x}"
   else
      xstr="${x}"
   fi

   name="${nameB}_${xstr}"
   dir="${mcwd}/run_dirs/${name}"

   cp ${mcwd}/convert.pro ${dir}
   cp ${mcwd}/convert.idl ${dir}
   cd ${dir}
   idl convert.idl
   echo "done ${dir}" 

   x=$[$x+1]
done

exit 0 
