#!/bin/bash -f

if [ -f groupfile ];
then
  rm groupfile
fi

nrep=`wc temperatures.dat | awk '{print $1}'`
echo $nrep
count=0
for TEMP in `cat temperatures.dat`
do
  let COUNT+=1
  REP=`printf "%03d" $COUNT`
  echo "TEMPERATURE: $TEMP K ==> FILE: equilibrate.mdin.$REP"
  sed "s/XXXXX/$TEMP/g" equilibrate.mdin > temp
  sed "s/RANDOM_NUMBER/$RANDOM/g" temp > equilibrate.mdin.$REP
  echo "-O -rem 0 -i equilibrate.mdin.$REP -o equilibrate.mdout.$REP -c min.rst -r equilibrate.rst.$REP -x equilibrate.mdcrd.$REP -inf equilibrate.mdinfo.$REP -p triad_thf_bent_GR.prmtop " >> equilibrate.groupfile
  
  rm -f temp
done
echo "#" >> groupfile

echo "N REPLICAS  = $nrep"
echo " Done."
