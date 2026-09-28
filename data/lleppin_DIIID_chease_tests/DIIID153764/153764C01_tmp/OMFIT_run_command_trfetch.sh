#!/bin/bash -l
module purge
module load ntcc

var=$(tr_fetch_r9 153764C01 D3D 13)
echo $var
if [[ ! $var =~ "Job retrieved successfully" ]]; then
  exit
fi
while true
do
  if [ -e D3D.*153764C01.tar.gz ]; then
    break
  else
    echo Waiting
    sleep 3
  fi
done
var=D3D.*153764C01.tar.gz
echo $var
tar -xzf $var
echo File untarred
if [ ! -f 153764C01.CDF ]; then
  echo netCDF File NOT Recovered!
  exit
fi
if [ ! -f /fusion/projects/codes/transp/nelsonand/153764C01_tmp//D3D/153764/C01/153764C01.CDF ]; then
  mkdir -p -m 2775 /fusion/projects/codes/transp/nelsonand/153764C01_tmp//D3D/153764
  mkdir -p -m 2775 /fusion/projects/codes/transp/nelsonand/153764C01_tmp//D3D/153764/C01
  echo Results directory created
else
  echo Warning!!!
  echo Results CDF for this run ID already exists!
  echo Will not overwrite netcdf file in results directory.
  exit
fi
cp ./153764C01.CDF /fusion/projects/codes/transp/nelsonand/153764C01_tmp//D3D/153764/C01/
ln -s -f /fusion/projects/codes/transp/nelsonand/153764C01_tmp//153764C01.CDF /fusion/projects/codes/transp/nelsonand/153764C01_tmp//D3D/153764/C01/C01.CDF
if [ -f 153764C01.DATA1 ]; then
  mv ./153764C01.DATA* /fusion/projects/codes/transp/nelsonand/153764C01_tmp//D3D/153764/C01/
fi
echo Run successfully recovered
