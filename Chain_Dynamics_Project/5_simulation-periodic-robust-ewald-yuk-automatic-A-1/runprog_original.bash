#!/bin/bash
#
# Sagnik Singha
# Created: Fri Aug 28 15:35:41 CDT 2015
#
#
#
sourcedir=$PWD

#Nx=5 #number of particles in x direction
#Ny=5 #number of particles in x direction
af=0.17
#af=0.22
#0.196
numpar=25 #number of particles
seed_val=1 #seed for the randomize function
ka=5.52
m_yuk=10
#p_sf_x=2 #particle spacing x dir
#p_sf_y=2 #particle spacing y dir
#factor=0.1 #factor to control random displacement of particles
#config=1 #1:square 2:linear 3:rectangle
#time_stepper=1 #1:normal 2:adaptive
dirname="test-tmp"
#dirname="0.7A-check-equivalence-af-${af}-numpar=${numpar}-$ka-ka-$m_yuk-m"
#dirname="hack_af-${af}-numpar=${numpar}-softer-$ka-ka-$m_yuk-m-equivalence-test-2-multiply"
#

cd local_simulations
if [ ! -d  "$dirname" ]; then
  mkdir $dirname
fi

echo SOURCE DIRECTORY:
echo $sourcedir
echo
echo SIMULATION DIRECTORY NAME:
echo $dirname
echo

cd $dirname
echo SIMULATION DIRECTORY:
pwd
echo

#rm *dat

echo COPYING FILES
cp ${sourcedir}/*.mod .
cp ${sourcedir}/*.run .
#ls
echo
echo 'EXECUTING THE SIMULATION'
echo
echo '-----------BASH OUTPUT ENDS-----------'

./main.run<<inputs
$af
$numpar
$seed_val
$ka
$m_yuk
inputs
#
#
#
#
nosteps=$(head frames.dat)
#
echo "###">plotchannel.gplot
echo "set xrange [-5:80]">> plotchannel.gplot
echo "set yrange [-5:35]">> plotchannel.gplot
echo "set size ratio -1">> plotchannel.gplot
echo "set isosample 10000">> plotchannel.gplot
echo "set nokey">> plotchannel.gplot

i=1
while [ $i -le $nosteps ]
do
   filename=A$(printf "%05d" $i).dat
   echo "plot '$filename' u 2:3:4 with circles">> plotchannel.gplot
   echo "pause 0">> plotchannel.gplot
   i=$[$i+1]
done
#gnuplot -geometry 900x900 plotchannel.gplot
