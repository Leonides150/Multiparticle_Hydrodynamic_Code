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

#af=0.3
Lx=23.04
#89.679859462366863
#36.611649314550760
#23.04
Ly=3.84
#89.679859462366863
#36.611649314550760
#3.84
numpar=2 #number of particles
seed_val=1 #seed for the randomize function
ka=6
m_yuk=8
#Aswap=3.325
#3.55

drho[1]=1.35d0
#4.
drho[2]=0.


#dirname="two_particle_swap_traj_drho1=${drho[1]}_drho2=${drho[2]}_Aswap=${Aswap}_L-${Lx}_testing"
#dirname="quadrupole_2_particles_Lx=23.04_Ly=3.84_ini_xsep_4_n2m10_yukoff_shortdistancefixed_ewald_derives_stop_test_position_velocity_matching"
#dirname="two_particle-ka-${ka}-m_yuk-${m_yuk}-Lx_Ly-${Lx}_${Ly}-swap"
#dirname="to_check_order_magnitude_ewald_sum_near_contact"
#dirname="modify_range_of_yukawa"
#dirname="dipolar_two_particle_large_cell_stopping_distance_2"
dirname="quadrupolar_two_particle_stopping_distance_Quals-2"


cd local_simulations_dobara
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
cp $sourcedir/*.f90 .
cp $sourcedir/runprog2.bash .

echo "###">test.txt
echo "af=$af">>test.txt
echo "ka=$ka">>test.txt
echo "numpar=$numpar">>test.txt
echo "m_yuk=$m_yuk">>test.txt
echo 'TEXT FILE CREATED'
#ls
echo
echo 'EXECUTING THE SIMULATION'
echo
echo '-----------BASH OUTPUT ENDS-----------'

./main.run<<inputs
$Lx
$Ly
$numpar
$seed_val
$ka
$m_yuk
$Aswap
${drho[1]}
${drho[2]}
inputs
#
#
#
#
nosteps=$(head frames.dat)
#
#touch unit_call_val.txt

# echo "###">> plotchannel_jb.gplot
# echo "">> plotchannel_jb.gplot
# echo "set term png font "dejavu/DejaVuSerif-Bold.ttf" 26 size 1000,1000 enhanced">> plotchannel_jb.gplot
# echo "">> plotchannel_jb.gplot
# echo"l=126.82647352339235">> plotchannel_jb.gplot
# echo"frac=0.25">> plotchannel_jb.gplot
# echo"">> plotchannel_jb.gplot
# echo "set xrange [-frac*l:(1.+frac)*l]">> plotchannel_jb.gplot
# echo "set yrange [-frac*l:(1.+frac)*l]">> plotchannel_jb.gplot
# echo "">> plotchannel_jb.gplot
# echo "set size ratio -1">> plotchannel_jb.gplot
# echo "set isosample 10000">> plotchannel_jb.gplot
# echo "set nokey">> plotchannel_jb.gplot
# echo"">> plotchannel_jb.gplot
# echo"rad=1">> plotchannel_jb.gplot
# echo"">> plotchannel_jb.gplot
# echo"#BOX">> plotchannel_jb.gplot
# echo"set object 1 rect from 0,0 to l,l">> plotchannel_jb.gplot
# echo"">> plotchannel_jb.gplot
# echo"imax=2750">> plotchannel_jb.gplot
# echo"">> plotchannel_jb.gplot
# echo"incr=4">> plotchannel_jb.gplot
# echo"">> plotchannel_jb.gplot
# echo"do for [i=1:imax]{">> plotchannel_jb.gplot
# echo"">> plotchannel_jb.gplot
# echo"set label "1024  particles-Periodic Quadrupolar Interaction"  at graph 0,1.02 left">> plotchannel_jb.gplot
# echo"">> plotchannel_jb.gplot
# echo"j=(i-1)*incr+1">> plotchannel_jb.gplot
# echo"set title sprintf("j=%i",j) font "Times-Roman,24"">> plotchannel_jb.gplot
# echo"">> plotchannel_jb.gplot
# echo"imagefile='M'.sprintf("%5.5i",i).'.png'">> plotchannel_jb.gplot
# echo"datafile='A'.sprintf("%5.5i",j).'.dat'">> plotchannel_jb.gplot
# echo"">> plotchannel_jb.gplot
# echo"print i">> plotchannel_jb.gplot
# echo"print imagefile">> plotchannel_jb.gplot
# echo"print datafile">> plotchannel_jb.gplot
# echo"">> plotchannel_jb.gplot
# echo"set output imagefile">> plotchannel_jb.gplot
# echo"plot datafile u 2:3:(rad*$6) with circles lc rgb "black", \"
# echo"datafile u ($2-$7):($3-$7):(rad*$6) with circles lc rgb "red", \"
# echo"datafile u ($2-$7):3:(rad*$6) with circles lc rgb "red", \"
# echo"datafile u ($2-$7):($3+$7):(rad*$6) with circles lc rgb "red", \"
# echo"datafile u 2:($3-$7):(rad*$6) with circles lc rgb "red", \"
# echo"datafile u 2:($3+$7):(rad*$6) with circles lc rgb "red", \"
# echo"datafile u ($2+$7):($3-$7):(rad*$6) with circles lc rgb "red", \"
# echo"datafile u ($2+$7):3:(rad*$6) with circles lc rgb "red", \"
# echo"datafile u ($2+$7):($3+$7):(rad*$6) with circles lc rgb "red""
# echo"pause 0"
# echo"}"
# echo""

# i=1
# while [ $i -le $nosteps ]
# do
#    filename=A$(printf "%05d" $i).dat
#    echo "plot '$filename' u 2:3:4 with circles">> plotchannel.gplot
#    echo "pause 0">> plotchannel.gplot
#    i=$[$i+1]
# done
# gnuplot -geometry 900x900 plotchannel.gplot
