#!/bin/bash
# Sagnik Singha
# Created: Tue Mar 26 14:47:10 CDT 2019
#
#
#
sourcedir=$PWD

#Nx=5 #number of particles in x direction
#Ny=5 #number of particles in x direction

Lx=48.
Ly=48.

dt=0.05

numpar=1 #number of particles after numpar_extra ###so Total base particles are assumed to be 50
alpha=1.
seed_val=1 #seed for the randomize function
ka=6
m_yuk=8
#Aswap=3.325

drho[1]=3.
drho[2]=0.

drho_2[1]=3.
#1.5
drho_2[2]=0.

num_high=8

dirname="${numpar}_dx-${drho[1]}_Lx-${Lx}_Ly_${Ly}_step_alpha-${alpha}_step_high-density_${num_high}"
#dirname="test_step"

cd local_simulations_chautha/sims_for_paper_fig_PRL 
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
cp $sourcedir/runprog4_varying_linear_backup_2_step_initial.bash .

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
runname=main_${numpar}_dt-${dt}_Lx_${Lx}_dx_${drho[1]}_${alpha}_step${num_high}
mv main.run ${runname}.run
ls $runname

./${runname}.run<<inputs
$Lx
$Ly
$dt
$numpar
$alpha
$seed_val
$ka
$m_yuk
${drho[1]}
${drho[2]}
${drho_2[1]}
${drho_2[2]}
$num_high

inputs
#
#
#
#
nosteps=$(head frames.dat)
#
