#!/bin/bash
#
# Sagnik Singha
# Created: Tue Mar 26 14:47:10 CDT 2019
#
#
#
sourcedir=$PWD

#Nx=5 #number of particles in x direction
#Ny=5 #number of particles in x direction

Lx=150.
#225
#75.61
#150
#75.61
#100
### 769.
### 60.
#150
#69.12
#100
#1000
Ly=15.38
#7.69
### 60.
#150
#17.28
#100
#1000
dt=0.05
### 0.05
### 0.025
numpar=200  ### 20  #44 #number of particles after numpar_extra ###so Total base particles are assumed to be 50
#i.e. numpar-numpar_extra=50
seed_val=1 #seed for the randomize function
ka=6
m_yuk=8
#Aswap=3.325
#3.55

### drho[1]=3.
#1.5
#2.
### drho[2]=0.
# pair_num=10 #number of particles
#pair_num = no. of particles pairs from the centre you want to evacuate to fill newer particles
# numpar_extra=10 #number of particles to be added
#i.e. numpar-numpar_extra=50

g_stren=0.05
#dirname="two_particle_swap_traj_drho1=${drho[1]}_drho2=${drho[2]}_Aswap=${Aswap}_L-${Lx}_testing"
#dirname="quadrupole_2_particles_Lx=23.04_Ly=3.84_ini_xsep_4_n2m10_yukoff_shortdistancefixed_ewald_derives_stop_test_position_velocity_matching"
#dirname="two_particle-ka-${ka}-m_yuk-${m_yuk}-Lx_Ly-${Lx}_${Ly}-swap"
#dirname="to_check_order_magnitude_ewald_sum_near_contact"
#dirname="modify_range_of_yukawa"
##dirname="follow_up_varying_linear_56_periodic_particle_3"
#dirname="linear_array_diffusion_quadrupoles_match"
#dirname="dipolar_quadrupolar_compressed_periodic_array_test_numpar-${numpar}_dx-${drho[1]}_Lx-${Lx}"
#dirname="dipolar_quadrupolar_swap_linear_array_sinusoidal_displacement_numpar-${numpar}_3"
### dirname="${numpar}-particle_dipolar-quadrupolar_2-swap"
###dirname="${numpar}-particle_${Lx}-Lx_${Ly}_Ly_${dt}-dt_gstren-${g_stren}_cell-frac0.5dipolar_self-term-0_side_wall-front_vel-trailing_end-final-production_2"
dirname="${numpar}-particle_${Lx}-Lx_${Ly}_Ly_${dt}-dt_gstren-${g_stren}_cell-frac0.5dipolar_self-term-0_side_wall-front_vel-trailing_end_single-piston-final"
###dirname="${numpar}-particle_${Lx}-Lx_${Ly}_Ly_${dt}-gstren-${g_stren}-dipolar_chain_investigations_4_4"
#dirname="${numpar}-particle_${Lx}-Lx_${Ly}_Ly_${dt}-dt_dipolar_noswap_helf-cell_self-term-0_9"
### dirname="test_1"

cd local_simulations_baiji
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
cp $sourcedir/runprog4_varying_linear.bash .

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
runname=main_${numpar}_${Lx}_${Ly}_${dt}_${g_stren}_frac_0.5
mv main.run ${runname}.run
echo $runname

./${runname}.run<<inputs
$Lx
$Ly
$dt
$numpar
$seed_val
$ka
$m_yuk
$g_stren
inputs
#
#
#
#
nosteps=$(head frames.dat)
#
