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

Lx=60.
#60.
#33.
#60.
#75
#23.04
#20
#100
#1000
Ly=60.
#60.
#33.
#60.
#75
#3.84
#20
#100
#1000
numpar=20
#20
#11
#20 #number of particles after numpar_extra ###so Total base particles are assumed to be 50
#i.e. numpar-numpar_extra=50
seed_val=1 #seed for the randomize function
ka=6
m_yuk=8
#Aswap=3.325
#3.55

drho[1]=3.
#2.5
#2.
drho[2]=0.
#pair_num=10 #number of particles
#pair_num = no. of particles pairs from the centre you want to evacuate to fill newer particles
#numpar_extra=18 #number of particles to be added
#i.e. numpar-numpar_extra=50


#dirname="two_particle_swap_traj_drho1=${drho[1]}_drho2=${drho[2]}_Aswap=${Aswap}_L-${Lx}_testing"
#dirname="quadrupole_2_particles_Lx=23.04_Ly=3.84_ini_xsep_4_n2m10_yukoff_shortdistancefixed_ewald_derives_stop_test_position_velocity_matching"
#dirname="two_particle-ka-${ka}-m_yuk-${m_yuk}-Lx_Ly-${Lx}_${Ly}-swap"
#dirname="to_check_order_magnitude_ewald_sum_near_contact"
#dirname="modify_range_of_yukawa"
##dirname="follow_up_varying_linear_56_periodic_particle_3"
#dirname="linear_array_diffusion_quadrupoles_match"
#dirname="dipolar_periodic_array_test_numpar-${numpar}_dx-${drho[1]}_Lx-${Lx}"
#dirname="quadrupolar_two_particle_stopping_distance_Quals-4_rectangle"
# dirname="linear_array_sinusoidal_displacement_numpar-${numpar}_10"
#dirname="di_correction_2_perturb_sim_Lx-${Lx}_low_amp_1_freq3_phase_shift90_full_sinusoid_even-pairs"
#dirname="${numpar}-larger_sys_smaller_amp_perturb-more_waves-6"
dirname="${numpar}-particle_dipolar_2-swap"

cd local_simulations_dobara/testing_zones_quadrupolar_swapping
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

./main.run<<inputs
$Lx
$Ly
$numpar
$seed_val
$ka
$m_yuk
${drho[1]}
${drho[2]}

inputs
#
#
#
#
nosteps=$(head frames.dat)
#
