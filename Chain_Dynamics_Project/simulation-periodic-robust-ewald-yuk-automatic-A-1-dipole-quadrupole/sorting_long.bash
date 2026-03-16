
#path='/home/ssingha/PROJECTS/simulation-periodic-quadrupole-robust-sigma-corrected-softer-potential-automatic-A-2/local_simulations_quagga/hack_seed-1-af-0.05-numpar=1024-softer-3-ka-4-m'
nosteps=14801
i=1
while [ $i -le $nosteps ]
do
   filename=A$(printf "%05d" $i).dat
   echo "'$filename'"
   cp $path_from/$filename $path_to/assorted/long/.
   
   i=$[$i+148]
done
