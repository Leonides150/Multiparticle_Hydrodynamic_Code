#!/bin/bash

rsync -anv /media/ssingha/678A-D3CC/report/ report

echo "---This is a DRY RUN---"
echo "---Want to continue actual rsync? Type yes or no:"
read response
echo "---You selected---" $response
echo "---Initiating rsync---"

if [ $response == 'yes' ]; then
    rsync -av /media/ssingha/678A-D3CC/report/ report
    echo "sync successful"
elif [ $response == 'no' ]; then
     echo "exiting"
fi
     
#rsync -ruvva /media/ssingha/678A-D3CC/report /home/ssingha/PROJECTS/simulation-periodic-robust-ewald-yuk-automatic-A-1/compare_ewald_swapping_trajectory_from_fortran/quals/.
