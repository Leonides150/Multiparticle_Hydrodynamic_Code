#!/bin/bash

rsync -anv report/ /media/ssingha/678A-D3CC/report

echo "---This is a DRY RUN---"
echo "---Want to continue actual rsync? Type yes or no:"
read response
echo "---You selected---" $response
echo "---Initiating rsync---"

if [ $response == 'yes' ]; then
    rsync -av report/ /media/ssingha/678A-D3CC/report
    echo "sync successful"
elif [ $response == 'no' ]; then
     echo "exiting"
fi
