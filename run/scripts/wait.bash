#!/bin/bash

#BSUB -J Wait
#BSUB -q normal
#BSUB -n 1
#BSUB -R rusage[mem=1G]
#BSUB -R span[ptile=1]

echo hello wait for 1 seconds
sleep 1



