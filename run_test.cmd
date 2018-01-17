#!/bin/csh
#$ -cwd
#$ -o /dev/null
#$ -e /dev/null
#$ -l h_data=1G,h_rt=24:00:00
./test \
12 \
5 \
3000
