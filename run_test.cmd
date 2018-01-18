#!/bin/csh
#$ -cwd
#$ -o /dev/null
#$ -e /dev/null
#$ -l h_data=2G,h_rt=24:00:00
./test2 \
19 \
2 \
2 \
1e-3
