#!/bin/bash
KVAL=18
NUMMSG=2000
NUMSEEDS=100
STARTSEED=1
sed -i "7s/[0-9][0-9]*/$KVAL/" run_test.cmd
sed -i "9s/[0-9][0-9]*/$NUMMSG/" run_test.cmd
#sleep 0.2
for ((SEED=$STARTSEED; SEED<=$STARTSEED+$NUMSEEDS; SEED++))
do
	sed -i "8s/[0-9][0-9]*/$SEED/" run_test.cmd
	sleep 0.1
	qsub run_test.cmd
done
