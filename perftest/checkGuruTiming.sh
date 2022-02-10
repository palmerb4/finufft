#!/bin/bash
# by Andrea Malleo, summer 2019.

srcpts=1e7
tolerance=1e-6
debug=1
modes[0]=1e6
modes[1]=1
modes[2]=1
modes[3]=1
modes[4]=1

modes[5]=1e3
modes[6]=1e3
modes[7]=1
modes[8]=1
modes[9]=1

modes[10]=1e2
modes[11]=1e2
modes[12]=1e2
modes[13]=1
modes[14]=1

modes[15]=20
modes[16]=20
modes[17]=20
modes[18]=20
modes[19]=1

modes[20]=10
modes[21]=10
modes[22]=10
modes[23]=10
modes[24]=10

for dimension in 1 2 3 4 5
do
    for type in 1 2 3
    do
	for n_trials in 1 20 41
	do
	    declare -i row
	    row=${dimension}-1

	    declare -i index
	    index=row*3

	    declare -i modeNum
	    modeNum1=${modes[index]}
	    modeNum2=${modes[index+1]}
	    modeNum3=${modes[index+2]}
	    modeNum4=${modes[index+3]}
	    modeNum5=${modes[index+4]}

	    echo "./guru_timing_test ${n_trials} ${type} ${dimension} ${modeNum1} ${modeNum2} ${modeNum3} ${modeNum4} ${modeNum5} ${srcpts} ${tolerance} ${debug}"				 
	    ./guru_timing_test ${n_trials} ${type} ${dimension} ${modeNum1} ${modeNum2} ${modeNum3} ${modeNum4} ${modeNum5} ${srcpts} ${tolerance} ${debug}
	done
    done
done

			       
