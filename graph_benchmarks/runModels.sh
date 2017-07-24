#!/bin/bash
export KMP_AFFINITY=verbose,granularity=core,compact
ls models/ | while read line
do
	echo $line
	echo clean
	#./bin/benchmark jaccard degree models/$line > results/jac-deg-$line
	#./bin/benchmark jaccard none models/$line > results/jac-non-$line
	./bin/benchmark cluster_coeff degree models/$line > results/cc-non-$line
	#./bin/benchmark cluster-coeff degree models/$line > results/cc-deg-$line
done
