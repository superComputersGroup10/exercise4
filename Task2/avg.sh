#!/bin/bash 

rm avg_results.csv

for j in {0..6}
do
	current_size=$((2**(6+j)))
	matrix=$current_size'x'$current_size
	echo -e "FOR MATRICES " $matrix >> avg_results.csv
	cat result_$matrix.csv | awk -F'Computation time: ' '{sum1+=$2; } END { print "Avg Computation time: = " sum1/6 }' >> avg_results.csv
	cat result_$matrix.csv | awk -F'MPI time:         ' '{sum2+=$2; } END { print "Avg MPI Time: = " sum2/6 }' >> avg_results.csv
echo -e '\n' >> avg_results.csv
done
