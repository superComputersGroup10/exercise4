#!/bin/bash 
#@ wall_clock_limit = 00:025:00
#@ job_name = pos-cannon-mpi-ibm
#@ job_type = Parallel
#@ output = cannon_64_$(jobid).out
#@ error = cannon_64_$(jobid).out
#@ class = test
#@ node = 4
#@ total_tasks = 64
#@ node_usage = not_shared
#@ energy_policy_tag = cannon
#@ minimize_time_to_solution = yes
#@ notification = never
#@ island_count = 1
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh

	for k in {0..5}
	do	
		SUPERMUC_PHASE="Phase1"
		CURRENT_SIZE=$((2**(6+$k)))
		MATRIX=$CURRENT_SIZE'x'$CURRENT_SIZE
		WRITE_FILE="$SUPERMUC_PHASE/result_$MATRIX.csv"
		date
		for i in {1..20}
		do
			echo -e "test " $i " on A($MATRIX) and B($MATRIX)" >> $WRITE_FILE
			mpiexec -n 64 ./cannon ../cannon_matrices/$MATRIX-1.in ../cannon_matrices/$MATRIX-2.in  >> $WRITE_FILE
			echo -e '\n' >> $WRITE_FILE
		done
	done

date
for i in {1..20}
do
	echo -e "test " $i " on A(4096x4096) and B(4096x4096)" >> $SUPERMUC_PHASE/result_4096x4096.csv
	mpiexec -n 64 ./cannon ../cannon_matrices/4096x4096-1.in ../cannon_matrices/4096x4096-2.in >> $SUPERMUC_PHASE/result_4096x4096.csv
	echo -e '\n' >> $SUPERMUC_PHASE/result_4096x4096.csv
done

