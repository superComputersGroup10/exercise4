#!/bin/bash 
#@ wall_clock_limit = 00:25:00
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

	for k in {0..6}
	do	
		SUPERMUC_PHASE="Phase1"
		CURRENT_SIZE=$((2**(6+$k)))
		MATRIX=$CURRENT_SIZE'x'$CURRENT_SIZE
		WRITE_FILE="$SUPERMUC_PHASE/result_$MATRIX.csv"
		date
		for i in {1..20}
		do
			echo -e "test " $i " on A($MATRIX) and B($MATRIX)" >> $WRITE_FILE
			mpiexec -n 64 ./cannon ../cannon_matrices_bin/$MATRIX-1.in ../cannon_matrices_bin/$MATRIX-2.in /dev/null >> $WRITE_FILE
			echo -e '\n' >> $WRITE_FILE
		done
	done

