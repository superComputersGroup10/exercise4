#!/usr/bin/python3

import argparse
import re
import csv
import math


def get_data(filename):
    """Returns a list with computing and mpi time
    :filename: File with output results
    :returns: list with computing and mpi time
    """
    with open(filename, 'r') as f:
        text = f.read()
    r_computation_times = re.findall(r'Computation time:\s+(\d+\.\d+)', text)
    r_mpi_times = re.findall(r'MPI time:\s+(\d+\.\d+)', text)
    r_total_times = re.findall(r'Total time:\s+(\d+.\d+)', text)
    r_fence_times = re.findall(r'Fence time:\s+(\d+.\d+)', text)
    r_initialization_time = re.findall(r'Initialization time:\s+(\d+.\d+)', text)

    computation_times = [float(i) for i in r_computation_times]
    mpi_times = [float(i) for i in r_mpi_times]
    total_times = [float(i) for i in r_total_times]
    fence_times = [float(i) for i in r_fence_times]
    initialization_time = [float(i) for i in r_initialization_time]
    return computation_times, mpi_times, total_times, fence_times, initialization_time


def get_mean(data):
    return math.fsum(data) / len(data)


def get_var(data, mean):
    return math.fsum([(mean - val)*(mean - val) for val in data]) / len(data)

def main():
    parser = argparse.ArgumentParser(description="Generates a .csv with the time results")
    parser.add_argument("-f", "--files", nargs = '*', help="List of files with time results")
    parser.add_argument("-o", "--output", help="Output file")
    args = parser.parse_args()

    with open(args.output, "w") as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',')
        spamwriter.writerow(["Case", "Mean MPI", "Mean Comp", "Mean Total", "Mean Fence", "Variance MPI",
            "Variance Comp", "Variance Total", "Variance Fence", "Standard Deviation MPI",
            "Standard Deviation Comp", "Standard Deviation total", "Standar Deviation Fence"])
        for f in args.files:
            computation_times, mpi_times, total_times, fence_times, initialization_time = get_data(f)
            mean_computation = get_mean(computation_times)
            mean_mpi = get_mean(mpi_times)
            mean_total = get_mean(total_times)
            mean_initialization = get_mean(initialization_time)

            var_computation = get_var(computation_times, mean_computation)
            var_mpi = get_var(mpi_times, mean_mpi)
            var_total = get_var(total_times, mean_total)
            var_initialization = get_var(initialization_time, mean_initialization)

            sd_computation = math.sqrt(var_computation)
            sd_mpi = math.sqrt(var_mpi)
            sd_total = math.sqrt(var_total)
            sd_initialization = math.sqrt(var_initialization)

            if len(fence_times) == 0:
                mean_fence = "Not found"
                var_fence = "Not found"
                sd_fence = "Not found"
            else:
                mean_fence = get_mean(fence_times)
                var_fence = get_var(fence_times, mean_fence)
                sd_fence = math.sqrt(var_fence)
            spamwriter.writerow([f, mean_mpi, mean_computation, mean_total, mean_fence, mean_initialization, var_mpi,
                var_computation, var_total, var_initialization, var_fence, sd_mpi, sd_computation, sd_total, 
                sd_initialization, sd_fence])


if __name__ == "__main__":
    main()
