#!/bin/bash

# Define the number of iterations
num_iterations=571

# Define the file paths
csv_dir="../../../data/mrmr/random"
output_dir="../../../data/mrmr/random/output"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Run mRMR for each iteration
for ((i=356; i<=num_iterations; i++)); do

    # Run mRMR with the shuffled CSV file and save the output
    ./mrmr -i "$csv_dir/copd_ctrl_train_$i.csv" -n 500 -t 1 -v 16235 -s 225 > "$output_dir/output_$i.txt"
done
