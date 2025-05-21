#!/bin/bash

# Check if correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_file.tsv> <output_directory> <number_of_chunks>"
    exit 1
fi

# Assigning arguments to variables
input_file=$1
output_dir=$2
num_chunks=$3

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Determine number of lines in input file
total_lines=$(wc -l < "$input_file")

# Calculate approximate lines per chunk
lines_per_chunk=$(( ($total_lines + $num_chunks - 1) / $num_chunks ))

# Split the input file into chunks
split -l $lines_per_chunk "$input_file" "$output_dir/chunk"

# Prepend header to each chunk
for file in "$output_dir"/*; do
    # Get the header from the original file
    header=$(head -n 1 "$input_file")
    # Prepend the header to the chunk file
    sed -i "1i $header" "$file"
done

echo "File has been chunked into $num_chunks files in $output_dir"
