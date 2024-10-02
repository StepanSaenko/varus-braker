#!/bin/bash

# Function to process a single genome.fa file
process_genome_file() {
    local file=$1

    # Extract headers and truncate them
    headers=$(grep "^>" "$file" | cut -d' ' -f1)

    # Get lengths of headers and sort them
    header_lengths=$(echo "$headers" | awk '{ print length($0) " " $0 }' | sort -n)

    # Check if there are exactly two shortest headers
    shortest_headers=$(echo "$header_lengths" | awk '{print $1}' | uniq -c | awk '$1 == 2 && NR == 1 { print $2 }')

    if [ ! -z "$shortest_headers" ]; then
        # Get the two headers with the shortest length
        two_shortest=$(echo "$header_lengths" | awk -v len="$shortest_headers" '$1 == len { print $2 }')
        echo "$two_shortest"
    else
        echo ">MT"
    fi
}

# Check if a file is provided as an argument
if [ $# -ne 1 ]; then
    echo "Usage: $0 <genome.fa>"
    exit 1
fi

# Process the input genome file
process_genome_file "$1"
