#!/bin/bash

# Check if filename is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

filename=$1

# Check if file exists
if [ ! -f "$filename" ]; then
    echo "File not found: $filename"
    exit 1
fi

# Process the file using awk
awk '{
    for (i = 2; i <= NF; i+=2) {
        printf "%s|%s\t", $i, $(i+1)
    }
    printf "\n"
}' "$filename"
