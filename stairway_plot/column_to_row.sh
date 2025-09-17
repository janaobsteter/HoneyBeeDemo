#!/bin/bash

# Check if filename is provided as argument
if [ $# -eq 0 ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

# Extract filename
filename=$1

# Read the file into an array
values=($(< "$filename"))

# Join array elements with space
row=$(IFS=' '; echo "${values[*]}")

# Write the row to a new file
echo "$row" > "${filename%.*}_one_row.txt"

echo "Converted $filename to one row in ${filename%.*}_one_row.txt"
