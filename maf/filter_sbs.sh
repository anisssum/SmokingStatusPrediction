#!/bin/bash

# Check if a filename argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

# Check if the file exists
if [ ! -f "$1" ]; then
    echo "File $1 not found"
    exit 1
fi

# Read the file line by line and print only those lines
# where the values in the 4th and 5th columns are the same
awk 'BEGIN{FS=OFS="\t"} NR==1 || $4 == $5 {print}' "$1"