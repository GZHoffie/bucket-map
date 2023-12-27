#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

filename="$1"

# Check if the file exists
if [ ! -e "$filename" ]; then
    echo "File not found: $filename"
    exit 1
fi

# Loop through each line in the file
while IFS= read -r line; do
    # Check if the line starts with '@'
    if [[ $line == @* ]]; then
        # Replace all spaces with underscores
        modified_line=$(echo "$line" | tr ' ' '_')
        echo "$modified_line"
    else
        # If the line doesn't start with '@', print it as is
        echo "$line"
    fi
done < "$filename"