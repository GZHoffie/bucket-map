#!/bin/bash

# Input SAM file
input_file="./output/graphmap_map.sam"

# Output SAM file
output_file="./output/graphmap_adjusted_map.sam"

# Use sed to remove strings within lines starting with "ZE:f"
sed 's/\tZE:f:[^\t]*//g' "$input_file" > "$output_file"

echo "Done! Result saved in $output_file"