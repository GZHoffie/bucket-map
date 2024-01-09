#!/bin/bash

file_name=$1
output_path=$2

# Delete all 'N' or 'n' in the file
sed '/^>/!s/[Nn]//g;/^$/d' $file_name > $output_path