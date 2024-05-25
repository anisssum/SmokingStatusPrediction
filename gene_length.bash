#!/bin/bash

# Use awk to extract the required columns from the tables and merge them
awk 'NR>1 {print $2}' data/gene_annotation_table.txt > data/temp1.txt
awk '{print $6}' data/human.gene_sums.R109.gtf > data/temp2.txt

# Merge the data into one table
paste data/temp1.txt data/temp2.txt > data/gene_length.txt

# Remove temporary files
rm data/temp1.txt data/temp2.txt