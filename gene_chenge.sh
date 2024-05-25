#!/bin/bash

# Input table
input_table="tcga.gene_sums.LUAD.R109"

# Second table with gene ID mappings
mapping_table="gene_annotation_table.txt"

# Output table with replaced values
output_table="tcga.gene_sums_refseq.csv"

# Create or clear the output table if it already exists
echo "gene_id\tgene_name" > "$output_table"

# Process the input table
while IFS= read -r line; do
    # Extract gene_id from each line
    gene_id=$(echo "$line" | grep -o -P 'gene_id "[^"]+"' | sed 's/gene_id "\(.*\)"/\1/')

    # If gene_id is found, replace it with the corresponding value from the mapping table
    if [ -n "$gene_id" ]; then
        replacement=$(grep "$gene_id" "$mapping_table" | awk '{print $2}')
        if [ -n "$replacement" ]; then
            line=$(echo "$line" | sed "s/$gene_id/$replacement/")
        fi
    fi

    # Append the line to the output table
    echo -e "$line" >> "$output_table"
done < "$input_table"