#!/bin/bash

# Define the list of base names (without extensions)
samples=(
"TCL_Bcell3660_2_Counts"
"TCL_Bcell3662_4_Counts"
"TCL_Bcell3666_6_Counts"
"TCL_Bcell3668_8_Counts"
"TCL_Bcell3670_10_Counts"
"TCL_Bcell3672_12_Counts"
"TCL_Bcell3658_1_Counts"
"TCL_Bcell3661_3_Counts"
"TCL_Bcell3665_5_Counts"
"TCL_Bcell3667_7_Counts"
"TCL_Bcell3669_9_Counts"
"TCL_Bcell3671_11_Counts")

# Loop through each sample and generate scatter plots
for sample in "${samples[@]}"; do
    # Define the .cnr and .cns file names
    cnr_file="${sample}.cnr"
    cns_file="${sample}.cns"

    # Check if both files exist
    if [[ -f "$cnr_file" && -f "$cns_file" ]]; then
        echo "Generating scatter plot for $sample..."
        
        # Run the cnvkit.py scatter command
        #cnvkit.py scatter -s "$cns_file" "$cnr_file" -o "${sample}_scatter.pdf"
        cnvkit.py call --drop-low-coverage  -m clonal --purity 1 --sample-sex Female  $cns_file -o "${sample}_call.cns"
        # Check if the scatter plot was successfully generated
        if [[ $? -eq 0 ]]; then
            echo "Successfully generated scatter plot for $sample"
        else
            echo "Error generating scatter plot for $sample"
        fi
    else
        echo "Missing files for $sample. Skipping..."
    fi
done

echo "All scatter plots processed!"

