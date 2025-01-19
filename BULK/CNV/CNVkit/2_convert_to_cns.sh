# Define the list of .cnr files
files=(
"TCL_Bcell3660_2_Counts.cnr"
"TCL_Bcell3662_4_Counts.cnr"
"TCL_Bcell3666_6_Counts.cnr"
"TCL_Bcell3668_8_Counts.cnr"
"TCL_Bcell3670_10_Counts.cnr"
"TCL_Bcell3672_12_Counts.cnr"
"TCL_Bcell3658_1_Counts.cnr"
"TCL_Bcell3661_3_Counts.cnr"
"TCL_Bcell3665_5_Counts.cnr"
"TCL_Bcell3667_7_Counts.cnr"
"TCL_Bcell3669_9_Counts.cnr"
"TCL_Bcell3671_11_Counts.cnr"
)

# Loop through each file and convert it to .cns
for file in "${files[@]}"; do
    # Extract the base name without extension
    base_name="${file%.cnr}"

    # Run the cnvkit.py segment command
    echo "Processing $file..."
    cnvkit.py segment -m cbs "$file" -o "${base_name}.cns"

    # Check if the conversion was successful
    if [[ $? -eq 0 ]]; then
        echo "Successfully converted $file to ${base_name}.cns"
    else
        echo "Error converting $file"
    fi
done

echo "All files processed!"
