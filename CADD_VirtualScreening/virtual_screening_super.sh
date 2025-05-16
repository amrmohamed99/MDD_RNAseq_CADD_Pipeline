#!/usr/bin/bash

echo "========================================================"
echo "Welcome, Amr Al-Hefnawy!"
echo "AAlhfnawy@nu.edu.eg"
echo "======================================================"
echo "Virtual Screening using Vina"
echo "======================================================"

echo "Paste the folder destination containing pdbqt files: "

read -r folder
cd $folder || exit  

# Loop through each file with .pdbqt extension in the folder
for f in *.pdbqt
do
    echo "Processing $f"

    # Get base name of the file (without the .pdbqt extension)
    base_name="${f%.pdbqt}"

    # Run AutoDock Vina and output uniquely named files for each ligand
    /home/amr/CADD/autodock_vina_1_1_2_linux_x86/bin/vina \
        --config ../config.txt \
        --ligand "$folder/$f" \
        --out "$folder/${base_name}_out.pdbqt" \
        --log "$folder/${base_name}_log.log"

done

