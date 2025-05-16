#!/bin/bash

echo "========================================================"
echo "Welcome, Amr Al-Hefnawy!"
echo "AAlhfnawy@nu.edu.eg"
echo "======================================================"
echo "RDKit Molecule Minimization Script Initialized"
echo "========================================================"
echo

# Ask user for directory
read -p "Enter the directory containing molecule files (.sdf or .pdb): " mol_dir
if [[ ! -d "$mol_dir" ]]; then
  echo "❌ Directory not found!"
  exit 1
fi

# Ask for force field
echo
echo "Select force field:"
echo "   1) UFF"
echo "   2) MMFF94"
echo "   3) MMFF94s"
read -p "Enter force field (UFF / MMFF94 / MMFF94s): " ff

# Validate force field
if [[ ! "$ff" =~ ^(UFF|MMFF94|MMFF94s)$ ]]; then
  echo "❌ Invalid force field. Choose from UFF, MMFF94, MMFF94s"
  exit 1
fi

# Create output directory
out_dir="rdkit_minimized"
mkdir -p "$out_dir"

# Run RDKit minimization in Python
python3 <<EOF
import os
from rdkit import Chem
from rdkit.Chem import AllChem

input_dir = "$mol_dir"
output_dir = "$out_dir"
ff = "$ff"

for fname in os.listdir(input_dir):
    if not fname.lower().endswith(('.sdf', '.pdb')):
        continue

    file_path = os.path.join(input_dir, fname)

    # Read molecule
    mol = None
    if fname.lower().endswith('.sdf'):
        suppl = Chem.SDMolSupplier(file_path, removeHs=False)
        mol = suppl[0] if suppl and suppl[0] is not None else None
    else:
        mol = Chem.MolFromPDBFile(file_path, removeHs=False)

    if mol is None:
        print(f"❌ Skipping {fname}: unable to read")
        continue
        
    # Prpare the molecules
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())

    try:
        if ff == "UFF":
            AllChem.UFFOptimizeMolecule(mol)
        else:
            if not AllChem.MMFFHasAllMoleculeParams(mol):
                print(f"⚠️ Skipping {fname}: missing MMFF parameters")
                continue
            props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="s" if ff == "MMFF94s" else "")
            ff_obj = AllChem.MMFFGetMoleculeForceField(mol, props)
            ff_obj.Minimize()
    except Exception as e:
        print(f"❌ Error minimizing {fname}: {e}")
        continue

    out_file = os.path.join(output_dir, os.path.splitext(fname)[0] + ".pdb")
    Chem.MolToPDBFile(mol, out_file)

    # Add author metadata as REMARKs
    try:
        with open(out_file, "r") as f:
            pdb_content = f.readlines()

        with open(out_file, "w") as f:
            f.write("REMARK  AUTHOR: Amr Al-Hefnawy - Nile University\n")
            f.write("REMARK  SOURCE: RDKit Minimized Structure\n")
            f.write("REMARK  LICENSE: For academic use only\n")
            f.writelines(pdb_content)
    except Exception as e:
        print(f"⚠️ Failed to insert metadata into {fname}: {e}")

    print(f"✅ Saved minimized with metadata: {out_file}")
EOF

echo
echo "🎉 All molecules processed successfully."
echo "Output directory: $out_dir/"

