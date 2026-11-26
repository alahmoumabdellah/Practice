import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw

mols = Chem.SDMolSupplier("pubchem_3mols.sdf")
print(f"The supplier has unreaden mols: {None in mols}")
#Droping the None ones just in case:
mols_list = [m for m in mols if m is not None]
#mols_list = list(mols)
print(f"The number of molecules read is: {len(mols_list)}")

print(f"The properties names are: {list(mols[0].GetPropNames())}")
print("The values of the first 5 properties:\n")
for prop in list(mols[0].GetPropNames())[:5]:
    print(prop,mols[0].GetProp(prop), sep=" : ")

#Inspecting one molecule:
num_atoms = mols[1].GetNumAtoms()
num_atoms_Hs = mols[1].GetNumAtoms(onlyExplicit=False)
num_bonds = mols[1].GetNumBonds()
num_rings = mols[1].GetRingInfo().NumRings()
print(
    f"The total number of atoms: {num_atoms}\n",
    f"The total number of atoms including Hydrogen atoms: {num_atoms_Hs}\n",
    f"The total number of bonds: {num_bonds}\n",
    f"The total number of rings: {num_rings}"
)

for mol in mols_list:
    print({mol.GetProp('PUBCHEM_IUPAC_OPENEYE_NAME')})
    for atom in mol.GetAtoms():
        nbr_labels = [
            f"{nbr.GetSymbol()}({nbr.GetIdx()})" for nbr in atom.GetNeighbors()
        ]
        print(
            f"Atome: {atom.GetSymbol()}({atom.GetIdx()}), neighbors:\n {", ".join(nbr_labels)}"
        )

#Drawing the molecules:
legends = [mol.GetProp("PUBCHEM_IUPAC_OPENEYE_NAME") for mol in mols]
img1 = Draw.MolsToImage(mols_list, legends=legends)
img1.save("Molecules.png")