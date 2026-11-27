from rdkit import Chem
from rdkit.Chem import Descriptors, Draw

mol_supp = Chem.SDMolSupplier("benzodiazepine.sdf")
benzodiazepine = mol_supp[0]

molecular_weight = Descriptors.MolWt(benzodiazepine)
print(f"The molecular weight: {molecular_weight:.2f}")
log_p = Descriptors.MolLogP(benzodiazepine)
print(f"The lipophilicity factor: {log_p:.2f}")
tpsa = Descriptors.TPSA(benzodiazepine)
print(f"The TPSA: {tpsa}")
num_HBD = Descriptors.NumHDonors(benzodiazepine)
print(f"The number of hydrogen bond donors: {num_HBD}")
num_HBA = Descriptors.NumHAcceptors(benzodiazepine)
print(f"The number of hydrogen bond acceptors: {num_HBA}")
num_rot_bonds = Descriptors.NumRotatableBonds(benzodiazepine)
print(f"The number of rotatable bonds: {num_rot_bonds}")
num_atoms = benzodiazepine.GetNumAtoms()
print(f"The total number of atoms: {num_atoms}")

img1= Draw.MolToImage(benzodiazepine)
img1.show()
img1.save("benzodiazepine.png")
print("Viscaaaaaaaaaaaaaaa barcaaaaaaaaaaaaa\nY visca catalunya")

#print(list(benzodiazepine.GetPropNames()))
print(f"IUPAC name: {benzodiazepine.GetProp("PUBCHEM_IUPAC_NAME")}")
print(f"IUPAC systematic name: {benzodiazepine.GetProp("PUBCHEM_IUPAC_SYSTEMATIC_NAME")}")
print(f"Molecular formula: {benzodiazepine.GetProp("PUBCHEM_MOLECULAR_FORMULA")}")
print(f"MW: {benzodiazepine.GetProp("PUBCHEM_MOLECULAR_WEIGHT")}")
print(f"HBA: {benzodiazepine.GetProp("PUBCHEM_CACTVS_HBOND_ACCEPTOR")}")
print(f"HBD: {benzodiazepine.GetProp("PUBCHEM_CACTVS_HBOND_DONOR")}")
print(f"TPSA: {benzodiazepine.GetProp("PUBCHEM_CACTVS_TPSA")}")
print(f"Def stereo count: {benzodiazepine.GetProp("PUBCHEM_ATOM_DEF_STEREO_COUNT")}")
print(f"UDef stereo count: {benzodiazepine.GetProp("PUBCHEM_ATOM_UDEF_STEREO_COUNT")}")
