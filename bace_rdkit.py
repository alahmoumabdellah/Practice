import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
import warnings
#warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

bacedf = pd.read_csv("bace_df_revised.csv")

bacedf["ROMol"] = bacedf["mol"].apply(Chem.MolFromSmiles)
print(bacedf["ROMol"].head())
#droping rows with invalid smiles 
bacedf = bacedf.dropna(subset=["ROMol"])
print(f"The number of molecules before canonicalization is: {bacedf.shape[0]}")
bacedf = bacedf.copy()
#Canonicalization and duplicates check:
bacedf["Canon_SMILES"] = bacedf["ROMol"].apply(Chem.MolToSmiles)
n_duplicates = bacedf["Canon_SMILES"].duplicated().sum()
print(f"The number of molecules with same canonical SMILES is: {n_duplicates}")
#Dropping duplicates with the same canon smiles:
bacedf = bacedf.drop_duplicates(subset=["Canon_SMILES"])


print(f"The final number of molecules is: {bacedf.shape[0]}")
#drawing the first 8 molecules:
mols = bacedf["ROMol"].iloc[:8].to_list()
img1 = Draw.MolsToImage(mols)
img1.show()

bacedf = bacedf.copy()

"""Computing the values of some descriptors for the molecules
and comparing their values with the data's original values"""
bacedf["Molecular_weight"] = bacedf["ROMol"].apply(Descriptors.MolWt)
bacedf["LogP"] = bacedf["ROMol"].apply(Descriptors.MolLogP)
print(bacedf[["MW", "Molecular_weight", "AlogP", "LogP"]].head())
bacedf["HBDonors"] = bacedf["ROMol"].apply(Descriptors.NumHDonors)
bacedf["HBAcceptors"] = bacedf["ROMol"].apply(Descriptors.NumHAcceptors)
print(bacedf[["HBD", "HBDonors", "HBA", "HBAcceptors"]].head())

#Seeing how many molecules pass the Lipinski rule of 5:
bacedf["Lipinski"] = (
    (bacedf["MW"]<=500)& 
    (bacedf["AlogP"]<=5)& 
    (bacedf["HBD"]<=5)&
    (bacedf["HBA"]<=10)
)
print(f"Total passing Lipinski filter: {bacedf["Lipinski"].sum()}")

#seeing which condition fails most:
failmw = (bacedf["MW"] > 500).sum()
failalogp = (bacedf["AlogP"] > 5).sum()
failhbd = (bacedf["HBD"] > 5).sum()
failhba = (bacedf["HBA"] > 10).sum()
print(f"Failed Molecular weight (>500): {failmw}")
print(f"Failed LogP (>5): {failalogp}")
print(f"Failed HBDonnors (>5): {failhbd}")
print(f"Failed HBAcceptors (>10): {failhba}")

pass_rates = bacedf.groupby("Class")["Lipinski"].mean()*100
print(f"The pass rate of each class is:\n{pass_rates}")

#Converting the molecular weight column into a numpy array 
mwarr = bacedf["Molecular_weight"].to_numpy()
mwmean = mwarr.mean()
mwstd = np.std(mwarr, ddof=1)
print(f"The average molecular weight is:{mwmean:.2f}",  f"While the standard deviation is:{mwstd:.2f}")

#Computing the number of rotatable bonds in order to use it for veber's rules filtering:
bacedf["Rot_Bonds_Num"] = bacedf["ROMol"].apply(Descriptors.NumRotatableBonds)

#Seeing how many molecule pass the veber rules using the data original value of PSA:
bacedf["Veber_org_pass"] = ((bacedf["Rot_Bonds_Num"]<= 10) & (bacedf["PSA"]<= 140))

#Seeing how many molecule pass the veber rules with the computed value of PSA:
bacedf["TPSA"] = bacedf["ROMol"].apply(Descriptors.TPSA)
bacedf["Veber_rdkit_pass"] = ((bacedf["Rot_Bonds_Num"]<= 10) & (bacedf["TPSA"]<= 140))

#Comparing the number of molecules that passed each filter:
print(f"Number of passing mols using original PSA :{bacedf["Veber_org_pass"].sum()}", f"While the number of mols that passed using the computed PSA is:{bacedf["Veber_rdkit_pass"].sum()}")

veber_pass_rates = bacedf.groupby("Class")["Veber_org_pass"].mean()*100
print(f"The pass rate for each class is:\n {veber_pass_rates}")