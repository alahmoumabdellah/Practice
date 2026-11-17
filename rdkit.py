import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

bacedf = pd.read_csv("bace_processed.csv")

bacedf["ROMol"] = bacedf["mol"].apply(Chem.MolFromSmiles)
print(bacedf["ROMol"].head())
print(bacedf.shape[0])
#droping rows with invalid smiles 
bacedf = bacedf.dropna(subset=["ROMol"])
print(bacedf.shape[0])

bacedf["Molecular_weight"] = bacedf["ROMol"].apply(Descriptors.MolWt)
bacedf["LogP"] = bacedf["ROMol"].apply(Descriptors.MolLogP)
print(bacedf[["ROMol", "Molecular_weight", "LogP"]].head())

#Converting the molecular weight column into a numpy array 
mwarr = bacedf["Molecular_weight"].to_numpy()
mwmean = mwarr.mean()
mwstd = np.std(mwarr, ddof=1)
print(mwmean, mwstd)

filtered_bacedf = bacedf.loc[(bacedf["LogP"]<5)&(bacedf["Molecular_weight"]< 500)]
print(filtered_bacedf.shape[0])