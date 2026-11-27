import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw

oli_df = pd.read_csv("polymers_dataset.csv", usecols=["smiles", "label"])
print(f"This dataframe has {oli_df.shape[0]} rows and {oli_df.shape[1]} columns.")
oli_df = oli_df.copy()
oli_df["MolObj"] = oli_df["smiles"].apply(Chem.MolFromSmiles)
oli_df = oli_df.copy()
oli_df = oli_df.dropna(subset=["MolObj"])
print(oli_df.shape)
oli_df = oli_df.copy()
oli_df["Canon_SMILES"] = oli_df["MolObj"].apply(Chem.MolToSmiles)
n_duplicates = oli_df["Canon_SMILES"].duplicated().sum()
oli_df = oli_df.drop("smiles", axis=1)
print(f"The number of duplicated canon SMILES is: {n_duplicates}")
oli_df = oli_df.drop_duplicates(subset=["Canon_SMILES"])
oli_df = oli_df.copy()
print(oli_df.shape)

#Computing some descriptors:
oli_df["MW"] = oli_df["MolObj"].apply(Descriptors.MolWt).round(2)
oli_df["LogP"] = oli_df["MolObj"].apply(Descriptors.MolLogP).round(2)
oli_df = oli_df.copy()
oli_df["HBD"] = oli_df["MolObj"].apply(Descriptors.NumHDonors)
oli_df["HBA"] = oli_df["MolObj"].apply(Descriptors.NumHAcceptors)
oli_df = oli_df.copy()
oli_df["TPSA"] = oli_df["MolObj"].apply(Descriptors.TPSA).round(2)
oli_df["CSP3_fraction"] = oli_df["MolObj"].apply(Descriptors.FractionCSP3).round(2)
oli_df = oli_df.copy()
oli_df["Num_Rot_Bonds"] = oli_df["MolObj"].apply(Descriptors.NumRotatableBonds)
oli_df["Num_Heavy_Atoms"] = oli_df["MolObj"].apply(lambda x: x.GetNumHeavyAtoms())
oli_df["Num_Atoms"] = oli_df["MolObj"].apply(lambda x: x.GetNumAtoms(onlyExplicit=False))
oli_df = oli_df.copy()

#Printing the new columns:
print(oli_df[["label", "MW", "LogP", "TPSA", "HBD", "HBA", "CSP3_fraction", "Num_Rot_Bonds", "Num_Heavy_Atoms", "Num_Atoms"]].head())
print(oli_df[["label", "MW", "LogP", "TPSA", "HBD", "HBA", "CSP3_fraction", "Num_Rot_Bonds", "Num_Heavy_Atoms", "Num_Atoms"]].tail())

#Summary statistics per label:
summary = oli_df.groupby("label")["MW"].agg(["mean", "median", "min", "max", "std"]).round(2)
print(f"Summary statistics for molecular weight per label:\n {summary}")

#Comparing TPSA and LogP across labels:
mean_tpsa = oli_df.groupby("label")["TPSA"].mean().round(2)
mean_logp = oli_df.groupby("label")["LogP"].mean().round(2)
print(f"The average TPSA per label:\n {mean_tpsa}")
print(f"The average LogP per label:\n {mean_logp}")

#Making a patterns dictionnary holding necessary smarts:
patterns = {
    "amide": Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]"),
    "hydroxyl": Chem.MolFromSmarts("[OX2H]"),
    "benzene": Chem.MolFromSmarts("c1ccccc1"),
    "halogene": Chem.MolFromSmarts("[F,Cl,Br,I]"),
    "carboxylic_acid": Chem.MolFromSmarts("[CX3](=O)[OX2H1]"),
    }

for name, patt in patterns.items():
    oli_df[name] = oli_df["MolObj"].apply(lambda mol: mol.HasSubstructMatch(patt))
    oli_df = oli_df.copy()

print(oli_df[list(patterns.keys())].head(10))

#Computing the mean of oligomers having each pattern per label:
for patt in patterns.keys():
    print(f"The fraction of molecules having the {patt} group is:\n{oli_df.groupby('label')[patt].mean().round(2)}")

#Filtering the df into a special df for each label and sorting by MW:
oligosaccharide_df = oli_df[oli_df["label"]=="oligosaccharide"].sort_values(by="MW")
print(f"The oligosaccharide df has {oligosaccharide_df.shape[0]} molecules.")
peptide_df = oli_df[oli_df["label"]=="peptide"].sort_values(by="MW")
print(f"The peptide df has {peptide_df.shape[0]} molecules.")
plastic_df = oli_df[oli_df["label"]=="plastic"].sort_values(by="MW")
print(f"The plastic df has {plastic_df.shape[0]} molecules.")

#Visualizing the first and last molecule of each label and saving it in one .png file:
label_dfs = {
    "oligosaccharide": oligosaccharide_df,
    "peptide": peptide_df,
    "plastic": plastic_df,
}
mols_list = []
for df in label_dfs.values():
    mols_list.append(df["MolObj"].iloc[0])
    mols_list.append(df["MolObj"].iloc[-1])
oligo_img = Draw.MolsToImage(mols_list)
oligo_img.save("Oligomers.png")

#Droping the mol objects columns and saving each label's df into a seperate .csv file:
for name, df in label_dfs.items():
    df.drop("MolObj", axis=1).to_csv(f"{name}.csv", index=False)