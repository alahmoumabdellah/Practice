from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit import DataStructs
 
smiles_list = [
    "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
    "c1ccccc1",                  # Benzene
    "CC(=O)O",                   # Acetic acid
    "CCO"                        # Ethanol
]
 
molecules = [Chem.MolFromSmiles(s) for s in smiles_list]
filterpass = [j for j in molecules if Descriptors.MolWt(j)<200]
fgplist = [Chem.RDKFingerprint(mol) for mol in filterpass]
benzenefgp = Chem.RDKFingerprint(molecules[1])
similarities = DataStructs.BulkTanimotoSimilarity(benzenefgp, fgplist)
print(similarities)