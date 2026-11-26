import pandas as pd
import numpy as np

bacedf= pd.read_csv("bace_processed.csv")
print(bacedf.head())
print(bacedf.tail())
print(bacedf.shape)
print(bacedf.dtypes.value_counts())
print(bacedf.columns[:20])
print(bacedf.describe())


print(bacedf["MW"])
print(bacedf[["MW", "AlogP", "pIC50"]])
print(bacedf.iloc[:10])
# I couldnt use just loc or just iloc because I only knew the indexes of rows from 5 to 15 and the names of columns
print(bacedf[["MW", "pIC50"]].iloc[5:16, ])
#renaming a column and setting the inplace False just to prove the capability of renaming:
bacedf = bacedf.rename(columns= {"ChiralCenterCount": "CCC"}, inplace=False)

nan_count = bacedf.isna().sum()
print(nan_count)
print(nan_count.sort_values(ascending=False).head(10))
bacedf.loc[bacedf.index[:10:2], "dCH2_Key"]= np.nan
print(bacedf["dCH2_Key"].iloc[:10])
bacedf.loc[bacedf.index[:10:2], "dCH2_Key"] = bacedf["dCH2_Key"].mean()
print(bacedf["dCH2_Key"].iloc[:10])
copydf = bacedf.copy()
copydf.loc[bacedf.index[:10:2], "dCH2_Key"]= np.nan
copydf = copydf.dropna(how= "any")
print(bacedf.shape)
print(copydf.shape)
print(bacedf.duplicated().sum())
bacedf.drop_duplicates(inplace=True)

copydf = copydf.loc[copydf['pIC50']> copydf["pIC50"].median()]
rowbefore = bacedf.shape[0]
copydf = copydf.loc[(copydf['Class']==1)&(copydf["MW"]> 350)]
# using | if we want one of the conditions to be met
rowafter = copydf.shape[0]
print(rowafter)
#task 21 :
print(f"fraction is {rowafter / rowbefore}")
copydf = copydf.loc[(copydf['HBA']<=1)&(copydf["HBA"]>= 5)]

print(bacedf.sort_values(by="MW", ascending=False).head(10))
bacedf["random_column"] = bacedf["MW"]/bacedf["pIC50"]
print(bacedf["random_column"].head())
conditions = [bacedf["pIC50"]>=9, (bacedf["pIC50"]<9) &(bacedf["pIC50"]>=7), bacedf["pIC50"]<7]
choices = ["High", "Medium", "Low"]
bacedf["Category_column"] = np.select(conditions, choices, default="Unkown")
print(bacedf["Category_column"].head())
print(bacedf["Category_column"].value_counts())

mean_values = bacedf.groupby("Category_column")["MW"].mean()
print(mean_values)
summary = bacedf.groupby("Category_column")["pIC50"].agg(["count", "mean", "max", "min"])
print(summary)
meancc = bacedf.groupby(["Class", "Category_column"])['MW'].mean()
print(meancc)
#seeing which grp of category column has the highest average and printing it:
print(mean_values.sort_values(ascending=False).head(1))

arr1 = bacedf["MW"].to_numpy()
print(arr1.shape)
#comparing the mean, std, max, and min between np and pd
meanarr = arr1.mean()
print(meanarr)
meanclm = bacedf['MW'].mean()
print(meanclm)
stdarr = np.std(arr1, ddof=1)
print(stdarr)
stdclm = bacedf["MW"].std()
print(stdclm)
maxarr1 = np.max(arr1)
print(maxarr1)
maxclm = bacedf["MW"].max()
print(maxclm)
minarr1 = np.min(arr1)
print(minarr1)
minclm = bacedf["MW"].min()
print(minclm)

print((arr1-meanarr)/stdarr)
print(np.sum(arr1>(meanarr+stdarr)))
arr2 = bacedf["AlogP"].to_numpy()

corrnp = np.corrcoef(arr1, arr2)
print(corrnp)
corrpd = bacedf["MW"].corr(bacedf["AlogP"])
print(corrpd)

bacedf.to_csv('claudetask.csv')