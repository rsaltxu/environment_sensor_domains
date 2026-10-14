#!/usr/bin/env python3

"""
Python script to perform ANCOM-BC and estimate response to aridity
through W statistic with P-value Bonferroni correction
"""

from skbio.stats.composition import ancombc, multi_replace
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


###Load relative abundance dataframe
data=pd.read_csv('/data/singlem/australia/australia_custom_abundances.csv',header=0,index_col=0)


samples=data.index.to_list()
taxa=data.columns.to_list()
#ANCOM can't deal with 0s, therefore we replace them with small non-zeros
data=multi_replace(data)

data = pd.DataFrame(
    data,
    index=samples, columns=taxa)

#Load metadata
metadata=pd.read_csv("/home/roman.sanz/scripts/JAE_HK_Aridity/data/Australia_merged_metadata_wphylum.csv",header=0,index_col=0)

#Drop NA values and remove samples in the semi-arid gap (AI 0.2-0.65)
metadata = metadata.dropna(subset=["Aridity.index...annual.mean"])
metadata = metadata[
    (metadata["Aridity.index...annual.mean"] > 0.65) |
    (metadata["Aridity.index...annual.mean"] < 0.2)
]

#Make both dataframes be on the same order
data = data.loc[metadata.index]
metadata["Aridity_categorical"] = np.where(
    metadata["Aridity.index...annual.mean"] < 0.65,
    "_Arid",
    "Non-arid"
)
metadata=metadata[["Aridity_categorical"]]

print(ancombc(data, metadata, formula='Aridity_categorical',p_adjust='bonf'))
