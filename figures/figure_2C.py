#!/usr/bin/env python3

import pandas as pd
import gzip
import numpy as np
import statsmodels.api as sm
import seaborn as sns
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Liberation Sans"

# Phyla to highlight
highlight_phyla = [
    "Actinomycetota",
    "Bacillota",
    "Desulfobacterota",
    "Alphaproteobacteria",
    "Gammaproteobacteria",
    "Cyanobacteriota"
]

# Colors
colors = {
    "Actinomycetota": "#E36414",
    "Bacillota": "#EC8517",
    "Desulfobacterota": "#50878E",
    "Cyanobacteriota": "#81B984",
        "Alphaproteobacteria" : "#CC8666",
        "Gammaproteobacteria" : "#3793AC"
}

parasite_palette={
    "Actinomycetota": "#E36414",
    "Bacillota": "#D63D00EB",
    "Desulfobacterota": "#50878E",
    "Cyanobacteriota": "#009908",
        "Alphaproteobacteria" : "#C56900",
        "Gammaproteobacteria" : "#8F00D1"
}

###CORRELATION HK/MB- Genome size

gtdb_metadata=pd.read_csv('/home/roman.sanz/data/diamond_pblast/gtdb_reps/hk_variables.csv',header=0,index_col=0)
genome_id=gtdb_metadata['Genome'].to_list()

###parasytegenomes to be highlighted
highlight = {
    "Bacillota": [],
    "Alphaproteobacteria": [],
    "Gammaproteobacteria": [],
    "Cyanobacteriota":[]
}

with gzip.open("/support/databases/gtdb/genomic_files_reps_v226/aux/bac120_metadata.tsv.gz", 'rt', encoding='utf-8') as annotation:
    metadata=pd.read_table(annotation,sep='\t',index_col=0,header=0)
    metadata=metadata.loc[genome_id]
    raw_tax= metadata['gtdb_taxonomy'].to_list()
    for i in range(len(raw_tax)):
        if 'g__Rickettsia' in raw_tax[i]:
            highlight['Alphaproteobacteria'].append(genome_id[i])
        elif 'g__Buchnera' in raw_tax[i]:
            highlight['Gammaproteobacteria'].append(genome_id[i])
        elif 'g__Mycoplasma' in raw_tax[i]:
            highlight['Bacillota'].append(genome_id[i])
        elif 'g__Prochlorococcus' in raw_tax[i]:
                highlight['Cyanobacteriota'].append(genome_id[i])

with gzip.open("/support/databases/gtdb/genomic_files_reps_v226/aux/bac120_metadata.tsv.gz", 'rt', encoding='utf-8') as annotation:
    metadata=pd.read_table(annotation,sep='\t',index_col=0,header=0)
    metadata=metadata.loc[genome_id]
    taxonomies=[]
    raw_tax= metadata['gtdb_taxonomy'].to_list()
    for tax in raw_tax:
        if 'seudomonadota' in tax: #Split Pseudomonadota in Alphaproteobacteria and Gammaproteobacteria
            phyl=tax.split(';')[2][3:]
        else:
            phyl=tax.split(';')[1][3:]
            if '_' in phyl:
                phyl=phyl[:-2]
        taxonomies.append(phyl)
        
gtdb_metadata['Phylum_v2']=taxonomies

x = gtdb_metadata["Genome_size"].to_numpy()
y = gtdb_metadata["HK/Mb"].to_numpy()
mask = np.isfinite(x) & np.isfinite(y)

x = x[mask]
y = y[mask]

# Add intercept
X = sm.add_constant(x)

# Fit model
model = sm.OLS(y, X).fit()

# Extract values
r2 = model.rsquared
p_value = model.pvalues[1]   # p-value for slope

###Make the plot

fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True, sharey=True)

for ax, phylum in zip(axes.flat, highlight_phyla):

    # Background (other phyla)
    others = gtdb_metadata[gtdb_metadata["Phylum_v2"] != phylum]
    sns.scatterplot(
        data=others,
        x="Genome_size",
        y="HK/Mb",
        color="lightgrey",
        s=5,
        alpha=0.25,
        linewidth=0,
        ax=ax,
        legend=False
    )
    sns.regplot(x=x, y=y, fit_reg=True,scatter=False,)

    # Highlight phylum
    subset = gtdb_metadata[gtdb_metadata["Phylum_v2"] == phylum]
    sns.scatterplot(
        data=subset,
        x="Genome_size",
        y="HK/Mb",
        color=colors[phylum],
        s=12,
        ax=ax,
        legend=False
    )

    if phylum in highlight:
        print(phylum)
        highlighted = subset[subset["Genome"].isin(highlight[phylum])]
        print(highlighted.shape)
        sns.scatterplot(
            data=highlighted,
            x="Genome_size",
            y="HK/Mb",
            color=parasite_palette[phylum],
            s=15,              # larger points
            edgecolor="black", # optional
            linewidth=0.7,
            zorder=10,         # draw on top
            legend="full",
            ax=ax
        )

    sns.regplot(x=x, y=y, fit_reg=True,scatter=False,ax=ax,color="#333333")
    ax.set_title(phylum, fontsize=20)
    ax.tick_params(labelsize=16)
    ax.set_ylabel('')
    ax.set_xlabel('')

plt.tight_layout()

plt.show()