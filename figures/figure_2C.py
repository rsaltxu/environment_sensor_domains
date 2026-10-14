#!/usr/bin/env python3

import pandas as pd
import gzip
#import polars as pl
import numpy as np



####LOAD FAMILY DICTIONARY
genome_id=[]
with open ('/home/roman.sanz/data/diamond_pblast/gtdb_reps/bac120_files_list.txt') as bac_list:
    for line in bac_list:
        genome_id.append(line[:-1])
print('genome list loaded')

family_dict={}
phylum_dict={}

with gzip.open('/support/databases/gtdb/genomic_files_reps_v226/aux/bac120_metadata.tsv.gz','rt',encoding='utf-8') as metadata_file:
    metadata = pd.read_table(metadata_file,sep='\t',index_col=0)
    metadata = metadata.loc[genome_id]
    print(genome_id[:5])
    print(metadata.index.to_list()[:5])
    taxonomies = []
    raw_tax = metadata['gtdb_taxonomy'].to_list()
    #for tax in raw_tax:
    for i in range(len(raw_tax)):  
        family=raw_tax[i].split(';')[4][3:]
        if family not in family_dict:
            family_dict[family]=[genome_id[i]]
        else:
            family_dict[family].append(genome_id[i])

    ###check phyla
    main_phyla=["Gammaproteobacteria","Actinomycetota",
    "Bacillota","Alphaproteobacteria","Bacteroidota",
    "Planctomycetota","Chloroflexota","Desulfobacterota",
    "Acidobacteriota","Cyanobacteriota","Verrucomicrobiota"]
    
    for i in range(len(raw_tax)):
        for phylum in main_phyla:
            if phylum in raw_tax[i]:
                if phylum not in phylum_dict:
                    phylum_dict[phylum]=[genome_id[i]]
                else:
                    phylum_dict[phylum].append(genome_id[i])                

def count_total_proteins(prodigal_annotation):
    count=0
    with gzip.open(prodigal_annotation,'rt',encoding='utf-8') as fin:
        for line in fin:
            if line.startswith('>'):
                count+=1
    return count


###CORRELATION HK/MB- Genome size
import numpy as np
import statsmodels.api as sm
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import font_manager
plt.rcParams["font.family"] = "Liberation Sans"

gtdb_metadata=pd.read_csv('/home/roman.sanz/data/diamond_pblast/gtdb_reps/hk_variables.csv',header=0,index_col=0)
genome_id=gtdb_metadata['Genome'].to_list()
print(gtdb_metadata[(gtdb_metadata["HK/Mb"] > 5)&(gtdb_metadata["Genome_size"] < 1)])
###parasytes to highlight
highlight = {
    "Bacillota": [],
    "Alphaproteobacteria": [],
    "Gammaproteobacteria": [],
    "Cyanobacteriota":[]
}

with gzip.open("/support/databases/gtdb/genomic_files_reps_v226/aux/bac120_metadata.tsv.gz", 'rt', encoding='utf-8') as annotation:
    metadata=pd.read_table(annotation,sep='\t',index_col=0,header=0)
    metadata=metadata.loc[genome_id]
    taxonomies=[]
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
        if 'seudomonadota' in tax:
            phyl=tax.split(';')[2][3:]
        else:
            phyl=tax.split(';')[1][3:]
            if '_' in phyl:
                phyl=phyl[:-2]
        taxonomies.append(phyl)
print(set(taxonomies))
gtdb_metadata['Phylum_v2']=taxonomies
#gtdb_metadata.to_csv('/home/roman.sanz/data/diamond_pblast/gtdb_reps/hk_variables.csv')
print(min(gtdb_metadata['Genome_size']))

x = gtdb_metadata["Genome_size"].to_numpy()
y = gtdb_metadata["HK/Mb"].to_numpy()
mask = np.isfinite(x) & np.isfinite(y)

x = x[mask]
y = y[mask]
""" x=x[y!=1]
y=y[y!=1] """
""" hk_coph_values=df["HK_coph"].to_numpy()
hk_coph_values=hk_coph_values[np.isfinite(hk_coph_values)]
hk_eucl_values=df["HK_eucl"].to_numpy()
hk_eucl_values=hk_eucl_values[np.isfinite(hk_eucl_values)] """
print(np.min(x), np.max(x))
print(np.min(y), np.max(y))
# Add intercept
X = sm.add_constant(x)

# Fit model
model = sm.OLS(y, X).fit()
print(model.params)
# Extract values
print(model.summary())
r2 = model.rsquared
p_value = model.pvalues[1]   # p-value for slope

print("R²:", r2)
print("p-value:", p_value) 

plt.figure(figsize=(14, 7))

ax=sns.scatterplot(
    data=gtdb_metadata,
    x="Genome_size",
    y="HK/Mb",
    s=6,
    palette={
        "Acidobacteriota" : "#686F98",
        "Actinomycetota" : "#E36414",
        "Bacteroidota" : "#8E9A87",
        "Chloroflexota" : "#577590",
        "Cyanobacteriota" : "#81B984",
        "Desulfobacterota" : "#50878E",
        "Bacillota" : "#EC8517",
        "Gemmatimonadota" : "#6289BD",
        "Myxococcota" : "#6A4C93",
        "Alphaproteobacteria" : "#CC8666",
        "Gammaproteobacteria" : "#3793AC",
        "Planctomycetota" : "#FCCA46",
        "Verrucomicrobiota" : "#43AA8B",
        'Campylobacterota':"grey",
          'Spirochaetota':"grey",
            'Bacteroidota_A':"#8E9A87", 
            'Pseudomonadota':"grey"},
    hue='Phylum'
)
sns.regplot(x=x, y=y, fit_reg=True,scatter=False)
#plt.title(f"GTDB vs {metric} distances")
plt.xlabel("Genome size (Mb)",fontsize=18)
plt.ylabel("HK/Mb",fontsize=18)
plt.xlim(right=1.5,)
plt.ylim(top=6)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()
plt.savefig(f'/home/roman.sanz/data/figures/HKdens_genomesize_regr_zoom.png',dpi=300)

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

plt.savefig(
    "/home/roman.sanz/data/figures/HKdens_genomesize_regr_2x3_v2.png",
    dpi=300
)
