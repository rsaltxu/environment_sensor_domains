#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import gzip
import numpy as np

plt.rcParams["font.family"] = "Liberation Sans"

color_dict ={
"Acidobacteriota" : "#686F98",
"Actinomycetota" : "#E36414",
"Armatimonadota" : "#4888AE",
"Bacteroidota" : "#8E9A87",
"Bdellovibrionota" : "#D53E4F",
"Chloroflexota" : "#577590",
"Crenarchaeota" : "#9D4EDD",
"Cyanobacteriota" : "#81B984",
"Dependentiae" : "#F4A728",
"Desulfobacterota" : "#50878E",
"Pseudomonadota" : "#956A7E",
"FCPU426" : "#1B4965",
"Fibrobacterota" : "#49988D",
"Bacillota" : "#EC8517",
"Gemmatimonadota" : "#6289BD",
"Common_pool" : "#AD162E",
"Methylomirabilota" : "#2EC4B6",
"Myxococcota" : "#6A4C93",
"Nitrospirota" : "#43AA8B",
"Alphaproteobacteria" : "#CC8666",
"Gammaproteobacteria" : "#3793AC",
"Patescibacteria" : "#B2231A",
"Planctomycetota" : "#FCCA46",
"Verrucomicrobiota" : "#43AA8B",
"x_Other" : "grey"}


def count_total_proteins(prodigal_annotation):
    count=0
    with gzip.open(prodigal_annotation,'rt',encoding='utf-8') as fin:
        for line in fin:
            if line.startswith('>'):
                count+=1
    return count

def top_11_frequent(strings):
    freq = {}
    for s in strings:
        freq[s] = freq.get(s, 0) + 1

    return [k for k, v in sorted(freq.items(), key=lambda x: x[1], reverse=True)[:11]]

genome_id=[]
with open ('/data/bac120_files_list.txt') as bac_list:
    for line in bac_list:
        genome_id.append(line[:-1])

with gzip.open("/support/databases/gtdb/genomic_files_reps_v226/aux/bac120_metadata.tsv.gz", 'rt', encoding='utf-8') as annotation:
    metadata=pd.read_table(annotation,sep='\t',index_col=0,header=0)
    metadata=metadata.loc[genome_id]
    taxonomies=[]
    raw_tax= metadata['gtdb_taxonomy'].to_list()
    for tax in raw_tax:
        #We divide Pseudomonadota into Alphaproteobacteria and Gammaproteobacteria
        if 'seudomonadota' in tax:
            phyl=tax.split(';')[2][3:]
        else:
            phyl=tax.split(';')[1][3:]
            if '_' in phyl:
                phyl=phyl[:-2]
        taxonomies.append(phyl)

top_11=top_11_frequent(taxonomies)#Only analyze 11 most frequent taxa

HK_counts=[]
phyla=[]
print('reading matrix...')
gene_count_matrix = pd.read_csv("/home/roman.sanz/data/diamond_pblast/gtdb_reps/HK_abundance_matrix_prot_90cov_50score.csv", index_col=0)
for id in genome_id:
        print(f'{genome_id.index(id)}/{len(genome_id)}')
        if taxonomies[genome_id.index(id)] in top_11:
            prodigal_annotation=f'/support/databases/gtdb/genomic_files_reps_v226/protein_faa_reps/bacteria/{id}_protein.faa.gz'

            total_prot=count_total_proteins(prodigal_annotation=prodigal_annotation)
            
            HK_count = gene_count_matrix.loc[id].sum()*total_prot
            
            Mb_value = metadata.loc[id]['genome_size']/1000000
            
            HK_counts.append(HK_count)
            
            phyla.append(taxonomies[genome_id.index(id)])
df = pd.DataFrame({
    "#HK": HK_counts,
    "Phylum": phyla
}) 

#Retrieve genome size
genome_size=[]
for id in genome_id:
        print(f'{genome_id.index(id)}/{len(genome_id)}')
        if taxonomies[genome_id.index(id)] in top_11:
            Mb_value = metadata.loc[id]['genome_size']/1000000
            genome_size.append(Mb_value)
df['genome_size']=genome_size

#Calculate sensor domain density
hk_norm=[]
for i in range(len(df['genome_size'].to_list())):
    hk_norm.append(df.iloc[i]['#HK']/df.iloc[i]['genome_size'])
df['HK_per_Mb']=hk_norm
phyla=df['Phylum'].to_list()
list_phylum=list(set(phyla))
means={}

for phylum in list_phylum:
    means[phylum]=np.median(df[df['Phylum']==phylum]['HK_per_Mb'].to_list())
means=sorted(means, key=means.get, reverse=True)
order=[phyl for phyl in means]

plt.figure(figsize=(12, 8))

sns.violinplot(
    data=df,
    x="Phylum",
    y="HK_per_Mb",
    order=order,
    palette=color_dict

)

plt.ylabel("Sensor domains per Mb",fontsize=24)
plt.ylim(top=50,bottom=-8)
plt.yticks(fontsize=22)
plt.xticks(rotation=45,fontsize=24,ha='right',rotation_mode='anchor')
plt.tight_layout()
plt.show()