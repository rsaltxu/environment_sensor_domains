#!/usr/bin/env python3
"""
Python script to process output from TMHMM and extract a list of extracellular protein IDs
"""
#Import statements
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

"""
Previously to this, we ran a BLAST search to get the whole protein sequence where the sensor domains
were found, on which we ran TMHMM
"""

tmhmm={}
loc_dict={'outside':'E',
          'TMhelix':'H',
          'inside':'I'}

#For each protein, we build a string indicating their location, e.g. EEEEEEHHHHHHIIIIIIIIII
with open('/data/tmhmm/HK_TMHMM_results.txt') as tmhmm_results:
    for line in tmhmm_results:
        if not line.startswith('#') and 'slurmstepd' not in line:
            if line.strip()=='':
                break
            #print(line)
            line=line.split()
            prot_id=line[0]
            loc=line[2]
            start=int(line[3])
            end=int(line[4])
            if prot_id not in tmhmm:
                tmhmm[prot_id]=loc_dict[loc]*(end-start+1)
            else:
                tmhmm[prot_id]=tmhmm[prot_id]+loc_dict[loc]*(end-start+1)


match_df=pd.read_table('/data/HK_related/HK_db_BLAST/HK_db_BLASTp_unified.txt',sep='\t',header=None)
'qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore'#BLAST output format

hk_list=list(set(match_df[0].to_list()))


in_count=0
memb_count=0
out_count=0

###Now we will build a list with the sensor domains classified as extracellular
with open('/data/tmhmm/outer_HKs.txt','w') as outer_hks:
    for hk in hk_list:
        protein_match=match_df[match_df[0]==hk][1].values[0]
        ###for every protein, look for the match and start/end position
        start_pos=match_df[match_df[0]==hk][8].values[0]
        end_pos=match_df[match_df[0]==hk][9].values[0]
        
        ###protein sequence stats:
        seq=tmhmm[protein_match][start_pos:end_pos+1] #we retrieve the location string for the region comprehending the sensor domain only
        #we quantify the proportion for each cell location
        tmhmm_in=seq.count('I')/len(seq)
        tmhmm_memb=seq.count('H')/len(seq)
        tmhmm_out=seq.count('E')/len(seq)
        if max([tmhmm_in,tmhmm_memb,tmhmm_out])==tmhmm_in:
            in_count+=1
        if max([tmhmm_in,tmhmm_memb,tmhmm_out])==tmhmm_memb:
            memb_count +=1
        if max([tmhmm_in,tmhmm_memb,tmhmm_out])==tmhmm_out: #if the majority is predicted as E (extracellular) we add it in the list
            out_count += 1
            outer_hks.write(f'{hk}\n')

print(f'Inner: {in_count}\nMembrane: {memb_count}\nOuter: {out_count}')#for general overview purposes
