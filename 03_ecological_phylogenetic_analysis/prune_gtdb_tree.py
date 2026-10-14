#!/usr/bin/env python3

"""
Python script to prune the bacterial tree of life from GTDB
"""

from ete4 import Tree
import gzip

###Load GTDB tree as Tree object

with gzip.open("/support/databases/gtdb/genomic_files_reps_v226/aux/bac120.tree.gz", "rt",encoding='utf-8') as f:
    t = Tree(f.read(), parser=1)#GTDB format has internal node labels, needs format=1 to be read

#load genome IDs to keep
with open("/home/roman.sanz/data/diamond_pblast/gtdb_reps/family_reps_top_phyla_list.txt") as f:
    keep = set(line.strip() for line in f)

# prune

t.prune(keep, preserve_branch_length=True)

# write pruned tree
t.write(outfile="/data/phylomodel/family_reps_top_phyla_gtdb.tree",parser=0)