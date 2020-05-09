#!/bin/bash -ue
mafft --anysymbol --retree 0 --treeout --parttree --reorder seatoxin.fa

t_coffee -other_pg seq_reformat -in seatoxin.fa.tree -in2 seatoxin.fa -input newick -action +mafftnewick2newick > seatoxin.parttreednd0.dnd
