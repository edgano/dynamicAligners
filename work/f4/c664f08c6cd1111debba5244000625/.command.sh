#!/bin/bash -ue
t_coffee -other_pg aln_compare              -al1 seatoxin.ref              -al2 seatoxin.reg_align.1000.dynamic.10000.dynamic.with.parttreednd0.tree.aln             -compare_mode sp             | grep -v "seq1" | grep -v '*' |             awk '{ print $4}' ORS="	"             > "seatoxin.reg_align.1000.dynamic.with.parttreednd0.tree.sp"

t_coffee -other_pg aln_compare              -al1 seatoxin.ref              -al2 seatoxin.reg_align.1000.dynamic.10000.dynamic.with.parttreednd0.tree.aln             -compare_mode tc             | grep -v "seq1" | grep -v '*' |             awk '{ print $4}' ORS="	"             > "seatoxin.reg_align.1000.dynamic.with.parttreednd0.tree.tc"

t_coffee -other_pg aln_compare              -al1 seatoxin.ref              -al2 seatoxin.reg_align.1000.dynamic.10000.dynamic.with.parttreednd0.tree.aln             -compare_mode column             | grep -v "seq1" | grep -v '*' |               awk '{ print $4}' ORS="	"             > "seatoxin.reg_align.1000.dynamic.with.parttreednd0.tree.col"
