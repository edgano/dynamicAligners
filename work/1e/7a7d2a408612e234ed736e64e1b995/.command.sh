#!/bin/bash -ue
## echo 'psicoffee_msa 20' > config.txt
## echo 'clustalo_msa 5000' >> config.txt
## echo 'famsa_msa 100000000' >> config.txt
## -dynamic_config /Users/edgargarriga/CBCRG/dynamicAligners/config.txt   
export blast_server_4_CLTCOFFEE=LOCAL
export protein_db_4_CLTCOFFEE=/users/cn/egarriga/datasets/db/uniref50.fasta


t_coffee -reg -reg_method dynamic_msa          -seq seatoxin.fa          -reg_nseq 1000          -dynamic 10000          -reg_homoplasy          -outfile seatoxin.reg_align.1000.dynamic.10000.dynamic.with.codnd.tree.aln
