#!/usr/bin/env nextflow

/*
 * Copyright (c) 2017-2018, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'XXXXXX'.
 *
 *   XXXXXX is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   XXXXXX is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with XXXXXX.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * Main XXX pipeline script
 *
 * @authors
 * Edgar Garriga

 */

/*
 * defaults parameter definitions
 */

seq2improve="cryst,blmb,rrm,subt,ghf5,sdr,tRNA-synt_2b,zf-CCHH,egf,Acetyltransf,ghf13,p450,Rhodanese,aat,az,cytb,proteasome,GEL"

//params.seqs ="${baseDir}/test/*.fa"
params.seqs ="/users/cn/egarriga/datasets/homfam/combinedSeqs/{$seq2improve}.fa"

// input reference sequences aligned in
//params.refs ="${baseDir}/test/*.ref"
params.refs ="/users/cn/egarriga/datasets/homfam/refs/*.ref"

// input guide trees in Newick format. Or `false` to generate trees
//params.trees = "/users/cn/egarriga/datasets/homfam/trees/*.{FAMSA,CLUSTALO,MAFFT_PARTTREE}.dnd"
params.trees =""

//aligner and tree generation
params.align_method = "dynamic"
params.tree_method = "codnd,parttreednd0,famsaSL"

// generate regressive alignments ?
params.regressive_align = true

// create progressive alignments ?
params.progressive_align = false

//run dynamic alignments
params.dynamic_align = true
params.dynamicAlnList=["psicoffee_msa", "clustalo_msa", "famsa_msa"]
params.dynamicMaxNseqList=[20, 10000, 1000000]
params.dynamicX='10000'
//params.db="${baseDir}/db/uniref50.fasta"
params.db='/users/cn/egarriga/datasets/db/uniref50.fasta'


// evaluate alignments ?
params.evaluate = true

// bucket sizes for regressive algorithm
params.buckets= '1000'

// output directory
//defined in nextflow.config

log.info """\
         F  A  M  S  A    P  i  p  e  l  i  n  e  ~  version 0.1"
         ======================================="
         Input sequences (FASTA)                        : ${params.seqs}
         Input references (Aligned FASTA)               : ${params.refs}
         Input trees (NEWICK)                           : ${params.trees}
         Alignment methods                              : ${params.align_method}
         Tree methods                                   : ${params.tree_method}
         Generate progressive alignments                : ${params.progressive_align}
         Generate regressive alignments (DPA)           : ${params.regressive_align}
         Bucket Sizes for regressive alignments         : ${params.buckets}
         Perform evaluation? Requires reference         : ${params.evaluate}
         Output directory (DIRECTORY)                   : ${params.outdir}
         """
         .stripIndent()


// Channels containing sequences
if ( params.seqs ) {
  Channel
  .fromPath(params.seqs)
  .map { item -> [ item.baseName, item] }
  .into { seqsCh; seqs2 }
}

if ( params.refs ) {
  Channel
  .fromPath(params.refs)
  .map { item -> [ item.baseName, item] }
  .set { refs }
}

// Channels for user provided trees or empty channel if trees are to be generated [OPTIONAL]
if ( params.trees ) {
  Channel
    .fromPath(params.trees)
    .map { item -> [ item.baseName.tokenize('.')[0], item.baseName.tokenize('.')[1], item] }
    .set { trees }
}
else { 
  Channel
    .empty()
    .set { trees }
}

tree_methods = params.tree_method
align_methods = params.align_method

if ( params.dynamicAlnList ) {
  Channel
  .from(params.dynamicAlnList)
  //.view()
  .set { aligners }
}

if ( params.dynamicMaxNseqList ) {
  Channel
  .from(params.dynamicMaxNseqList)
  //.view()
  .set { maxnseq }
}

aligners
    .merge( maxnseq )
    //.view()
    .set {dynamicMap}
    
/*
 * GENERATE GUIDE TREES USING MEHTODS DEFINED WITH "--tree_method"
 *
 * NOTE: THIS IS ONLY IF GUIDE TREES ARE NOT PROVIDED BY THE USER
 * BY USING THE `--trees` PARAMETER
 */

process generate_trees {
    tag "${id}.${tree_method}"
    publishDir "${params.outdir}/guide_trees", mode: 'copy', overwrite: true
   
    input:
    set val(id), \
         file(seqs) \
         from seqsCh

    each tree_method from tree_methods.tokenize(',')

   output:
     set val(id), \
       val(tree_method), \
       file("${id}.${tree_method}.dnd") \
       into treesGenerated

   when:
     !params.trees

   script:
    template "tree/generate_tree_${tree_method}.sh"
}

treesGenerated
  .mix ( trees )
  .combine ( seqs2, by:0 )
  .set {seqsAndTreesForDynamicAlignment }

process dynamic_msa {
    tag "${id}-${align_method}-${bucket_size}-${dynamic_size}-${tree_method}"
    publishDir "${params.outdir}/alignments", mode: 'copy', overwrite: true

    input:
      set val(id), \
        val(tree_method), \
        file(guide_tree), \
        file(seqs) \
        from seqsAndTreesForDynamicAlignment
      
      val (dynamic_size) from params.dynamicX

      each bucket_size from params.buckets.tokenize(',')
      each align_method from align_methods.tokenize(',') 

    when:
      params.dynamic_align

    output:
      set val(id), \
        val("${align_method}"), \
        val(tree_method), \
        val("reg_align"), \
        val(bucket_size), \
        file("${id}.reg_align.${bucket_size}.dynamic.${dynamic_size}.${align_method}.with.${tree_method}.tree.aln") \
        into dynamicOut

  script:
  """
  ## echo 'psicoffee_msa 20' > config.txt
  ## echo 'clustalo_msa 5000' >> config.txt
  ## echo 'famsa_msa 100000000' >> config.txt
  ## -dynamic_config $PWD/config.txt \
  
  t_coffee -reg -reg_method dynamic_msa \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -dynamic ${dynamic_size} \
         -reg_homoplasy \
         -blast_server LOCAL \
         -protein_db=${params.db} \
         -outfile ${id}.reg_align.${bucket_size}.dynamic.${dynamic_size}.${align_method}.with.${tree_method}.tree.aln
  """
}

refs
  .cross (dynamicOut )
  .map { it -> [it[0][0], it[1][1], it[1][2], it[1][3], it[1][4], it[1][5], it[0][1]] }
  .into { toEvaluate; toEvaluate2}
process evaluation {
    tag "${id}.${align_method}.${tree_method}.${align_type}.${bucket_size}"
    publishDir "${params.outdir}/individual_scores", mode: 'copy', overwrite: true

    input:
      set val(id), \
          val(align_method), \
          val(tree_method), \
          val(align_type), \
          val(bucket_size), \
          file(test_alignment), \
          file(ref_alignment) \
          from toEvaluate

    output:
      set val(id), \
          val(tree_method), \
          val(align_method), \
          val(align_type), \
          val(bucket_size), \
          file("*.sp"), \
          file("*.tc"), \
          file("*.col") \
          into scores

    when:
      params.evaluate

     script:
     """
       t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode sp \
            | grep -v "seq1" | grep -v '*' | \
            awk '{ print \$4}' ORS="\t" \
            > "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.sp"

       t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode tc \
            | grep -v "seq1" | grep -v '*' | \
            awk '{ print \$4}' ORS="\t" \
            > "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.tc"

       t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode column \
            | grep -v "seq1" | grep -v '*' | \
              awk '{ print \$4}' ORS="\t" \
            > "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.col"

    """
}

workflow.onComplete {
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' } runName: ${workflow.runName}"
}
