#!/usr/bin/env nextflow

process METHYLKIT {
  time 2.h
  memory 40.GB
  cpus 4
  
  tag "${generation} & ${treatment} & ${genomic_feature}"

  //publishDir "${params.outdir}/stub",  mode: 'copy', pattern: "*.Rds"

  container 'library://flerpan01/singularity-r/bulk-rrbs:1.1'

  input:
  path metadata
  path coveragefiles
  path cpgislands_GRCm39
  path refseq_UCSC_GRCm39
  each generation
  each treatment
  each genomic_feature
  val diffmeth_test
  val overdisp
  val eff
  val multiple_test_corr

  output:
  path "diffmeth_${genomic_feature}_${generation}_DBP${treatment}.Rds"

  """
  methylKit.R \
    ${generation} \
    ${treatment} \
    ${genomic_feature} \
    ${diffmeth_test} \
    ${overdisp} \
    ${eff} \
    ${multiple_test_corr}
  """
}