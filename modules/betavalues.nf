#!/usr/bin/env nextflow

process BETAVALUES {
  time 2.h
  memory 20.GB
  cpus 4
  
  tag "${genomic_feature} in gen: ${generation}"

  //publishDir "${params.outdir}/data", mode: 'copy'
  container 'library://flerpan01/singularity-r/bulk-rrbs:1.1'
  
  input:
  path coveragefiles
  path metadata
  path cpgislands_GRCm39
  path refseq_UCSC_GRCm39
  each generation
  each genomic_feature

  output:
  path "betavalues_${genomic_feature}_${generation}.Rds", emit: BETAVALUES

  """
  betavalues.R ${generation} ${genomic_feature}
  """
}