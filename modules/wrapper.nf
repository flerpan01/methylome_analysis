#!/usr/bin/env nextflow

process WRAPPER {
  time 2.h
  memory 8.GB

  publishDir "${params.outdir}/data", mode: 'copy', pattern: '*.Rds'
  container 'library://flerpan01/singularity-r/bulk-rrbs:1.1'
  
  input:
  path DMR_TABLES
  path DMG_TABLES
  path BETAVALUES
  path GO_TABLES

  output:
  path '*.Rds'

  """
  wrapper.R
  """
}