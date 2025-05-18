#!/usr/bin/env nextflow

process DMR_TABLE {
  time 2.h
  memory 8.GB
  container 'library://flerpan01/singularity-r/bulk-rrbs:1.1'

  tag "${genomic_feature}"

  //publishDir "${params.outdir}/data",   mode: 'copy', pattern: "*.Rds"
  publishDir "${params.outdir}",        mode: 'copy', pattern: "*.txt"
  
  input:
  path DMR_files
  path ensembl_dataset
  path cpgislands_GRCm39
  path refseq_UCSC_GRCm39
  each genomic_feature

  output:
  path "*.Rds",   emit: DMR_TABLES
  path "*.txt",   emit: DMR

  """
  DMR_table.R ${genomic_feature}
  """
}

process DMG_TABLE {
  time 8.h
  memory 8.GB
  container 'library://flerpan01/singularity-r/bulk-rrbs:1.1'
  scratch false

  tag "filetype: ${DMG_table_output} & ${genomic_feature}"

  //publishDir "${params.outdir}/data", mode: 'copy', pattern: "*.Rds"
  publishDir "${params.outdir}/manuscript/supplementary/", mode: 'copy', pattern: "*.xlsx"

  input:
  path DMR_TABLE
  path ensembl_dataset
  each DMG_table_output
  each genomic_feature

  output:
  path "*.Rds",   emit: DMG_TABLES,  optional: true
  path "*.xlsx",  optional: true

  """
  DMG_table.R ${DMG_table_output} ${genomic_feature}
  """
}