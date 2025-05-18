#!/usr/bin/env nextflow

process SUPPLEMENTARY_PLOTS {
  time 12.h
  //memory 16.GB

  tag "${genomic_feature}"

  //script takes to long time and scratch disk gets busy/unavailble after awhile
  //https://www.nextflow.io/docs/latest/reference/process.html#scratch
  scratch false

  container 'library://flerpan01/singularity-r/bulk-rrbs:1.1'

  publishDir "${params.outdir}/manuscript/supplementary/", mode: 'copy', pattern: "*.pdf"

  input:
  path metadata
  path DMR_TABLES
  path BETAVALUES
  each genomic_feature

  output:
  path "*.pdf"

  """
  supplementary_plots.R ${genomic_feature}
  """
}

process SUPPLEMENTARY_EXCEL {
  time 2.h
  //memory 16.GB

  tag "${genomic_feature}"

  //script takes to long time and scratch disk gets busy/unavailble after awhile
  //https://www.nextflow.io/docs/latest/reference/process.html#scratch
  scratch false

  container 'library://flerpan01/singularity-r/bulk-rrbs:1.1'

  publishDir "${params.outdir}/manuscript/supplementary/", mode: 'copy', pattern: "*.xlsx"

  input:
  path DMR_TABLES
  each genomic_feature

  output:
  path "*.xlsx"

  """
  supplementary_excel.R ${genomic_feature}
  """
}