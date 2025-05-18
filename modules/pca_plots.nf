#!/usr/bin/env nextflow

process PCA_PLOTS {
  time 2.h
  memory 40.GB
  cpus 4
  
  tag "${generation} & ${treatment}"

  publishDir "${params.outdir}/manuscript/supplementary/", mode: 'copy', pattern: "*.pdf"
  
  container 'library://flerpan01/singularity-r/bulk-rrbs:1.1'

  input:
  path metadata
  path coveragefiles
  each generation
  each treatment

  output:
  path "*.pdf"

  """
  PCA_plots.R \
    ${generation} \
    ${treatment}

  rm Rplots.pdf
  """
}