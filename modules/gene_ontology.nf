#!/usr/bin/env nextflow

process GENE_ONTOLOGY {
  time 2.h
  memory 8.GB
  container 'library://flerpan01/singularity-r/gene-ontology:1.0'
  tag "${genomic_feature}"

  //publishDir "${params.outdir}/data", mode: 'copy', pattern: "*.Rds"
  
  input:
  path DMR_TABLES
  path ensembl_dataset
  each genomic_feature

  output:
  path "*.Rds", emit: GO_TABLES

  """
  gene_ontology.R ${genomic_feature}
  """
}