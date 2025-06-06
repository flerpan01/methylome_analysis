process ENSEMBL {
  tag "Ref genome: ${reference_genome}"
  //label
  
  conda "${moduleDir}/environment.yml"
  container 'library://flerpan01/singularity-r/bulk-rrbs:1.2'
  
  time 1.h
  memory 8.GB
  cpus 1

  input:
  val reference_genome

  output:
  path 'ensembl_dataset.csv.gz', emit: ENSEMBL_DATASET

  script:
  //def args = task.ext.args ?: ''
  //args =+ // HÄR ÄR DU
  """
  # set environment variables for biomart cache
  export BIOMART_CACHE="/scratch/.cache"

  ensembl.R ${reference_genome}
  """
}