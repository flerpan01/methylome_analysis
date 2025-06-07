process METHYLKIT {
  tag "Diff. meth. for ${generation}, ${treatment}, ${genomic_feature}"
  //label
  
  conda "${moduleDir}/environment.yml"
  container 'library://flerpan01/singularity-r/bulk-rrbs:1.2'

  time 2.h
  memory 40.GB
  cpus 4
  
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