#!/usr/bin/env nextflow

// ~~ Import processes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
include { PCA_PLOTS           } from './modules/pca_plots.nf'
include { BETAVALUES          } from './modules/betavalues.nf'
include { METHYLKIT           } from './modules/methylKit.nf'
include { DMR_TABLE           } from './modules/diffmeth_tables.nf'
include { DMG_TABLE           } from './modules/diffmeth_tables.nf'
include { GENE_ONTOLOGY       } from './modules/gene_ontology.nf'
include { SUPPLEMENTARY_EXCEL } from './modules/supplementary.nf'
include { SUPPLEMENTARY_PLOTS } from './modules/supplementary.nf'
include { WRAPPER             } from './modules/wrapper.nf'

// ~~ Channels ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
ch_coverage_files     = Channel.fromPath(params.coveragefiles)
ch_metadata           = Channel.fromPath(params.metadata)
ch_ensembl_dataset    = Channel.fromPath(params.ensembl_dataset)
ch_cpgislands_GRCm39  = Channel.fromPath(params.cpgislands_GRCm39)
ch_refseq_UCSC_GRCm39 = Channel.fromPath(params.refseq_UCSC_GRCm39)
generations           = Channel.of( 'F0', 'F1', 'F2' )
treatments            = Channel.of( '10', '100' )
DMG_table_output      = Channel.of( 'excel', 'Rds')
genomic_features      = Channel.of( 'cpg' , 'islands' , 'promoter')

log.info \
  """
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Outdir                    : ${params.outdir}

    methylKit parameters
    --------------------
    Diff. meth. test          : ${params.diffmeth_test}
    Overdispertion            : ${params.overdisp}
    Effect                    : ${params.eff}
    Multiple test. corr.      : ${params.multiple_test_corr}
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  """
  .stripIndent(true)

// ~~ Workflow ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *

workflow {
  
  PCA_PLOTS (
    ch_metadata,
    ch_coverage_files,
    generations,
    treatments
  )

  BETAVALUES (
    ch_coverage_files,
    ch_metadata,
    ch_cpgislands_GRCm39,
    ch_refseq_UCSC_GRCm39,
    generations,
    genomic_features
  )

  METHYLKIT (
    ch_metadata,
    ch_coverage_files,
    ch_cpgislands_GRCm39,
    ch_refseq_UCSC_GRCm39,
    generations,
    treatments,
    genomic_features,
    params.diffmeth_test,
    params.overdisp,
    params.eff,
    params.multiple_test_corr
  )

  DMR_TABLE (
    METHYLKIT.out.collect(),
    ch_ensembl_dataset,
    ch_cpgislands_GRCm39,
    ch_refseq_UCSC_GRCm39,
    genomic_features
  )

  DMG_TABLE (
    DMR_TABLE.out.DMR_TABLES.collect(),
    ch_ensembl_dataset,
    DMG_table_output,
    genomic_features
  )

  GENE_ONTOLOGY (
    DMR_TABLE.out.DMR_TABLES.collect(),
    ch_ensembl_dataset,
    genomic_features
  )

  SUPPLEMENTARY_EXCEL (
    DMR_TABLE.out.DMR_TABLES.collect(),
    genomic_features
  )

  SUPPLEMENTARY_PLOTS (
    ch_metadata,
    DMR_TABLE.out.DMR_TABLES.collect(),
    BETAVALUES.out.BETAVALUES.collect(),
    genomic_features
  )

  WRAPPER ( // concatinate all tables w/ genomic features
    DMR_TABLE.out.DMR_TABLES.collect(),
    DMG_TABLE.out.DMG_TABLES.collect(),
    BETAVALUES.out.BETAVALUES.collect(),
    GENE_ONTOLOGY.out.GO_TABLES.collect()
  )
}

workflow.onComplete {
  def msg = """\
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Pipeline execution summary
    --------------------------
    Completed at     : ${workflow.complete}
    Duration         : ${workflow.duration}
    Success          : ${workflow.success}
    workDir          : ${workflow.workDir}
    exit status      : ${workflow.exitStatus}

    methylKit options
    --------------------------
    Diff. meth. test : ${params.diffmeth_test}
    Overdispertion   : ${params.overdisp}
    Effect           : ${params.eff}
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """
    .stripIndent()

  sendMail (
    to: "${params.email}",
    subject: 'RRBS DBP Spleen analysis',
    body: msg
  )
}