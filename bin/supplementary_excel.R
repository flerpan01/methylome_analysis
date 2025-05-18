#!/usr/bin/env Rscript

arg <- commandArgs(trailingOnly = TRUE)
# arg1 = genomic_feature:       islands, promoter, cpg
genomic_feature=arg[1]

cat(paste(
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
  "Genomic feature:\t", genomic_feature, "\n",
  " > Make excel w/ all DMRs and their metadata",
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
))

suppressPackageStartupMessages({
  library(scales)
  library(openxlsx)
  library(gtools)
})

# ~~ Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

extract_dmr_info <- function(dmr){

  n <- length(dmrs) - which(dmrs %in% dmr)
  cat(paste("\n", comma(n), dmr)) # for testing

  tmp <- subset(meth, dmr_id %in% dmr & sign %in% TRUE)
  cols <- c("pvalue", "qvalue", "meth.diff")
  tmp[, cols] <- lapply(tmp[, cols], signif, 3)

  cols <- c("dmr_id", "region", "cgi", "gene_name", "gene_id", "gene_info", "gene_type")

  data.frame(
    tmp[1, cols],
    
    gen = paste(tmp$gen, collapse = "\n"),
    sign = paste(tmp$treatment, collapse = "\n"),
    pvalue = paste(tmp$pvalue, collapse = "\n"),
    qvalue = paste(tmp$qvalue, collapse = "\n"),
    meth.diff = paste(tmp$meth.diff, collapse = "\n")
  )
}

# ~~ Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

files <- list.files(
  ".",
  pattern = paste0("DMR_table_", genomic_feature),
  recursive = TRUE,
  full.names = TRUE
)

# Meth table
meth <- readRDS(files)

# ~~ Code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Only use significant DMRs
dmrs <- unique(meth$dmr_id[meth$sign])   # 8,947
dmrs <- mixedsort(dmrs)

# Make an excel table
excel <- lapply(dmrs, extract_dmr_info)
excel <- Reduce(function(x,y) rbind(x,y), excel)

filename_excel <- paste0("DMR_", genomic_feature, ".xlsx")
write.xlsx(excel, file = filename_excel)

cat(paste(
  "\n\n~~ supplementary_excel.R complete! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n",
  "Generated files:", filename_excel,
  "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
))