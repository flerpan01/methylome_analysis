#!/usr/bin/env Rscript

# expects argument for type of file to save
arg <- commandArgs(trailingOnly = TRUE)
output=arg[1];genomic_feature=arg[2]
#output="Rds";genomic_feature="cpg"

cat(paste(
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
  "Generating table based on DMRs located inside genes\n",
  " > Each row is a unique gene\n",
  " > Using genomic feature:", genomic_feature, "\n",
  " > Saving output as:", output, "file",
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
))

suppressPackageStartupMessages({
  if (output == "excel") library(openxlsx)
})

# ~~ Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

make_genetable <- function(gene, output){
  # output = Rds, excel

  # get gene info
  cols <- c("gene_name", "chr", "gene_id", "gene_info", "gene_type")
  geneinfo <- ens[ens$gene_id %in% gene, cols]

  # region of interests = promoter, exon
  region <- c("promoter", "exon", "intron", "intergenic")

  # get DMR info
  rows <- DMR$gene_id %in% gene & DMR$sign & DMR$region %in% region
  dmrinfo <- DMR[rows, ]

  # make gene table
  out <- data.frame(geneinfo, num_DMR = nrow(dmrinfo))

  if (output == "excel"){
    if (nrow(dmrinfo) > 0){
      out$meth <- paste(dmrinfo$type, collapse = "\n")
      out$coef <- paste(paste0(dmrinfo$gen, "_DBP", dmrinfo$treatment), collapse = "\n")
      out$promoter <- paste(dmrinfo$dmr_id[dmrinfo$region %in% "promoter"], collapse = "\n")
      out$exon <- paste(dmrinfo$dmr_id[dmrinfo$region %in% "exon"], collapse = "\n")
      out$intron <- paste(dmrinfo$dmr_id[dmrinfo$region %in% "intron"], collapse = "\n")
      out$intergenic <- paste(dmrinfo$dmr_id[dmrinfo$region %in% "intergenic"], collapse = "\n")
      out$cgi <- paste(dmrinfo$dmr_id[!is.na(dmrinfo$cgi)], collapse = "\n")
    } else {
      out$meth <- ""
      out$coef <- ""
      out$promoter <- ""
      out$exon <- ""
      out$intron <- ""
      out$intergenic <- ""
      out$cgi <- ""
    }
  }

  if (output == "Rds"){
    if (nrow(dmrinfo) > 0){
      out$meth <- list(dmrinfo$type)
      out$coef <- list(unique(paste0(dmrinfo$gen, "_DBP", dmrinfo$treatment)))
      out$promoter <- list(dmrinfo$dmr_id[dmrinfo$region %in% "promoter"])
      out$exon <- list(dmrinfo$dmr_id[dmrinfo$region %in% "exon"])
      out$intron <- list(dmrinfo$dmr_id[dmrinfo$region %in% "intron"])
      out$intergenic <- list(dmrinfo$dmr_id[dmrinfo$region %in% "intergenic"])
      out$cgi <- list(dmrinfo$dmr_id[!is.na(dmrinfo$cgi)])
    } else {
      out$meth <- NA
      out$coef <- NA
      out$promoter <- NA
      out$exon <- NA
      out$intron <- NA
      out$intergenic <- NA
      out$cgi <- NA
    }
  }

  return(out)
}

# ~~ Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

files <- list.files(
  ".",
  pattern = paste0("DMR_table_", genomic_feature),
  recursive = TRUE,
  full.names = TRUE
)

DMR <- readRDS(files)
ens <- read.csv("ensembl-dataset.csv.gz")
ens <- ens[!duplicated(ens$gene_id), ]    # only keep unique gene IDs

# ~~ Variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Only use genes w/ DMRs
genes <- unique(ens$gene_id)
genes <- genes[genes %in% na.omit(unique(DMR$gene_id[DMR$sign]))]

# ~~ Code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

out <- lapply(genes, \(gene) make_genetable(gene, output))
out <- Reduce(function(x,y) rbind(x,y), out)

if (output == "excel"){
  filename <- paste0("DMG_", genomic_feature, ".xlsx")
  write.xlsx(out, file = filename)
}

if (output == "Rds"){
  filename <- paste0("DMG_", genomic_feature, ".Rds")
  saveRDS(out, file = filename)
}

cat(paste(
  "\n~~ DMG_table.R complete ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
  filename, "generated",
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
))