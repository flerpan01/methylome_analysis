#!/usr/bin/env Rscript

arg <- commandArgs(trailingOnly = TRUE)

# methylKit parameters
# arg1 = genomic_feature:       islands, promoter, cpg
genomic_feature=arg[1]

cat(paste(
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
  "Genomic feature:", genomic_feature, "\n",
  "Generating two tables\n",
  " > table of all the DMRs (Holy Grail)\n",
  " > Table based on each gene with the DMR info included in a column",
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
))

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(genomation)
})

# ~~ Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

build_DMRtable <- function(meth){
  # make genomic ranges objects
  cols <- c("chr", "start", "end", "strand")
  gr_meth <- as(meth[, cols], "GRanges")
  gr_ens <- as(ens[, cols], "GRanges")
  
  # add 2000bp upstream each TSS to include the promoter region in the gene
  rows <- as.character(strand(gr_ens)) == "+"
  start(gr_ens)[rows] <- start(gr_ens)[rows] - 2000
  rows <- as.character(strand(gr_ens)) == "-"
  end(gr_ens)[rows] <- end(gr_ens)[rows] + 2000

  # overlap meth data and ensembl gene data
  rows <- data.frame(findOverlaps(gr_ens, gr_meth))

  # add gene ID
  meth$gene_id <- NA
  meth$gene_id[rows$subjectHits] <- ens$gene_id[rows$queryHits]

  # add gene name
  meth$gene_name <- NA
  meth$gene_name[rows$subjectHits] <- ens$gene_name[rows$queryHits]

  # add gene info
  meth$gene_info <- NA
  meth$gene_info[rows$subjectHits] <- ens$gene_info[rows$queryHits]

  # add gene type
  meth$gene_type <- NA
  meth$gene_type[rows$subjectHits] <- ens$gene_type2[rows$queryHits]

  # overlap meth data with gene feature data
  cgi <- readGeneric("cpgislands_GRCm39.bed", keep.all.metadata = TRUE)

  # assign each cgi an ID, numeric
  values(cgi) <- DataFrame(id = 1:length(cgi), values(cgi))

  # https://genome-euro.ucsc.edu/cgi-bin/hgTables
  anno <- readTranscriptFeatures("refseq_UCSC_GRCm39.bed")
  anno[["cgi"]] <- cgi

  # add info about genome regions
  # promoter > exon > intron > intergenic
  
  # intergenic
  meth$region <- "intergenic"

  # introns
  rows <- data.frame(findOverlaps(anno$intron, gr_meth))
  meth$region[rows$subjectHits] <- "intron"

  # exons
  rows <- data.frame(findOverlaps(anno$exon, gr_meth))
  meth$region[rows$subjectHits] <- "exon"

  # promoters
  rows <- data.frame(findOverlaps(anno$promoters, gr_meth))
  meth$region[rows$subjectHits] <- "promoter"

  # intergenic
  #rows <- is.na(meth$gene_id)
  #meth$region[rows] <- "intergenic"

  # add cpgisland data, include as ID
  rows <- data.frame(findOverlaps(anno$cgi, gr_meth))
  meth$cgi <- NA
  meth$cgi[rows$subjectHits] <- anno$cgi$id[rows$queryHits]

  return(meth)
}

# ~~ Variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

methdiff <- ifelse(genomic_feature %in% c("islands", "promoter"), 10, 30)

# ~~ Code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# import files
pattern <- paste0(genomic_feature, "_F[0-9].+.Rds")
files <- list.files(
  ".", # expects an Rds-file from each generation and treatment
  pattern = pattern,
  recursive = TRUE,
  full.names = TRUE
)
meth <- lapply(files, readRDS)
meth <- Reduce(function(x,y) rbind(x,y), meth)

# add signifiance threshold
meth$sign <- meth$qvalue <= 0.05 & abs(meth$meth.diff) >= methdiff

ens <- read.csv("ensembl-dataset.csv.gz")
ens <- ens[!duplicated(ens$gene_id), ]    # only keep unique gene IDs

# add genomic features
meth <- build_DMRtable(meth)

filename <- paste0("DMR_table_", genomic_feature, ".Rds")
saveRDS(meth, file = filename)

tmp <- subset(meth, sign == TRUE)
dat <- table(tmp$region, tmp$gene_type)
out <- as.matrix.data.frame(dat)
out <- data.frame(out)
colnames(out) <- colnames(dat)
out$region <- factor(rownames(dat), levels = c("promoter", "exon", "intron", "intergenic"))
cols <- factor(colnames(out), c("region", "protein_coding", "lncRNA", "ncRNA",
   "miRNA", "rRNA", "pseudogene", "ribozyme", "TEC"))
out <- out[order(out$region), ]
out <- out[, order(cols)]

dmr_filename <- paste0("DMR_", genomic_feature, ".txt")

write.table(
  out,
  file = dmr_filename,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

cat(paste(
  "\n~~ DMR_table.R complete ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
  filename, "generated",
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
))