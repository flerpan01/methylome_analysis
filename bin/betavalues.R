#!/usr/bin/env Rscript

# arg1 = generations
# arg2 = genomic feature
args <- commandArgs(trailingOnly = TRUE)
gen=args[1];genomic_feature=args[2]

cores=4;treatment=c(10,100)

cat(paste(
	"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
	"Generating betavalues matrix\n",
	"\n % % % % % Setup % % % % %",
	"\n",
	"\n Generation:\t\t", gen,
	"\n Treatment:\t\t", paste(treatment, collapse = " & "),
	"\n Genomic feature:\t", genomic_feature,
	"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
))

suppressPackageStartupMessages({
  library(methylKit)
  library(gtools)
  library(scales)

  if ( genomic_feature %in% c("islands", "promoter") ) library(genomation)
})

# ~~ Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

build_methylraw_obj <- function(gen, treatment, genomic_feature){
	
	metadata <- meta[meta$gen %in% gen & meta$treatment %in% c(0, treatment), ]
	files <- coverage_files[names(coverage_files) %in% metadata$id]
	
	cat("\n\n Reading samples...\n")

	obj <- methRead(
	  location = split(files, seq_along(files)),
	  sample.id = split(names(files), seq_along(files)),
	  assembly = "GRCm39",

	  # "amp", "bismark","bismarkCoverage", "bismarkCytosineReport"
	  pipeline = "bismarkCoverage",
	  header = FALSE,
	  treatment = metadata$treatment,
	)

	cat("\n\n Filtering and Normalising...\n")

	# Filter reads
	min_coverage <- ifelse(genomic_feature %in% c("islands", "promoter"), 10, 20)
	
	obj <- filterByCoverage(
		obj, 
		lo.count = min_coverage,
		hi.perc = 99.9 						# the top 99th percentile (PCR duplicates)
	)

	# Normalisation, scaling factor between samples based on differences
	# between median of coverage distribution
	obj <- normalizeCoverage(
	  obj,
	  method = "median"
	)

	# add prefix 'chr' to chromosome column
	for (i in seq_along(obj)){
		obj[[i]]$chr <- paste0("chr", obj[[i]]$chr)
	}
	
	return(obj)
}

build_methylation_matrix <- function(obj, genomic_feature){
	if (genomic_feature == "islands"){
		cat("\n\n Converting genomic regions to CpG-islands\n")
		
		# overlap meth data with gene feature data
	  cgi <- readGeneric("cpgislands_GRCm39.bed", keep.all.metadata = TRUE)
		obj <- regionCounts(obj, regions = cgi)

	} else if (genomic_feature == "promoter"){
		cat("\n\n Converting genomic regions to promoters\n")
		
		# overlap meth data with gene feature data
	  refseq <- readTranscriptFeatures("refseq_UCSC_GRCm39.bed")
		obj <- regionCounts(obj, regions = refseq$promoters)
		
	} else {
		cat("\n\n Using CpG-sites\n")

	}

	#N <- ceiling(min(table(attributes(obj)$treatment)) / 2)
	N <- 1

	cat("\n\n Minimum number of samples per group:", N, "\n\n")

	meth <- methylKit::unite(
		obj, 
		min.per.group = as.integer(N),
		mc.cores = cores
	)

	cat("\n\n Saving betavalues...\n")
  
  dmr_id <- paste0(getData(meth)$chr, ":", getData(meth)$start)

  percMat <- percMethylation(meth)
  rownames(percMat) <- dmr_id

  return(percMat)
}

# ~~ Code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# locate all coverage files
coverage_files <- list.files(
	path = ".",
	pattern = "*.cov.gz",
	recursive = TRUE,
	full.names = TRUE
)

# extract sample id
names(coverage_files) <- sapply(coverage_files, \(path){
	id <- basename(path)
	sub("(F[0-9]_[0-9]+).*", "\\1", id)
})

# order files
coverage_files <- coverage_files[mixedorder(names(coverage_files))]

meta <- read.csv("metadata.csv")

obj <- build_methylraw_obj(gen, treatment, genomic_feature)
percMat <- build_methylation_matrix(obj, genomic_feature)

filename <- paste0("betavalues_", genomic_feature, "_", gen, ".Rds")

saveRDS(percMat, file = filename)

cat(paste(
	"\n~~ betavalues.R complete ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
	filename, "generated",
	"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
))