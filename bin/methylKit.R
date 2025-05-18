#!/usr/bin/env Rscript

# ~~ arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
arg <- commandArgs(trailingOnly = TRUE)

# experimental setup
# arg1 = generations
# arg2 = treatment
# arg3 = genomic_feature: 			islands, promoter, cpg
gen=arg[1];treatment=arg[2];genomic_feature=arg[3]

# methylKit parameters
# arg4 = differantial test: 		F, Chisq
# arg5 = overdispersion: 				none, MN, shrinkMN
# arg6 = effect: 								wmean, mean, predicted
# arg7 = multiple test. corr.
diffmeth_test=arg[4];overdisp=arg[5];eff=arg[6];multiple_test_corr=arg[7]

cores=4

# ~~ testing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#gen="F1";treatment=10;cores=1;diffmeth_test="Chisq";
#overdisp="MN";eff="mean";multiple_test_corr="BH";genomic_feature="islands"

cat(paste(
	"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
	"Calculating differentially methylated regions w/ methylKit\n",
	"\n % % % % %  ~  A R G U M E N T S  ~  % % % % %\n",
	"\n Generation:\t\t", gen,
	"\n Treatment:\t\t", treatment,
	"\n Genomic feature:\t", genomic_feature,
	"\n",
	"\n Using # cores:\t\t", cores,
	"\n Diff. meth. test:\t", diffmeth_test,
	"\n Overdispersion:\t", overdisp,
	"\n Mean meth. diff:\t", eff,
	"\n Multiple test. corr.:\t", multiple_test_corr,
	"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
))

suppressPackageStartupMessages({
  library(methylKit)
  library(gtools)
  library(scales)
  
  if ( genomic_feature %in% c("islands", "promoter") ) library(genomation)
})

# ~~ Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

check_outliers <- function(x, diffsize = 3){
	# takes in a vector
	# runs scales
	# removes any scaled abs(value) > 3
	# reruns again
	# returns vector w/o outliers

	while ( any( (abs(scale(x))) > diffsize, na.rm = TRUE ) ){
		#x <- x[abs(scale(x)) < diffsize]
		row <- which( abs(scale(x)) > diffsize )
		x[row] <- NA
	}
	return(x)
}

build_methylraw_obj <- function(gen, treatment, genomic_feature) {
	
	metadata <- meta[meta$gen %in% gen & meta$treatment %in% c(0, treatment), ]
	files <- coverage_files[names(coverage_files) %in% metadata$id]
	
	cat("\n Reading samples...\n\n")

	obj <- methRead(
	  location = split(files, seq_along(files)),
	  sample.id = split(names(files), seq_along(files)),
	  assembly = "GRCm39",

	  # "amp", "bismark","bismarkCoverage", "bismarkCytosineReport"
	  pipeline = "bismarkCoverage",
	  header = FALSE,
	  treatment = metadata$treatment
	)

	# Un-annotated chromosomes contains . (dot), remove these
	# add prefix 'chr' to chromosome column
	for (i in seq_along(obj)){
		rows <- grepl("\\.", obj[[i]]$chr)
		obj[[i]] <- obj[[i]][!rows,]

		obj[[i]]$chr <- paste0("chr", obj[[i]]$chr)
	}

	cat("\n Filtering and Normalising...\n\n")

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

	return(obj)
}

diff_meth_analysis <- function(obj, genomic_feature){
	
	if (genomic_feature == "islands"){
		cat("\n\n Converting genomic regions to CpG-islands\n")
		
		# overlap meth data with gene feature data
	  cgi <- readGeneric("cpgislands_GRCm39.bed", keep.all.metadata = TRUE)
		obj <- regionCounts(obj, regions = cgi)

		methdiff <- 10
	} else if (genomic_feature == "promoter"){
		cat("\n\n Converting genomic regions to promoters\n")
		
		# overlap meth data with gene feature data
	  refseq <- readTranscriptFeatures("refseq_UCSC_GRCm39.bed")
		obj <- regionCounts(obj, regions = refseq$promoters)
		
		methdiff <- 10
	} else {
		cat("\n\n Using CpG-sites\n")

		methdiff <- 25
	}

	N <- ceiling(min(table(attributes(obj)$treatment)) / 2)

	cat("\n\n Minimum number of samples per group:", N, "\n\n")

	meth <- methylKit::unite(
		obj, 
		min.per.group = as.integer(N),
		mc.cores = cores
	)
	
	# remove regions with low variation, std < 2
	mat <- percMethylation(meth)
	std <- matrixStats::rowSds(mat, na.rm = TRUE)
	meth <- meth[std > 2]

	cat(paste(
		"\n",
		"number regions post filtering:", 
		round(sum(std > 2) / length(std) * 100, 2), "%",
		"(", comma(sum(std > 2)), "/", comma(length(std)), ")",
		"\n"
	))

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
	# dmr="chr1:4258007"
	#mat <- percMethylation(meth)
	#rownames(mat) <- paste0("chr", getData(meth)$chr, ":", getData(meth)$start)
	#mat2 <- percMethylation(meth)[1:100,]
	#mat2[] <- apply(mat2, 1, check_outliers)
	# 1. Get the data frame
	#mb_df <- getData(methylBase_obj)

	# 2. Set values to NA (example: set all values in sample2 with coverage < 5 to NA)
	#mb_df$numCs2[mb_df$coverage2 < 5] <- NA
	#mb_df$numTs2[mb_df$coverage2 < 5] <- NA

	# 3. Reconstruct the object
	#new_methylBase <- methylKit::methylKitObject(
	#  mb_df,
	#  sample.ids = methylBase_obj@sample.ids,
	#  assembly = methylBase_obj@assembly,
	#  context = methylBase_obj@context,
	#  treatment = methylBase_obj@treatment
	#)

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

	#tmp <- pool(meth, sample.ids = c("0", "10"))
	#q <- percMethylation(tmp)
	#rownames(q) <- paste0("chr", tmp$chr, ":", tmp$start)
	
	#foo <- function(treatment){
	#	cols <- which(meth@treatment %in% treatment)
	#	cols <- lapply(cols, \(col) grep(paste0("num[a-zA-Z]+", col,"$"), names(meth), value = TRUE))
	#	unlist(cols)
	#}
	#sum(getData(meth)[1, foo(0)], na.rm=T)
	#sum(getData(meth)[1, foo(10)], na.rm=T)

	#dat <- calculateDiffMeth(
  #	tmp,
  #  overdispersion = overdisp,
  #  adjust = multiple_test_corr,
  #  mc.cores = cores
  #)


	cat("\n\n Calculating DMR...\n\n")

	dat <- calculateDiffMeth(
  	meth,
    overdispersion = overdisp,
    adjust = multiple_test_corr,
    effect = eff,
    test = diffmeth_test,
    mc.cores = cores
  )

  # convert to dataframe
  data <- getData(dat)

  # add type and cpg_id
  data$type <- ifelse(data$meth.diff > 0, "hyper", "hypo")
  data$dmr_id <- paste0(data$chr, ":", data$start)
  data$feature <- genomic_feature
  data$gen <- gen
  data$treatment <- treatment

  cat(
		"\n # DMRs found:", 
		comma(sum(data$qvalue < 0.05 & abs(data$meth.diff) > methdiff)),
		"\n # CpG tested:", 
		comma(nrow(data)),
		"\n"
	)

  return(data)
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
meth <- diff_meth_analysis(obj, genomic_feature)

filename <- paste0("diffmeth_", genomic_feature, "_", gen, "_DBP", treatment, ".Rds")

saveRDS(meth, file = filename)

cat(paste(
	"\n~~ methylKit.R complete ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
	filename, " generated",
	"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
))