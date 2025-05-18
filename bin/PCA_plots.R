#!/usr/bin/env Rscript

# ~~ arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
arg <- commandArgs(trailingOnly = TRUE)

# experimental setup
# arg1 = generations
# arg2 = treatment
gen=arg[1];treatment=arg[2]

#gen="F0";treatment=10;cores=1
cores=4

cat(paste(
	"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
	"Generating PCA plots\n",
	"\n % % % % %  ~  A R G U M E N T S  ~  % % % % %",
	"\n Generation:\t\t", gen,
	"\n Treatment:\t\t", treatment,

	"\n Using # cores:\t\t", cores,
	"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
))

suppressPackageStartupMessages({
  library(methylKit)
  library(gtools)
  library(scales)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  library(cowplot)
  library(patchwork)
})

# ~~ Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

build_methylraw_obj <- function(gen, treatment){
	metadata <- meta[meta$gen %in% gen & meta$treatment %in% c(0, treatment), ]
	files <- coverage_files[names(coverage_files) %in% metadata$id]
	
	cat("\n\n Reading samples...\n\n")

	obj <- methRead(
	  location = split(files, seq_along(files)),
	  sample.id = split(names(files), seq_along(files)),
	  assembly = "GRCm39",

	  # "amp", "bismark","bismarkCoverage", "bismarkCytosineReport"
	  pipeline = "bismarkCoverage",
	  header = FALSE,
	  treatment = metadata$treatment,
	)

	cat("\n Filtering and Normalising...\n\n")

	# Filter reads
	obj <- filterByCoverage(
		obj, 
		lo.count = 30,
		hi.perc = 99.9 	# the top 99th percentile (PCR duplicates)
	)

	# Normalisation, scaling factor between samples based on differences
	# between median of coverage distribution
	obj <- normalizeCoverage(
	  obj,
	  method = "median"
	)

	return(obj)
}

pca_plot <- function(obj, title = NULL, size = 12){
	res <- PCASamples(obj, obj.return = TRUE)

	data <- data.frame(res$x[,1:2], treatment = obj@treatment)

	data$PC1 <- as.numeric(data$PC1) / res$sdev[1] * sqrt(nrow(data))
	data$PC2 <- as.numeric(data$PC2) / res$sdev[2] * sqrt(nrow(data))
	data$treatment <- as.factor(data$treatment)
	data$id <- rownames(data)

	var_perc <- round(res$sdev^2 / sum(res$sdev^2) * 100, 1)

	# set colors
  colors <- brewer.pal(8, "Set1")[c(3,5,8)]
  names(colors) <- c("0", "10", "100")

	p <- ggplot(data, aes(PC1, PC2, fill = treatment, label = id)) +
		theme_linedraw(size) + background_grid() +
		scale_fill_manual(values = colors) +
	  geom_label_repel(size = 3, show.legend = FALSE) +
	  geom_point(shape = 21, size = 3) +
    labs(
    	title = title,
    	x = paste("PC1:", var_perc[1], "% var."),
      y = paste("PC2:", var_perc[2], "% var"),
    	fill = "")

	return(p)
}

hist_plot <- function(std){

	
	ggplot() +
		theme_linedraw(size) + background_grid() +
		geom_histogram(aes(std[std > 3]), fill = "green") +
		geom_histogram(aes(std), color = "black", alpha = 0.4)
}

# ~~ Variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

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

obj <- build_methylraw_obj(gen, treatment)

N <- ceiling(min(table(attributes(obj)$treatment)) / 2)
cat("\n Minimum number of samples per group:", N, "\n\n")

meth <- methylKit::unite(
	obj, 
	min.per.group = as.integer(N),
	mc.cores = cores
)

p1 <- pca_plot(meth, "Unfiltered")

mat <- percMethylation(meth)
std <- matrixStats::rowSds(mat, na.rm = TRUE)

meth2 <- meth[std > 2]
p2 <- pca_plot(meth2, "Filtered")

filename_geneplots <- paste0("PCA-", gen, "-", treatment, ".pdf")
pdf(file = filename_geneplots, width = 9, height = 4)

plot_grid(p1, p2, nrow = 1)

cat(paste(
	"\n~~ pca.R complete ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
	filename_geneplots, "generated",
	"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
))




#badsamples <- c("F0_17", "F0_7")
#samples <- obj2@sample.ids[!obj2@sample.ids %in% badsamples]
#treatments <- obj2@treatment[!obj2@sample.ids %in% badsamples]
#obj3 <- reorganize(obj2, sample.ids = samples, treatment = treatments)

#p3 <- pca_plot(obj3)
#plot_grid(p1, p2, p3, nrow=1)




build_metadata <- function(dose = "all", gen = NULL) {
	if (is.null(gen)) gen <- c("F0", "F1", "F2")

	files <- list.files(filespath, pattern = "*.gz", full.names = TRUE)
	files <- data.frame(
		files = files,
		sample_id = sub("_trimmed_bismark_bt2.bismark.cov.gz", "", files)
		)
	files$sample_id = sub(paste0(filespath, "/"), "", files$sample_id)
	files$generation <- sapply(files$sample_id, function(x) {
		strsplit(x, "_")[[1]][1]
		})

	metadata <- fread(file.path(dir, "data/dbp_phenotypes.csv"))[, 1:4]
	metadata <- merge(metadata, files)
	
	if (dose != "all") metadata <- subset(metadata, treatment %in% c(0, dose) & generation %in% gen)

	return(metadata)
}

build_methylObj <- function(metadata, tile = NULL, chromosomes, make_plots = FALSE, cores = 20){

	dmr_type <- ifelse(is.null(tile), "cpg", paste0("tile", tile))

	filelist <- split(
		as.character(metadata$files),
		seq_along(metadata$files)
	)

	samplelist <- split(
		as.character(metadata$sample_id),
		seq_along(metadata$sample_id)
	)

	obj <- methRead(
		location = filelist,
		sample.id = samplelist,
		treatment = metadata$treatment,
		header = FALSE,
		assembly = "GRCm39",

		# "amp", "bismark","bismarkCoverage", "bismarkCytosineReport"
		pipeline = "bismarkCoverage"
	)

	# sanity check position against the reference genome
	ref_cg <- fread("~/projdir/GRCm39/cg_motifs_GRCm39.csv.gz")

	ref_cg <- subset(ref_cg, chr %in% chromosomes)
	
								# quick fix
								ref_cg$chr <- sub("chrM", "chrMT", ref_cg$chr) 
								ref_cg$pos <- sub("chrM", "chrMT", ref_cg$pos)

	gr_ref <- as(ref_cg[, 1:3], "GRanges")

	# iterate ref cpgs over all samples
	for (i in seq_along(obj)){
		dat <- getData(obj[[i]])[, 1:3]
		dat$chr <- paste0("chr", dat$chr)

		gr_obj <- as(dat, "GRanges")

		rows <- findOverlaps(gr_ref, gr_obj)
		rows <- attributes(rows)$to

		print(paste(nrow(dat) - length(rows), "cpg-pos removed from sample", i))

		obj[[i]] <- obj[[i]][rows, ]
	}

	if (!is.null(tile)){
		print(paste("calculating tiles w/ size:", tile))
		
		obj <- tileMethylCounts(
			obj,
			win.size = tile,
			step.size = tile,
			cov.bases = 2, 
			mc.cores = 20
		) 
	}

	if (make_plots) {
		name <- paste0("getCoverageStats_", dmr_type, ".pdf")
		filename <- file.path(reportpath, name)

		pdf(filename)
		for (i in seq_along(obj)) {
			print(getCoverageStats(obj[[i]], plot = TRUE))
		}
		dev.off()

		name <- paste0("getMethylationStats_", dmr_type, ".pdf")
		filename <- file.path(reportpath, name)

		pdf(filename)
		for (i in seq_along(obj)) {
			print(getMethylationStats(obj[[i]], plot = TRUE))
		}
		dev.off()
	}

	# Filter reads
	obj <- filterByCoverage(
		obj, 
		lo.count = 10,	# < 10
		hi.perc = 99.9 	# the top 99th percentile (PCR duplicates)
	) 

	# Normalisation, scaling factor between samples based on differences
	# between median of coverage distribution
	obj <- normalizeCoverage(
	  obj,
	  method = "median"
	)

	return(obj)
}

build_methdata <- function(obj, tile = NULL, nsamples = 1L, make_plots = FALSE){

	dmr_type <- ifelse(is.null(tile), "cpg", paste0("tile", tile))

	if (!is.integer(nsamples)){
		nsamples <- as.integer(ceiling(table(metadata$treatment) %>% mean * nsamples))
	} 

	print(paste("Minimum number of samples per group:", nsamples))

	meth <- unite(
		obj, 
		min.per.group = nsamples
	)

	# calc % methylation scores & standard deviation
	mat <- percMethylation(meth)
	std <- matrixStats::rowSds(mat, na.rm = TRUE)

	if (make_plots) {
		# descriptive statistics before std-filtering
		pca_before <- pca_plot(meth)
		hist_before <- hist_plot(std)
	}

	# remove CpG-sites with low variation, std < 2
	meth <- meth[std > 2]

	print(paste("# of CpG-sites passed the filtering:", nrow(meth)))
	print(paste("# of CpG-sites with low variation:", sum(std < 2)))

	if (make_plots) {
		# descriptive statistics after std-filtering
		pca_after <- pca_plot(meth)
		hist_after <- hist_plot(std[std > 2])

		# The plot 2 x 2
		title1 <- ggdraw() +
			draw_label("PCA and hist of std", x = 0.05, hjust = 0, vjust = 1)

		title2 <- ggdraw() +
			draw_label("std < 2 CpG-sites removed", x = 0.05, hjust = 0, vjust = 1)

		top <- plot_grid(pca_before, hist_before, nrow = 1, labels = c("A", "B"))
		bottom <- plot_grid(pca_after, hist_after, nrow = 1, labels = c("C", "D"))

		p <- plot_grid(
			title1,
			top,
			title2,
			bottom,
			ncol = 1,
			rel_heights = c(0.1, 1, 0.1, 1)
		)

		name <- paste0("meth_stats_", dmr_type, ".pdf")
		filename <- file.path(reportpath, name)
		
		pdf(filename, width = 12)
		print(p)
		dev.off()
	}

	return(meth)
}