#!/usr/bin/env Rscript

cat(paste(
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
  "Table wrapper, concatinates all tables into a list\n",
  " > DMR tables\n",
  " > DMG tables\n",
  " > Betavalues tables",
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
))

# ~~ Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

table_wrapper <- function(tabletype){
	cat(" >", tabletype, "\n")

	# build list of betavalues matrix
	files <- list.files(
	  path = ".",
	  pattern = tabletype,
	  recursive = TRUE,
	  full.names = TRUE
	)

	if (tabletype == "betavalues"){
		
		genomic_features <- sapply(files, \(file){
			file <- sub(paste0(tabletype, "_"), "", basename(file))
			file <- sub("(_F[0-9]).*", "", file)
		})
		genomic_features <- unique(genomic_features)
		names(genomic_features) <- genomic_features

		dat <- lapply(genomic_features, \(genomic_feature){
			betavalue_files <- grep(genomic_feature, files, value = TRUE)
			mat <- lapply(betavalue_files, readRDS)
			gen <- sapply(mat, \(x){
		  	unique(sub("(F[0-9]).*", "\\1", colnames(x)))
			})
			names(mat) <- gen

			return(mat)
		})
	
	} else {

		dat <- lapply(files, readRDS)
		names(dat) <- sapply(dat, \(file){
			unique(file$feature)
		})
	}

	filename <- paste0(tabletype, ".Rds")

	saveRDS(dat, filename)
}

# ~~ Variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

tables <- c("DMR", "DMG", "betavalues", "gene_ontology")

# ~~ Code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

cat(" Building table for\n\n")
for (tabletype in tables) table_wrapper(tabletype)

cat("\n~~ wrapper.R complete ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")