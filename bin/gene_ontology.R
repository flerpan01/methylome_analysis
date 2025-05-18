#!/usr/bin/env Rscript

arg <- commandArgs(trailingOnly = TRUE)
genomic_feature=arg[1]
#genomic_feature="cpg"

cat(paste(
	"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
	"Gene ontology enrichment analysis\n\n",
  " > Genomic feature:", genomic_feature, "\n\n",
	"Divided between DMRs in:\n",
	" > (1) promoter, (2) exon, (3) CGI\n\n",

	" > hyper-/hypomethylated promoter\n",
	" > hyper-/hypomethylated exon\n",
	" > hyper-/hypomethylated CGI\n",
	" > differentailly methylated promoter\n",
	" > differentailly methylated exon\n",
	" > differentailly methylated CGI\n",
	"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
))

suppressPackageStartupMessages({
  # doi: 10.1039/C5MB00663E
  library(ReactomePA)
	
	# doi: 10.1089/omi.2011.0118
	# doi: 10.1016/j.xinn.2021.100141
	library(clusterProfiler) 
	
	library(org.Mm.eg.db)
	#library(biomaRt)
		
	library(dplyr)
})

# ~~ Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

go_analysis <- function(gen, treatment, meth, database, region = NULL, organism = "mouse"){
	#gen="F0";treatment=10

	# ~~ background genelist (universe) = genes tested (DMR table) ~~~~~~~~~~~~~
	
	ensembl <- DMR[DMR$gen %in% gen & DMR$treatment %in% treatment, "gene_id"]
	ensembl <- na.omit(unique(ensembl))
	entrez <- as.character(na.omit(ens[ens$gene_id %in% ensembl, "entrez_id"]))
  universe <- list(ensembl = ensembl, entrez = entrez)

	# ~~ genelist ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
  # if CpG-site resolution is used divide between overlapping gene regions
  if (!is.null(region)){
    if (length(region) == 1 && region == "cgi"){
      genes <- DMR[DMR$gen %in% gen & DMR$treatment %in% treatment & DMR$sign & DMR$type %in% meth & !is.na(DMR$cgi), "gene_id"]
    } else {
      genes <- DMR[DMR$gen %in% gen & DMR$treatment %in% treatment & DMR$sign & DMR$type %in% meth & DMR$region %in% region, "gene_id"]
    }
  }

  if (is.null(region)){
    genes <- DMR[DMR$gen %in% gen & DMR$treatment %in% treatment & DMR$sign & DMR$type %in% meth, "gene_id"]
    region <- genomic_feature
  }

	ensembl <- na.omit(unique(genes))
	entrez <- as.character(na.omit(ens[ens$gene_id %in% genes, "entrez_id"]))
	genelist <- list(ensembl = ensembl, entrez = entrez)

	# ~~ enrichment analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	cat(paste0(
		"\nRunning ", database, " analysis for ", gen, "DBP", treatment, " ( ",
		paste(meth, collapse = "/"), " & ", paste(region, collapse = "/")), ") : ",
    "\n"
	)

  if (database == "KEGG"){
    if (organism == "human") org <- "hsa"
    if (organism == "mouse") org <- "mmu"

    go <- clusterProfiler::enrichKEGG(
      gene = genelist$entrez,
      organism = org,
      universe = universe$entrez,
      qvalueCutoff = 0.2,
    )
    go <- data.frame(go)
    
    cat(nrow(go), "hits\n")

    if (nrow(go) == 0) return(NULL)

    go$geneID <- unlist(sapply(go$geneID, simplify = FALSE, \(gene){
      gene <- strsplit(gene, "/")[[1]]

      rows <- match(gene, ens$entrez_id)
      paste(ens$gene_name[rows], collapse = "/")
    }))

    go$ONTOLOGY <- database
  }

  if (database == "GO"){
    if (organism == "human") org <- "org.Hs.eg.db"
    if (organism == "mouse") org <- "org.Mm.eg.db"

    go <- clusterProfiler::enrichGO(
      gene = genelist$ensembl,
      OrgDb = org,
      keyType = "ENSEMBL",
      ont = "ALL",
      pAdjustMethod = "BH",
      universe = universe$ensembl,
      qvalueCutoff = 0.2,
      readable = TRUE
    )
    go <- data.frame(go)

    cat(nrow(go), "hits\n")

    if (nrow(go) == 0) return(NULL)
  }

  if (database == "Reactome"){
    go <- ReactomePA::enrichPathway(
      gene = genelist$entrez,
      organism = organism,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2,
      universe = universe$entrez,
      readable = TRUE
    )
    go <- data.frame(go)

    cat(nrow(go), "hits\n")

    if (nrow(go) == 0) return(NULL)

    go$ONTOLOGY <- database
  }

  # ~~ make go table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  go$gene_name <- sapply(go$geneID, simplify = FALSE, \(gene){
    strsplit(gene, "/")[[1]]
  })
  go$gene_id <- sapply(go$gene_name, simplify = FALSE, \(gene){
    rows <- match(gene, ens$gene_name)
    ens$gene_id[rows]
  })
  go$geneID <- NULL
  go$database <- database

  columns <- c(
    "database",
    "ONTOLOGY",
    "ID",
    "Description",
    "GeneRatio",
    "BgRatio",
    "pvalue",
    "p.adjust",
    "qvalue",
    "Count",
    "gene_name",
    "gene_id"
  )

  go <- go[, columns]
  go <- go[order(-go$Count, go$qvalue), ]
  rownames(go) <- NULL

  return(go)
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

generations <- unique(DMR$gen)
names(generations) <- generations

treatments <- unique(DMR$treatment)
names(treatments) <- paste0("DBP", treatments)

meths <- list("Both" = c("hyper", "hypo"), "Hyper" = "hyper", "Hypo" = "hypo")

# if CpG-site resolution is used divide between overlapping gene regions
if (genomic_feature == "cpg"){
  regions <- list(
    "All" = c("promoter", "exon"),
    "Promoter" = "promoter",
    "Exon" = "exon",
    "CGI" = "cgi"
  )
}

databases <- c("KEGG", "GO", "Reactome")
names(databases) <- databases

# ~~ Code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# make an array:
# gen $ treatment $ meth-type $ genomic-region
# F0	$ DBP10 		$ hyper 		$ promoter 				(example)

if (genomic_feature == "cpg"){
  go <- lapply(generations, \(gen){
    lapply(treatments, \(treatment){
      lapply(meths, \(meth){
        lapply(regions, \(region){
          
          res <- lapply(databases, \(database){
            go_analysis(gen, treatment, meth, database, region)
          })

          Reduce(function(x,y) rbind(x,y), res)

        })
      })
    })
  })
}

if (genomic_feature %in% c("promoter", "islands")){
  go <- lapply(generations, \(gen){
    lapply(treatments, \(treatment){
      lapply(meths, \(meth){
                  
        res <- lapply(databases, \(database){
          go_analysis(gen, treatment, meth, database)
        })

        Reduce(function(x,y) rbind(x,y), res)

      })
    })
  })
}


filename <- paste0("gene_ontology_", genomic_feature, ".Rds")

saveRDS(go, file = filename)

cat(paste(
  "\n~~ gene_ontology.R complete ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
  filename, "generated",
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
))