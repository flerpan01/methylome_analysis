#!/usr/bin/env Rscript

arg <- commandArgs(trailingOnly = TRUE)
# arg1 = genomic_feature:       islands, promoter, cpg
genomic_feature=arg[1]

cat(paste(
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
  "Genomic feature:\t", genomic_feature, "\n",
  " > Make plots of the DMRs, beta-values on y-axis\n",
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
))

suppressPackageStartupMessages({
  library(methylKit)
  library(ggplot2)
  library(ggsignif)
  library(cowplot)
  library(patchwork)
  library(RColorBrewer)
  library(scales)
  library(ggrepel)
  library(gtools)
})

# ~~ Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

meth_plotter <- function(dmr, size = 12, labeltext = TRUE){
  # set colors
  colors <- brewer.pal(8, "Set1")[c(3,5,8)]
  names(colors) <- c("0", "10", "100")

  data <- lapply(names(methMat), \(gen){
    dmrs <- rownames(methMat[[gen]])

    if (dmr %in% dmrs){
      data <- data.frame(betavalues = methMat[[gen]][dmr, ])
      data$id <- rownames(data)
      data$treatment <- meta$treatment[match(data$id, meta$id)]
    } else {
      data <- data.frame(
        betavalues = NA,
        id = NA,
        treatment = c(0, 10, 100)
      )
    }

    data$gen <- gen
    data$treatment <- factor(data$treatment)

    return(data)
  })
  names(data) <- names(methMat)
  
  # extract pvalue annotations
  dat <- subset(meth, dmr_id %in% dmr & sign %in% TRUE)
  dat$sign <- sapply(dat$treatment, \(x){
    gsub("treatment", "", x)
  })
  dat$annotation <- sapply(dat$qvalue, \(x){
    if (x < 0.05) out <- "*"
    if (x < 0.01) out <- "**"
    if (x < 0.001) out <- "***"
    return(out)
  })

  baseplot <- function(data, dat, gen, size, labeltext = TRUE){
    anno <- dat[dat$gen %in% gen, ]
    ymax <- ifelse(nrow(anno) > 1, 120, 110)
    
    expression <- nrow(data[[gen]]) > 3

    if (!expression){
      p <- ggplot(data[[gen]], aes(x = factor(treatment))) + 
        theme_linedraw(size) + background_grid() + 
        labs(title = gen, y = NULL, x = NULL)
    } else {
      p <- ggplot(data[[gen]], aes(treatment, betavalues, fill = treatment)) +
        scale_fill_manual(values = colors) +
        labs(x = NULL, y = "Betavalue", title = gen, fill = "[DBP] mg/kg") +
        theme_linedraw(size) + background_grid() +
        scale_y_continuous(
          limits = c(0, ymax),
          breaks = seq(0, 100, by = 20),
          labels = seq(0, 100, by = 20)) +
        theme(plot.title = element_text(hjust = 0.5)) +

        { if (gen == "F0") theme(legend.position = "none") } +
        { if (gen == "F1") theme(legend.position = "none", axis.title.y = element_blank()) } +
        { if (gen == "F2") theme(legend.position = "right", axis.title.y = element_blank()) } +

        { if (nrow(anno) > 0)
          geom_signif(
            comparisons = lapply(anno$treatment, \(x) c(0, x)),
            step_increase = 0.1,
            annotations = anno$annotation, map_signif_level = TRUE,
            textsize = 5) } +

        stat_boxplot(geom = "errorbar", width = 0.2) +
        geom_boxplot(outlier.shape = NA, width = 0.6, coef = 0) +

        { if (labeltext)
          geom_label_repel(
            aes(label = id),
            position = position_jitter(width = 0.2, seed = 1, height = 0),
            point.padding = 0.1,
            size = 3,
            show.legend = FALSE,
            min.segment.length = 0,
            max.overlaps = 2) } +

        geom_jitter(
          shape = 21,
          position = position_jitter(width = 0.2, seed = 1, height = 0),
          size = 2.5,
          show.legend = NULL)
    }
    
    return(p)
  }

  pltlist <- lapply(names(data), \(gen) baseplot(data, dat, gen, size))
  p <- plot_grid(plotlist = pltlist, nrow = 1, rel_widths = c(1, 1, 1.5))

  readableDMRid <- function(dmr){
    chrom <- strsplit(dmr, ":")[[1]][1]
    pos <- as.numeric(strsplit(dmr, ":")[[1]][2])
    paste0(chrom, ":", comma(pos))
  }

  # Add genomic info
  title <- readableDMRid(dmr)

  # add gene info
  geneID <- ifelse(dat$gene_name[1] == "", dat$gene_id[1], dat$gene_name[1])
  subtitle <- paste0("DMR location: ", dat$region[1])
  caption <- NULL
  if (!is.na(dat$cgi[1])) subtitle <- paste0(subtitle, " & CpG-island")
  if (!is.na(dat$gene_id[1])) caption <- paste0(geneID, " (", dat$gene_type[1], ") ", dat$gene_info[1])

  p <- p + 
    plot_annotation(
      title = title,
      subtitle = subtitle,
      caption = caption
    ) &

    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      plot.caption = element_text(size = 12))

  return(p)
}

# ~~ Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# build list of betavalues matrix
files <- list.files(
  path = ".",
  pattern = paste0("betavalues_", genomic_feature),
  recursive = TRUE,
  full.names = TRUE
)
methMat <- lapply(files, readRDS)
gen <- sapply(methMat, \(x){
  unique(sub("(F[0-9]).*", "\\1", colnames(x)))
})
names(methMat) <- gen

# metadata
meta <- read.csv("metadata.csv")

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
#dmrs <- meth[order(meth$qvalue), "dmr_id"][1:100]
dmrs <- unique(meth$dmr_id[meth$sign])
dmrs <- mixedsort(dmrs)

# Plot the DMRs
filename_geneplots <- paste0("DMR_", genomic_feature, ".pdf")
pdf(file = filename_geneplots, width = 10, height = 4)

for ( dmr in dmrs ){
  show(dmr)
  print(meth_plotter(dmr))
}

dev.off()

cat(paste(
  "\n\n~~ supplementary.R complete! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n",
  "DMRs plotted:\t\t", comma(length(dmrs)), "\n",
  "Generated files:\t", filename_geneplots,
  "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
))