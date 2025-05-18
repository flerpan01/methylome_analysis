# README

This repository holds a nextflow pipeline for analysing bulk RRBS studies. The experimental design is build on multiple generations (F0, F1, F2, etc) and doses (0, 10, 100, 1000, etc). The pipeline expect coverage files (`sample.cov.gz`) as input generated with Bismark and outputs **(1)** differentially methylated regions, **(2)** differentially methylated genes and **(3)** gene ontology analysis results, all in `Rds` file format containing all generations and doses in a single file, respectively. The pipeline also outputs the tables as `excel` files to be included as supplementary tables in a scientic report.

<details>
  <summary>Coverage files</summary>

To generate methylation coverage files from sequencing files refer to [nf-core/methylseq pipeline](https://nf-co.re/methylseq/latest/)

</details>

<details>
  <summary>Reference genome</summary>

This pipeline is per default setup for the [mouse genome (GRCm39)](https://www.ensembl.org/Mus_musculus/)

</details>

<details>
  <summary>Statistics</summary>

Differentially methylated regions are identified with the R-package [`methylKit`](https://bioconductor.org/packages/release/bioc/html/methylKit.html), using logistic regression test and with overdispertion correction and calulating the generic mean methylation between groups. `FDR < 0.01` was used for multiple testing correction (Benjamini-Hochberg qvalue).

</details> 

<details>
  <summary>Gene ontology analysis</summary>

To investigate if any biological functions, processes or pathways are enriched (over-represented) the _Over Representation Analysis (ORA)_ [Boyle et al., 2004](https://doi.org/10.1093/bioinformatics/bth456) method is used. ORA uses hypergeometric distribution and compares the differentially methylated genes with all genes in the dataset. The _p_-values are adjusted to _q_-values for multiple corretion (significance threshold `qvalue < 0.2`).

Enrichment is analysed in three databases; (1) Gene Ontology (**GO**), (2) Kyoto Encyclopedia of Genes and Genomes (**KEGG**), and **Reactome** pathways. GO and KEGG enrichment are tested with the R-package [`clusterProfiler`](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html), [Yu et al., 2012](https://doi.org/10.1089/omi.2011.0118), [Wu et al., 2021](https://doi.org/10.1016/j.xinn.2021.100141). The reactome pathways are tested with the R-package [`ReactomePA`](https://bioconductor.org/packages/release/bioc/html/ReactomePA.html), [Yu et al., 2016](https://doi.org/10.1039/C5MB00663E). 

</details>

## The pipeline

The nextflow pipeline produce the following:

+ PCA plots (for sanity check of experiment)
+ Betavalues matrix (for plotting methylation levels, in beta-values, over a tested region)
+ methylKit-objects
+ Differentially methylated regions table
+ Differentially methylated genes table
+ Gene ontology analysis
+ Supplementary files (plots, excel-tables)
+ Concatinated tables (for easy import for results report)

## Preparation

### Metadata

Setup the `data/metadata.csv` to look like the following:

```csv
gen,id,treatment,...
F0,F0_1,0,...
F0,F0_2,0,...
F0,F0_3,0,...
F0,F0_4,0,...
```

Each row represents a sample in the column order: generation, sample id and treatment/dose.

### Gene annotations

The pipeline expects three gene annotations files (1) CpG-islands positions (`data/cpgislands_GRCm39.bed`), (2) curated reference sequencing database (`data/refseq_UCSC_GRCm39.bed`) and ensembl gene annotations (`data/ensembl-dataset.csv.gz`)

## Parameters

The pipeline accepts seven parameters:

Experimental design:

+ generations (F0, F1, F2, etc)
+ doses (0, 10, 100, 1000, etc)
+ genomic features (CpG-sites, Promoters, CpG-islands)

`methylKit` options (defaults in parenthesis):

+ differential methylation test (`Chisq`)
+ overdispersion (`MN`, overdispersion accounted)
+ effect (`mean`)
+ multiple testing correction (`BH`, Benjamini-Hochberg, FDR)

## Containers

For reproducibility this pipeline uses two singularity containers, which can be downloaded from the [Cloud Library](https://cloud.sylabs.io/library). The `bulk-rrbs` holds most of the R-packages used in the analysis, while `gene-ontology` holds gene ontology related R-packages

```sh
# apptainer (instead of singulartiy) also works

IMAGE1='library://flerpan01/singularity-r/bulk-rrbs:1.1'
IMAGE2='library://flerpan01/singularity-r/gene-ontology:1.0'

singularity pull ${IMAGE1}
singularity pull ${IMAGE2}
```

<details>
  <summary>Run interactively</summary>

To run scripts manually with the containers use the `exec` flag or run the script interactively with `shell`.

```sh
# execute script
singularity exec ${IMAGE} <scriptfile>

# run script interactively
singularity shell ${IMAGE}
$ Rscript <scriptfile>
```

</details>