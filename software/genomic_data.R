# Installation for genomic analysis in R

# The packages below include tools from Bioconductor, a software suite designed for genomic analsyis in R. These packages include (relatively) large data files, as well as many dependencies (other packages that are installed with the main package). Please plan ahead for installing these tools, and be aware that the output from these installs is verbose and may include warnings. Feel free to show up early to lecture if you'd like assistance confirming the installs completed satisfactorily.

# install data.table from CRAN
install.packages("data.table")

# install package to manage Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# install bioconductor packages
BiocManager::install(c("VariantAnnotation", "GenomicRanges", "Rsamtools", "BSgenome.Hsapiens.UCSC.hg19", "GenomicFeatures", "TxDb.Hsapiens.UCSC.hg19.knownGene", "biomaRt"), update = TRUE, ask = FALSE)

# if you receive the following prompt in the console, hit Enter (confirming yes)
# Do you want to install from sources the package which needs compilation? (Yes/no/cancel) 
