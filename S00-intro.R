## --------------------------------
## Lecture notes

## https://rformassspectrometry.github.io/docs/

## --------------------------------
## Installation

## CRAN - ggplot, dplyr
## Bioconductor
## Github

## install.packages("BiocManager")

BiocManager::install("remotes")

BiocManager::install("ggplot2") ## CRAN
BiocManager::install("limma")   ## Bioconductor
BiocManager::install("ProtGenerics") ## Bioconductor (release)
BiocManager::install("RforMassSpectrometry/ProtGenerics") ## Github (devel)
BiocManager::install("RforMassSpectrometry/PSMatch")      ## Github


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("msdata")
BiocManager::install("mzR")
BiocManager::install("rhdf5")
BiocManager::install("rpx")
BiocManager::install("MsCoreUtils")
BiocManager::install("QFeatures")
BiocManager::install("Spectra")
BiocManager::install("RforMassSpectrometry/ProtGenerics")
BiocManager::install("RforMassSpectrometry/PSMatch")
BiocManager::install("RforMassSpectrometry/SpectraVis")
BiocManager::install("pheatmap")

## upon applying this command I am getting this message

## The downloaded source packages are in
##      ‘C:\Users\sa1361\AppData\Local\Temp\RtmpCQWkam\downloaded_packages’
## Installation paths not writeable, unable to update packages
##   path: C:/Program Files/R/R-4.1.3/library
##   packages:
##     MASS, Matrix, nlme, survival
## Old packages: 'fansi'
## Update all/some/none? [a/s/n]:

packageVersion("ProtGenerics")

## Bioconductor:
## ProtGenerics' 1.26.0 - release
##              ‘1.27.2’- devel

library("package")

remove.packages("msdata", "/opt/Rpatched/lib64/R/library")
library("msdata")
library("msdata", lib.loc = "/home/lgatto/R/x86_64-pc-linux-gnu-library/4.1")



## --------------------------------
## Getting data from PRIDE

library("rpx")

px <- PXDataset("PXD000001")
px


## Querying ProteomeXchange for PXD000001.
## Error in file(file, “rt”) :
##   cannot open the connection to ‘ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2012/03/PXD000001/README.txt’
## In addition: Warning message:
## In file(file, “rt”) :
##   cannot open URL ‘ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2012/03/PXD000001/README.txt’: FTP status was ‘550 Requested action not taken; file unavailable’

packageVersion("rpx")

BiocManager::install("lgatto/rpx")
unloadNamespace("rpx")

library("rpx")

px

pxfiles(px)
pxtax(px)
pxref(px)
pxprotocols(px)

## file name
fafile <- pxget(px,"erwinia_carotovora.fasta")
## file index
fafile <- pxget(px, 1)
## choose from menu
pxget(px)
## pxget(px, "all")


mzf <- pxget(px, 8)

mzf

## custom cache in my current project directory
my_cache <- BiocFileCache::BiocFileCache("/home/lgatto/wd/physalia")
px2 <- PXDataset("PXD000001", cache = my_cache)
pxfiles(px2)

## --------------------------------
## Getting data from msdata

## BiocManager::install("msdata")
library("msdata")

proteomics()

proteomics(full.names = TRUE)

## mzf <- proteomics(full.names = TRUE)[4]

ident()

quant()
