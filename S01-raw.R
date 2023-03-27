library("Spectra")

## https://rformassspectrometry.github.io/Spectra/

## I get this error when loading “Spectra” > library(“Spectra”)
## Error: package ‘S4Vectors’ required by ‘Spectra’ could not be found

## BiocManager::install("S4Vectors")

## data.frame()
## tibble::tibble()

spd <- DataFrame(msLevel = c(1L, 2L),
                 rtime = c(60, 61))

spd$mz <- list(c(100, 120, 130),
               c(100, 122))
spd$intensity <- list(c(100, 10, 123),
                      c(89, 345))

sp <- Spectra(spd)

sp

spectraVariables(sp)

spectraData(sp)

peaksData(sp)

peaksData(sp)[[1]]

peaksData(sp)[[2]]

sp[1]

sp <- Spectra(mzf)
sp

## Warning message:
## In fun(libname, pkgname) :
##   mzR has been built against a different Rcpp version (1.0.6)
## than is installed on your system (1.0.8.3). This might lead to errors
## when loading mzR. If you encounter such issues, please send a report,
## including the output of sessionInfo() to the Bioc support forum at
## https://support.bioconductor.org/. For details see also
## https://github.com/sneumann/mzR/wiki/mzR-Rcpp-compiler-linker-issue.


## > sp2<- Spectra(mzf)
## Error: BiocParallel errors
##   1 remote errors, element index: 1
##   0 unevaluated and other errors
##   first remote error: Can not open file C:\Users\cmidha\AppData\Local\R\cache\R\rpx\d9889c56f_TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML! Original error was: Error in pwizModule$open(filename): [SpectrumList_mzML::create()] Bad istream.
## In addition: Warning message:
## In serialize(data, node$con) :
##   'package:stats' may not be available when loading

## mzf2 <- msdata::proteomics(full.names = TRUE)[4]
## Spectra(mzf2)

## Accessors - see lecture notes or ?Spectra for more examples
rtime(sp)
tic(sp)
msLevel(sp)
sp$msLevel


## Produce the total ion chromatogram figure.
library(tidyverse)

tmp <- sp |>
    filterMsLevel(1) |>
    spectraData()

plot(tmp$rtime, tmp$totIonCurrent,
     type = "l")

sp |>
    spectraData() |>
    as_tibble() |>
    filter(msLevel == 1) |>
    ggplot(aes(x = rtime,
               y = totIonCurrent)) +
    geom_line()

sp |>
    filterMsLevel(1) |>
    filterRt(c(1000, 3300)) |>
    spectraData() |>
    as_tibble() |>
    ggplot(aes(x = rtime,
               y = totIonCurrent)) +
    geom_line()


msLevel(sp)
centroided(sp)

table(msLevel(sp),
      centroided(sp))


## Backends
## - https://bioconductor.org/packages/release/bioc/html/MsBackendRawFileReader.html
## - https://rformassspectrometry.github.io/MsBackendMgf/articles/MsBackendMgf.html
## - https://github.com/rformassspectrometry/MsBackendTimsTof

plotSpectra(sp[1:12])

sp_i <- sp[2709]

plotSpectra(sp_i)

plotSpectra(sp_i, xlim = c(125, 133))
?plotSpectra

mzLabel <- function(z) {
    z <- peaksData(z)[[1]]
    lbls <- format(z[, "mz"], digits = 5)
    lbls[z[, "intensity"] < 1e6] <- ""
    lbls
}


plotSpectra(sp_i,
            xlim = c(125, 133),
            labels = mzLabel,
            labelPos = 2,
            labelSrt = -30)



sp_ij <- sp[2709:2711]

plotSpectraOverlay(sp_ij,
                   col = c("red", "steelblue",
                           "green"))

plotSpectraMirror(sp[2709], sp[2710])

BiocManager::install("RforMassSpectrometry/PSMatch")

BiocManager::install("RforMassSpectrometry/SpectraVis", force = TRUE)

library(SpectraVis)

SpectraVis::browseSpectra(sp)
SpectraVis::browseSpectra(sp[2708:2718])
SpectraVis::plotlySpectra(sp_i)
