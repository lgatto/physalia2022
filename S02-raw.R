BiocManager::version()
BiocManager::valid()


## convert <- paste0("mono ~/bin/ThermoRawFileParser/1.4.0/bin/x64/Debug/ThermoRawFileParser.exe",
##                   " -i=", raws,
##                   " -b=", sub("raw", "mzML", basename(raws)),
##                   " -f=2", ## indexed mzML
##                   " -p=1") ## no peak picking for MS level 1
##
## sapply(convert, system)

sp

plotlySpectra(sp[2807])

plotlySpectra(pickPeaks(sp[2807]))


plotSpectra(sp[2807], xlim = c(494.3, 494.5), type = "h")

plotSpectraMirror(sp[2807],
                  filterIntensity(sp, c(1e7, Inf))[2807],
                  xlim = c(490, 500))

plotSpectraMirror(pickPeaks(sp[2807]),
                  pickPeaks(filterIntensity(sp, c(1e7, Inf))[2807]),
                  xlim = c(490, 500))

abline(h = 1e8, lty = "dotted")
abline(h = -1e8, lty = "dotted")

## -------------------------------
## Visualisation exercise

## Chromatogram
with(spectraData(filterMsLevel(sp, 1)),
     plot(rtime, totIonCurrent, type = "l"))
abline(v = rtime(sp)[2807], col = "red")

## Can MS1 and 10 MS2 group
ms_2 <- filterPrecursorScan(sp, 2807)
ms_2

## MS1 with precursor peaks
plotSpectra(sp[2807], xlim = c(400, 1000))
abline(v = precursorMz(ms_2)[-1], col = "grey")
abline(v = precursorMz(ms_2)[2], col = "red")

## Zoom in
plotSpectra(sp[2807], xlim = c(521.2, 522.5), type = "l")
abline(v = precursorMz(ms_2)[2], col = "red")

## Zoom in and compare peak picking
plotSpectraMirror(sp[2807], pickPeaks(sp[2807]),
                  xlim = c(521.2, 522.5))

## MS2 scans
plotSpectra(ms_2[-1])
