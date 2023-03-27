library("PSMatch")

idf <- msdata::ident(full.names = TRUE)

id <- PSM(idf)

id

names(id)

nrow(id)
ncol(id)

length(unique(id$spectrumID))

length(unique(id$sequence))

table(id$isDecoy)

data.frame(id)
DataFrame(id)
as_tibble(id)

## Filtering (1)

id2 <- filterPsmDecoy(id)
id2

## Keeping all matches

table(table(id2$spectrumID))

which(table(id2$spectrumID) == 4)

i <- which(id2$spectrumID == "controllerType=0 controllerNumber=1 scan=1774")

DT::datatable(as_tibble(id2[i, ]))

j <- which(id2$spectrumID == "controllerType=0 controllerNumber=1 scan=5490")

DT::datatable(as_tibble(id2[j, ]))

id2r <- reducePSMs(id2)


k <- which(id2r$spectrumID == "controllerType=0 controllerNumber=1 scan=5490")

DT::datatable(as_tibble(id2r[k, ]))


k2 <- which(id2r$spectrumID == "controllerType=0 controllerNumber=1 scan=1774")

DT::datatable(as_tibble(id2r[k2, ]))

## Filtering (cont)

id3 <- filterPsmDecoy(id) |>
    filterPsmRank() |>
    filterPsmShared()

id3 <- filterPSMs(id)

table(table(id3$spectrumID))

id3r <- reducePSMs(id3)

## ------------------------------------------
## combine MS and id data
sp2 <- joinSpectraData(sp, id3r,
                       by.x = "spectrumId",
                       by.y = "spectrumID")



id_filtered <- filterPsmDecoy(id) |>
    filterPsmRank()

id_filtered

## -------------------------------------
## exploring shared and unique peptides

describeProteins(id_filtered)

describePeptides(id_filtered)

am <- makeAdjacencyMatrix(id_filtered)

cc <- ConnectedComponents(am)

connectedComponents(cc, 1)

connectedComponents(cc, 527)

connectedComponents(cc, 38)

connectedComponents(cc, 920)

i <- which(nrows(cc) > 2 & ncols(cc) > 2)

dims(cc)[i, ]

cx <- connectedComponents(cc, 1082)

cx

plotAdjacencyMatrix(cx)


## https://rformassspectrometry.github.io/PSMatch/articles/AdjacencyMatrix.html

## Ex:
## Compare the distribution of raw idenfication scores
## of the decoy and non-decoy hits. Interpret the figure.

library(ggplot2)


msdata::ident(full.names = TRUE) |>
    PSM() |>
    as_tibble() |>
    ggplot(aes(x = MS.GF.RawScore,
               colour = isDecoy)) +
    geom_density()

## -----------------------------
## More visualisation with id

sp2 |>
    spectraData() |>
    as_tibble() |>
    filter(msLevel == 1) |>
    ggplot(aes(x = rtime,
               y = totIonCurrent)) +
    geom_point() +
    geom_line()


## BiocManager::install("RforMassSpectrometry/Spectra")

sp2 <- countIdentifications(sp2)

sp2 |>
    filterMsLevel(1) |>
    spectraData() |>
    as_tibble() |>
    ggplot(aes(x = rtime,
               y = totIonCurrent)) +
    geom_point(aes(
        colour = ifelse(countIdentifications == 0,
                        NA, countIdentifications)),
        alpha = 0.5) +
    geom_line() +
    labs(colour = "Number of ids")


d <- computeMzDeltas(sp)
plotMzDelta(d)

i <- which(sp2$MS.GF.RawScore > 100)[1]

plotSpectra(sp2[i])

addFragments(sp2[i])

plotSpectra(sp2[i],
            labels = addFragments,
            labelPos = 3,
            labelCol = "steelblue")

## --------------------
## Comparing spectra


sample(which(msLevel(sp2) == 2), 2)

compareSpectra(sp2[4557], sp2[4390])

compareSpectra(sp2[862:866])

## 1. Create a new Spectra object containing the MS2 spectra with sequences
##    "SQILQQAGTSVLSQANQVPQTVLSLLR" and "TKGLNVMQNLLTAHPDVQAVFAQNDEMALGALR".

i <- which(sp2$sequence %in% c("SQILQQAGTSVLSQANQVPQTVLSLLR",  "TKGLNVMQNLLTAHPDVQAVFAQNDEMALGALR"))

sp_i <- sp2[i]


## 2. Calculate the 5 by 5 distance matrix between all spectra using
##    compareSpectra. Draw a heatmap of that distance matrix


m <- compareSpectra(sp_i)
m

rownames(m) <-
    colnames(m) <-
    strtrim(sp_i$sequence, 2)

library(pheatmap)
pheatmap(m)

## 3. Compare the spectra with the plotting function seen previously.


filterIntensity(sp_i, 1e4) |>
    plotSpectra(main = sp_i$sequence)

plotSpectraMirror(
    filterIntensity(sp_i[1], 1e4),
    filterIntensity(sp_i[2], 1e4))



plotSpectraOverlay(
    filterIntensity(sp_i[3:5], 1e4),
    col = c("red", "steelblue", "green"))


par(mfrow = c(3, 1))

plotSpectraMirror(
    filterIntensity(sp_i[3], 1e4),
    filterIntensity(sp_i[4], 1e4))

plotSpectraMirror(
    filterIntensity(sp_i[3], 1e4),
    filterIntensity(sp_i[5], 1e4))

plotSpectraMirror(
    filterIntensity(sp_i[4], 1e4),
    filterIntensity(sp_i[5], 1e4))
