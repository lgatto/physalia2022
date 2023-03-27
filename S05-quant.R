## -------------------------------
## Suppl ex

PXD022816 <- PXDataset("PXD022816")

(mzids <- pxget(PXD022816, grep("mzID", pxfiles(PXD022816))[1:3]))
(mzmls <- pxget(PXD022816, grep("mzML", pxfiles(PXD022816))[1:3]))

sp <- Spectra(mzmls)
sp

id <- PSM(mzids)
id_filtered <- filterPSMs(id)
id_filtered

pkey <- paste0(sub("^.+_QEP", "", basename(dataOrigin(sp))),
               "::", sub("^.+scan=", "", sp$spectrumId))

pkey <- paste0(sub("^.+QEP", "QEP", basename(dataOrigin(sp))),
               sp$spectrumId)


head(paste0(sub("P.+QEP", "QEP", id$spectrumFile),
            id$spectrumID))

## is.vector(id_filtered) ## <<<<<<<<<<
## [1] FALSE
## > id_filtered$spectrumFile
## CharacterList of length 21748
## [["controllerType=0 controllerNumber=1 scan=10000"]] QEP2LC6_HeLa_50ng_251120_01-calib.mzML ...
## [["controllerType=0 controllerNumber=1 scan=10001"]] QEP2LC6_HeLa_50ng_251120_01-calib.mzML ...
## [["controllerType=0 controllerNumber=1 scan=10002"]]


library(QFeatures)
data(feat1)
se <- feat1[[1]]

## ---------------------------------
## The SummarizedExperiment class
se
dim(se)
assay(se)
rowData(se)
colData(se)

colData(se)$condition <- c("CTRL", "COND1")
colData(se)$sex <- c("M", "F")

colData(se)

se

colData(se[1:4, 1])

## ---------------------------------
## The QFeatures class
feat1

feat1[[1]]

rowData(feat1[[1]])
assay(feat1[[1]])
colData(feat1[[1]])

colData(feat1)
colData(feat1)$sex <- c("M", "F")
colData(feat1)

feat1[["psms"]]

## feature aggregation

feat1 <- aggregateFeatures(feat1,
                           i = "psms",
                           fcol = "Sequence",
                           name = "peptides",
                           fun = colMeans)

feat1

assay(feat1[[2]])
rowData(feat1[[2]])

## Ex: aggregate peptides into proteins using median

feat1 <- aggregateFeatures(feat1,
                           i = "peptides",
                           fcol = "Protein",
                           name = "proteins",
                           fun = colMedians)

feat1

assay(feat1[[3]])
rowData(feat1[[3]])

colData(feat1)

## Discussion: what is the default way to
##             aggregate quantitation data?

## input: matrix
## output: vector of length ncol(matrix)
topN <- function(x) {
    apply(x, 2,
          function(xx) sum(head(sort(xx, decreasing = TRUE), 3)))
}

m <- matrix(1:15, ncol = 3)
topN(m)

## Filtering and subsetting

feat1[, , 1]

feat1[[1]]

feat1[, 2, ]

feat1["ProtA", ]

filterFeatures(feat1, ~ pval < 0.05)

filterFeatures(feat1, ~ location == "Mitochondrion")

data(hlpsms)

readQFeatures(hlpsms,
              ecol = 1:10,
              name = "psms")


se <- readSummarizedExperiment(hlpsms,
                               ecol = 1:10)

ft <- QFeatures(list(psms = se))


f <- msdata::quant(full.names = TRUE)
f

## Ex: (1) Read these data in as either a SummerizedExperiment or a
##     QFeatures object and (2) annotated the experiment.

names(read.delim(f))

e <- grep("Intensity\\.", names(read.delim(f)))

cptac_se <- readSummarizedExperiment(f, ecol = e, sep = "\t")

colnames(cptac_se) <- sub("Intensity\\.", "", colnames(cptac_se))

colData(cptac_se)$condition <- rep(LETTERS[1:2], each = 3)
colData(cptac_se)$id <- rep(7:9, 2)

keep_var <- c("Sequence", "Proteins", "Leading.razor.protein", "PEP",
              "Score", "Reverse", "Potential.contaminant")
rowData(cptac_se) <- rowData(cptac_se)[, keep_var]


## psms <- readSummarizedExperiment(fpsms)
## peps <- readSummarizedExperiment(fpeps)
## prts <- readSummarizedExperiment(fprts)
##
## qf <- QFeatures(list(psms = psms, peptides = peps,
##                      proteins = prts))



## ------------------------------------------
## PSMatch/ConnectedComponents on quant data

prot <- rowData(cptac_se)$Proteins
names(prot) <- rowData(cptac_se)$Sequence

am <- makeAdjacencyMatrix(prot, split = ";")

am[1:10, 1:5]

cc <- ConnectedComponents(am)

(cctab <- prioritiseConnectedComponents(cc))

## BiocManager::install("factoextra")
library(factoextra)

fviz_pca(prcomp(cctab, scale = TRUE, center = TRUE))

plotAdjacencyMatrix(connectedComponents(cc, 1081))

rownames(connectedComponents(cc, 1081))

plotAdjacencyMatrix(connectedComponents(cc, 151))

pdf("cc223.pdf") ## png()
plotAdjacencyMatrix(connectedComponents(cc, 223))
dev.off()

plotAdjacencyMatrix(connectedComponents(cc, 1200))

colnames(connectedComponents(cc, 1200))

peps_to_remove <- c("GVLLYGPPGTGK", "IIMATNR")

which(rowData(cptac_se)$Sequence %in% peps_to_remove)

## --------------------------------
## Analysis pipeline

## - Exploratory data analysis (PCA)
## - Data cleaning
## - Transformation and normalisation
## - Aggregation
## - Downstream analysis: limma, MSqRob/msqrob2, ProDA

library("tidyverse")
library("ggplot2")
library("QFeatures")
library("limma")
library("factoextra")

## - Missing data management (filtering and/or imputation)

cptac_se

cptac_se <- zeroIsNA(cptac_se)

barplot(nNA(cptac_se)$nNAcols$pNA)

table(nNA(cptac_se)$nNArows$nNA)

cptac_se <- filterNA(cptac_se, pNA = 4/6)

## filterNA(cptac_se)

nNA(cptac_se)

impute(cptac_se, method = "knn")

## Lieven Clement - MSqRob

table(rowData(cptac_se)$Reverse)

table(rowData(cptac_se)$Potential.contaminant)

rowData(cptac_se) |>
    as_tibble() |>
    ggplot(aes(x = PEP,
               colour = Reverse)) +
    geom_density()

rowData(cptac_se) |>
    as_tibble() |>
    ggplot(aes(x = Score,
               colour = Reverse)) +
    geom_density()

cptac <- QFeatures(list(peptides = cptac_se))
colData(cptac) <- colData(cptac_se)

cptac <- cptac |>
    filterFeatures(~ Reverse != "+") |>
    filterFeatures(~ Potential.contaminant != "+") |>
    filterFeatures(~ PEP < 0.05)

cptac <- logTransform(cptac,
                      i = "peptides",
                      name = "log_peptides")

cptac <- normalize(cptac,
          i = "log_peptides",
          name = "lognorm_peptides",
          method = "center.median")


limma::plotDensities(assay(cptac[[1]]))

par(mfrow = c(1, 2))
limma::plotDensities(assay(cptac[[2]]))
limma::plotDensities(assay(cptac[[3]]))

boxplot(assay(cptac[[2]]))


cptac <- aggregateFeatures(cptac,
                           i = "lognorm_peptides",
                           name = "proteins_med",
                           fcol = "Leading.razor.protein",
                           fun = colMedians,
                           na.rm = TRUE)


pca_pep <- cptac[["lognorm_peptides"]] %>%
    filterNA() %>%
    assay() %>%
    t() %>%
    prcomp(scale = TRUE, center = TRUE) %>%
    fviz_pca_ind(habillage = cptac$condition, title = "Peptides")


pca_prot <-
    cptac[["proteins_med"]] %>%
    filterNA() %>%
    assay() %>%
    t() %>%
    prcomp() %>%
    fviz_pca_ind(habillage = cptac$condition,
                 title = "Proteins (median aggregation)")

library(patchwork)

pca_pep + pca_prot

cptac["P02787ups|TRFE_HUMAN_UPS", , c("lognorm_peptides", "proteins_med")]


longFormat(cptac["P02787ups|TRFE_HUMAN_UPS", ,
                 c("lognorm_peptides", "proteins_med")]) %>%
    as_tibble() %>%
    mutate(condition = ifelse(grepl("A", colname), "A", "B")) %>%
    ggplot(aes(x = colname, y = value, colour = rowname, shape = condition)) +
    geom_point(size = 3) +
    geom_line(aes(group = rowname)) +
    facet_grid(~ assay) +
    ggtitle("P02787ups|TRFE_HUMAN_UPS")

## %>% (magrittr) ~= |> (R since 4.0)
