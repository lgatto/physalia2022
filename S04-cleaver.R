download.file("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz",
              "UP000005640_9606.fasta.gz")

myp <- c("MDEGVNFKPTEAGR", "NGIPFTVYDNVR",
         "DYALNTDSAAGLLIR", "IDITDAETLSR",
         "SSCSSKPNLDTMCK", "IDITDAETLLLLR")

library(Biostrings)
library(cleaver)

pr <- readAAStringSet("UP000005640_9606.fasta.gz")

pe <- cleave(pr)

pe <- unique(unname(as.character(unlist(pe))))

myp %in% pe

setdiff(pe, myp)
