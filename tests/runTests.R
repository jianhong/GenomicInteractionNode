require("GenomicInteractionNodes") ||
  stop("unable to load Package:GenomicInteractionNodes")
require("TxDb.Hsapiens.UCSC.hg19.knownGene") || 
  stop("unable to load Package:TxDb.Hsapiens.UCSC.hg19.knownGene")
require("org.Hs.eg.db") || stop("unable to load Package:org.Hs.eg.db")
require("GO.db") || stop("unable to load Package:GO.db")
require("testthat") || stop("unable to load testthat")
test_check("GenomicInteractionNodes")