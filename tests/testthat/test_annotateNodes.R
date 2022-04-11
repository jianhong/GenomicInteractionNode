set.seed(123)
node_regions <- 
  suppressMessages(createRandomNodes(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                    upstream=500, downstream=500,
                                    maxDist=1e4))
node_regions <- 
  suppressMessages(
     annotateNodes(node_regions,
                   TxDb.Hsapiens.UCSC.hg19.knownGene,
                   org.Hs.eg.db))
enr <- enrichmentAnalysis(node_regions, org.Hs.eg.db, onto="BP")
test_that("annotateNodes works not correct", {
    gene <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
    p <- promoters(gene, upstream = 2000, downstream = 500)
    idx <- vapply(node_regions$gene_id, FUN = function(.ele){
      !is.na(.ele[1])
    }, FUN.VALUE=logical(1))
    ol <- findOverlaps(node_regions, p)
    expect_true(all(queryHits(ol) %in% which(idx)))
    expect_true(all(which(idx) %in% queryHits(ol)))
})

test_that("enrichmentAnalysis works not correct", {
  x <- enr[[1]][[1]]
  pvalue <- phyper(q=x$countInDataset-1,
                   m=x$countInGenome,
                   n=x$totalGeneInGenome-x$countInGenome,
                   k=x$totalGeneInDataset,
                   lower.tail = FALSE, log.p = FALSE)
  expect_true(all(pvalue==x$pvalue))
})