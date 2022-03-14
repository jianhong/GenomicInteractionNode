set.seed(123)
hub_regions <- 
  suppressMessages(createRandomHubs(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                    upstream=500, downstream=500,
                                    maxDist=1e4))
sel_comp_id <- sample(unique(hub_regions$comp_id), 3)
bck_regions <-
  hub_regions[!hub_regions$comp_id %in% sel_comp_id]
bck_regions <- reduce(bck_regions)
sel_regions <- 
  hub_regions[hub_regions$comp_id %in% sel_comp_id]
ol <- findOverlaps(sel_regions, drop.redundant=TRUE, drop.self=TRUE)
sel_regions <- sel_regions[!seq_along(sel_regions) %in%
                             c(queryHits(ol), subjectHits(ol))]
bck_1 <- sample.int(length(bck_regions), size = 300)
bck_2 <- sample.int(length(bck_regions), size = 300)
keep <- bck_1!=bck_2
bck <- Pairs(bck_regions[bck_1[keep]], bck_regions[bck_2[keep]])
sel <- split(sel_regions, sel_regions$comp_id)
sel <- lapply(sel, function(.ele){
  .ele <- reduce(.ele)
  n <- combn(seq_along(.ele), 2)
  Pairs(.ele[n[1, ]], .ele[n[2,]])
})
sel <- Reduce(c, sel)
pr <- c(bck, sel)
test_that("detectHubs works not correct", {
  dh <- detectHubs(pr)
  ol <- findOverlaps(dh$hub_regions, sel_regions)
  expect_true(all(seq_along(sel_regions) %in% subjectHits(ol)))
})
