#' Annotate hub regions
#' @description Assigne gene id and gene symbols to hub regions by interacted.
#' @param hub_regions GRanges object represent regions interacted with hubs. 
#' @param txdb An object of \link[GenomicFeatures:TxDb-class]{TxDb} to extract
#' gene information
#' @param orgDb An object of \link[AnnotationDbi:AnnotationDb-class]{OrgDb}
#' to extract gene symbols
#' @param upstream,downstream An integer(1) value indicating the number of bases
#' upstream or downstream from the transcription start site. For additional
#' details see \link[GenomicFeatures:transcripts]{promoters}.
#' @param ... parameter can be passed to
#' \link[GenomicFeatures:transcripts]{genes}
#' @return GRanges object with gene_id and symbols metadata.
#' @export
#' @importClassesFrom GenomicFeatures TxDb
#' @importClassesFrom AnnotationDbi OrgDb
#' @importMethodsFrom GenomicFeatures genes promoters
#' @importMethodsFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits SimpleList
#' @importFrom methods is
#' @examples 
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene) ## for human hg19
#' library(org.Hs.eg.db) ## used to convert gene_id to gene_symbol
#' set.seed(123)
#' hub_regions <- createRandomHubs(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' annoHubs(hub_regions, TxDb.Hsapiens.UCSC.hg19.knownGene, org.Hs.eg.db)

annoHubs <- function(hub_regions,
                     txdb, orgDb,
                     upstream=2000,
                     downstream=500, ...){
  stopifnot(is(txdb, "TxDb"))
  stopifnot(is(orgDb, "OrgDb"))
  stopifnot(is(upstream, "numeric"))
  stopifnot(is(downstream, "numeric"))
  ### annotation for promoter region only
  gene <- genes(txdb, ...)
  promoter <- promoters(gene,
                        upstream=upstream,
                        downstream = downstream)
  ol <- findOverlaps(hub_regions, promoter)
  ol_gene_id <- split(promoter[subjectHits(ol)]$gene_id, queryHits(ol))
  ol_gene_id <- SimpleList(ol_gene_id)
  hub_regions$gene_id <- SimpleList(NA)
  hub_regions[as.numeric(names(ol_gene_id))]$gene_id <- ol_gene_id
  ol_gene_symbol <- lapply(ol_gene_id, FUN=eg2symbol, orgDb=orgDb)
  hub_regions$symbols <- SimpleList(NA)
  hub_regions[as.numeric(names(ol_gene_id))]$symbols <- 
    SimpleList(ol_gene_symbol)
  hub_regions
}

## help function: add gene symbols
eg2symbol <- function(id, orgDb){
  stopifnot(is(orgDb, "OrgDb"))
  env <- get(sub(".db", "SYMBOL", orgDb$packageName))
  unlist(AnnotationDbi::mget(id, ## overlaod the mget from AnnotationDbi
                             envir = env,
                             ifnotfound = NA))
}
