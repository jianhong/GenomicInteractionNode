#' Gene ontology enrichment analysis
#' @description GO enrichment analysis for nodes
#' @param node_regions GRanges object represent regions interacted with nodes.
#' The object must be annotated by \link{annotateNodes} 
#' with comp_id and gene_id in the metadata. 
#' @param orgDb An object of \link[AnnotationDbi:AnnotationDb-class]{OrgDb}
#' to extract gene symbols.
#' @param onto Ontology category.
#' @param evidence The acceptable evidence code.
#' @param minGeneNum An integer(1) value indicating the minimal number of gene
#' to start the enrichment analysis. If total gene counts is smaller than
#' the `minGeneNum`, the NULL will be returned.
#' @param ... Not used.
#' @return A list with element enriched and enriched_in_compound.
#' Or NULL if total counts of gene is smaller than `minGeneNum`.
#' @export
#' @importClassesFrom AnnotationDbi OrgDb
#' @importClassesFrom GenomicRanges GRanges
#' @importMethodsFrom AnnotationDbi mappedkeys select
#' @importFrom stats phyper
#' @importFrom stats p.adjust
#' @import GO.db
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene) ## for human hg19
#' library(org.Hs.eg.db) ## used to convert gene_id to gene_symbol 
#' library(GO.db)
#' set.seed(123)
#' node_regions <- createRandomNodes(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' node_regions <- 
#'     annotateNodes(node_regions,
#'                   TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                   org.Hs.eg.db)
#' enr <- enrichmentAnalysis(node_regions, org.Hs.eg.db, onto="BP")

enrichmentAnalysis <- function(node_regions, orgDb, onto=c("BP", "CC", "MF"),
                      minGeneNum=3,
                      evidence=list("Experimental_evidence_codes"=
                                      c("EXP", "IDA", "IPI", "IMP", "IGI",
                                        "IEP", "HTP", "HDA", "HMP", "HGI",
                                        "HEP"),
                                    "Phylogenetically-inferred_annotations"=
                                      c("IBA", "IBD", "IKR", "IRD"),
                                    "Computational_analysis_evidence_codes"=
                                      c("ISS", "ISO", "ISA", "ISM", "IGC",
                                        "RCA"),
                                    "Author_statement_evidence_codes"=
                                      c("TAS", "NAS"),
                                    "Curator_statement_evidence_codes"=
                                      c("IC", "ND"),
                                    "Electronic_annotation_evidence_code"=
                                      c("IEA")),
                      ...){
  stopifnot(is(orgDb, "OrgDb"))
  check_node_region(node_regions)
  onto <- match.arg(onto, choices = c("BP", "CC", "MF"), several.ok = TRUE)
  
  all_gene_id <- unlist(node_regions$gene_id)
  all_gene_id <- unique(all_gene_id[!is.na(all_gene_id)])
  if(length(all_gene_id)<minGeneNum){
    return(NULL)
  }
  ti <- termInfo(orgDb, onto=onto, evidence=evidence)
  ##term info(ti) is a list, named as onto, with dataframe
  ### for all nodes
  enrich_all <- lapply(ti, hyperGT, gene_id=all_gene_id, orgDb=orgDb)
  ### for each node compound
  comp_gene_id <- split(node_regions$gene_id, node_regions$comp_id)
  comp_gene_id <- lapply(comp_gene_id, function(.ele){
    .ele <- unique(unlist(.ele))
    .ele[!is.na(.ele)]
  })
  comp_gene_id <- comp_gene_id[lengths(comp_gene_id)>=minGeneNum]
  enrich_comp <- lapply(ti, function(df_all){
    lapply(comp_gene_id, hyperGT, df_all=df_all, orgDb=orgDb)
  })
  list(enriched=enrich_all,
       enriched_in_compound=enrich_comp)
}

## enrichment analysis help functions
check_node_region <- function(node_regions, check_col=c("comp_id", "gene_id")){
  stopifnot(is(node_regions, "GRanges"))
  for(j in check_col){
    if(length(node_regions)!=length(mcols(node_regions)[, j])){
      stop("node_regions does not contain metadata",
           j)
    }
  }
}
## add ancestors
addAncestors <- function(go_eg, onto){
  onto <- match.arg(onto, choices = c("BP", "CC", "MF"), several.ok = FALSE)
  GOIDs <- go_eg$GO
  if(length(GOIDs)<0){
    return(go_eg)
  }
  GOIDs <- unique(GOIDs)
  env <- get(paste0("GO", onto, "ANCESTOR"), envir = GO.db)
  Ancestors <- AnnotationDbi::mget(GOIDs, envir = env, ifnotfound = NA)
  Ancestors <- data.frame(GO=rep(names(Ancestors), lengths(Ancestors)),
                          ANCESTORS=unlist(Ancestors))
  Ancestors <- Ancestors[!is.na(Ancestors$ANCESTORS), , drop=FALSE]
  Ancestors <- Ancestors[Ancestors$ANCESTORS!="all", , drop=FALSE]
  Ancestors <- merge(Ancestors, go_eg, by="GO")
  Ancestors <- Ancestors[, c("ANCESTORS", "ENTREZID")]
  colnames(Ancestors) <- colnames(go_eg)
  rownames(Ancestors) <- NULL
  go_eg <- rbind(go_eg, Ancestors)
  go_eg <- unique(go_eg)
  go_eg[order(go_eg$ENTREZID), , drop=FALSE]
}
## get term information
# Evidence code:
# Experimental evidence codes
# Inferred from Experiment (EXP)
# Inferred from Direct Assay (IDA)
# Inferred from Physical Interaction (IPI)
# Inferred from Mutant Phenotype (IMP)
# Inferred from Genetic Interaction (IGI)
# Inferred from Expression Pattern (IEP)
# Inferred from High Throughput Experiment (HTP)
# Inferred from High Throughput Direct Assay (HDA)
# Inferred from High Throughput Mutant Phenotype (HMP)
# Inferred from High Throughput Genetic Interaction (HGI)
# Inferred from High Throughput Expression Pattern (HEP)
# Phylogenetically-inferred annotations
# Inferred from Biological aspect of Ancestor (IBA)
# Inferred from Biological aspect of Descendant (IBD)
# Inferred from Key Residues (IKR)
# Inferred from Rapid Divergence (IRD)
# Computational analysis evidence codes
# Inferred from Sequence or structural Similarity (ISS)
# Inferred from Sequence Orthology (ISO)
# Inferred from Sequence Alignment (ISA)
# Inferred from Sequence Model (ISM)
# Inferred from Genomic Context (IGC)
# Inferred from Reviewed Computational Analysis (RCA)
# Author statement evidence codes
# Traceable Author Statement (TAS)
# Non-traceable Author Statement (NAS)
# Curator statement evidence codes
# Inferred by Curator (IC)
# No biological Data available (ND)
# Electronic annotation evidence code
# Inferred from Electronic Annotation (IEA)
termInfo <- 
  function(orgDb, onto,
           evidence=list("Experimental evidence codes"=
                           c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
                             "HTP", "HDA", "HMP", "HGI", "HEP"),
                         "Phylogenetically-inferred annotations"=
                           c("IBA", "IBD", "IKR", "IRD"),
                         "Computational analysis evidence codes"=
                           c("ISS", "ISO", "ISA", "ISM", "IGC", "RCA"),
                         "Author statement evidence codes"=
                           c("TAS", "NAS"),
                         "Curator statement evidence codes"=
                           c("IC", "ND"),
                         "Electronic annotation evidence code"=
                           c("IEA"))
  ){
    stopifnot(is(orgDb, "OrgDb"))
    stopifnot("The 'GO.db' package is required"=
                requireNamespace("GO.db", quietly = TRUE)) 
    onto <- match.arg(onto, choices = c("BP", "CC", "MF"), several.ok = TRUE)
    goAnn <- get(sub(".db", "GO", orgDb$packageName))
    mapped_genes <- mappedkeys(goAnn)
    totalN.genes <- length(unique(mapped_genes))
    message("Using `select()` to retreive the GO terms.")
    all.GO <- select(x=orgDb, keys = mapped_genes,
                     columns = c("GO", "ONTOLOGY"),
                     keytype = "ENTREZID")
    all.GO <- all.GO[all.GO$EVIDENCE %in% unlist(evidence),
                     c("GO", "ONTOLOGY", "ENTREZID"), drop=FALSE]
    all.GO <- unique(all.GO)## in case of duplicated items in database
    ## add Ancestor
    all.GO <- split(all.GO[, c("GO", "ENTREZID")], all.GO$ONTOLOGY)
    onto <- onto[onto %in% names(all.GO)]
    all.GO <- all.GO[onto]
    all.GO <- mapply(FUN=addAncestors,
                     all.GO, names(all.GO),
                     SIMPLIFY = FALSE)
    return(all.GO)
  }
## hypergeometric distribution test
hyperGT <- function(df_all, gene_id, orgDb){
  stopifnot(all(c("GO", "ENTREZID") %in% colnames(df_all)))
  df_sub <- df_all[df_all$ENTREZID %in% gene_id, , drop=FALSE]
  total_gene_count <- length(unique(df_all$ENTREZID))
  this_gene_count <- length(unique(df_sub$ENTREZID))
  each_term_gene_count <- table(df_all$GO)
  each_term_hits <- split(df_sub$ENTREZID, df_sub$GO)
  each_term_hits_count <- lengths(each_term_hits)
  GO <- names(each_term_hits)
  each_term_gene_count <- as.numeric(each_term_gene_count[GO])
  each_term_hits_symbol <- lapply(each_term_hits, FUN=eg2symbol, orgDb=orgDb)
  pvalue <- phyper(q=each_term_hits_count-1,
                   m=each_term_gene_count,
                   n=total_gene_count-each_term_gene_count,
                   k=this_gene_count,
                   lower.tail = FALSE, log.p = FALSE)
  data.frame(GO=GO, 
             pvalue=pvalue,
             fdr=p.adjust(pvalue, method = "BH"),
             countInDataset=each_term_hits_count,
             countInGenome=each_term_gene_count,
             totalGeneInDataset=rep(this_gene_count, length(GO)),
             totalGeneInGenome=rep(total_gene_count, length(GO)),
             geneInDataset=vapply(each_term_hits_symbol,
                                  FUN = paste,
                                  FUN.VALUE = character(1),
                                  collapse=";"))
}
