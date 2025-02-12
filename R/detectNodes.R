# help function getAnchor
getAnchor <- function(gr, region, ...){
  ol <- findOverlaps(gr, region, ...)
  names(region)[subjectHits(ol[match(seq_along(gr), queryHits(ol))])]
}

#' Detect the interaction node
#' @description Define the interaction node from input Pairs.
#' @param interaction An object of \link[S4Vectors:Pairs-class]{Pairs} to
#' represent interactions.
#' @param pval_cutoff Cutoff P value for interaction node by Poisson distribution
#' @param ... Not used.
#' @return A list of interaction nodes with elements:
#'  node_connection, Pairs object represent interactions interacted with nodes;
#'  nodes, GRanges object represent regions with maximal interactions involved in nodes;
#'  node_regions, GRanges object represent regions interacted with nodes.
#' @export
#' @importMethodsFrom S4Vectors first second mcols mcols<-
#' @importClassesFrom S4Vectors Pairs Hits
#' @importClassesFrom graph graphNEL
#' @importClassesFrom GenomicRanges GRangesList
#' @importMethodsFrom IRanges reduce findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom GenomicRanges GRangesList
#' @importFrom graph ugraph
#' @importFrom RBGL connectedComp
#' @importFrom stats ppois
#' @importFrom methods is new
#' @examples
#' library(rtracklayer)
#' p <- system.file("extdata", "WT.2.bedpe",
#'                  package = "GenomicInteractionNodes")
#' interactions <- import(con=p, format="bedpe")
#' nodes <- detectNodes(interactions)
detectNodes <- function(interaction, pval_cutoff=0.05, ...){
  stopifnot(is(interaction, "Pairs"))
  regions <- reduce(c(first(interaction), second(interaction)))
  stopifnot("Maximal regions of interaction should be less than 10G"=
              length(regions)<=1e10)
  names(regions) <- paste0("p", seq_along(regions))
  option_scipen <- options(scipen=10) ## maximal interaction number is 1e10
  if(is.list(option_scipen)){
    option_scipen <- option_scipen$scipen
  }
  if(!is.numeric(option_scipen)){
    option_scipen <- 0
  }
  on.exit(options(scipen=option_scipen))
  olm <- data.frame(l=getAnchor(first(interaction), regions),
                    r=getAnchor(second(interaction), regions))
  stopifnot("unexpect happend at check point 1: wrong anchor number"=
              length(interaction)==nrow(olm))
  olm_cp <- olm[, c("r", "l")]
  colnames(olm_cp) <- colnames(olm)
  olm <- rbind(olm, olm_cp)
  edgeL <- split(olm[,2], olm[,1])
  edgeL <- lapply(edgeL, unique)
  nodes <- names(regions)
  ## there is an issue, the nodes and edgeL have limits
  gR <- new("graphNEL", nodes=nodes, edgeL=edgeL)
  Merged <- connectedComp(ugraph(gR))
  len <- lengths(Merged)
  ## cut by p value, giant component from Erdos-Renyi random graph
  ## the G(n,p) model, due to Erdos and Renyi, has two parameters, n and p.
  ## Here n is the number of vertices of the graph and p is the edge probability.
  ## For each pair of distinct vertices, v and w, p is the probability that
  ## the edge (v, w) is present.
  ## p=d/n, d=d(n-1)/n, if d>1, there is giant component randomly
  n <- length(nodes)
  # length of the edges
  l <- sum(lengths(edgeL))
  # p = l/n
  p <- l/n
  p <- ppois(len, lambda=p, lower.tail=FALSE)
  components <- Merged[p<pval_cutoff]
  components <- lapply(components, function(.ele) regions[.ele])
  ## get the interactions in the components
  components_ids <- rep(names(components), lengths(components))
  node_connection <- unlist(GRangesList(components))
  node_connection$comp_id <- components_ids
  ol <- findOverlaps(node_connection, interaction)
  node_connection <- data.frame(connection_id=subjectHits(ol),
                                comp_id=node_connection$comp_id[queryHits(ol)])
  node_connection <- unique(node_connection)
  node_conn <- interaction[node_connection$connection_id]
  mcols(node_conn)$comp_id <- node_connection$comp_id
  ## find the node of the interactions
  sub_group <- c(first(node_conn), second(node_conn))
  sub_group$comp_id <- rep(mcols(node_conn)$comp_id, 2)
  node_regions <- unique(sub_group)
  ol <- findOverlaps(sub_group, node_regions, type = "equal")
  stopifnot("unexpect happend at check point 2: wrong node number"=
              length(ol)==length(sub_group))
  node_name <- data.frame(node_region_id=subjectHits(ol),
                          node_comp_id=sub_group$comp_id[queryHits(ol)])
  node_name_table <- table(node_name)
  node_name <- apply(node_name_table, 2, which.max, simplify = FALSE)
  node_name <-
    cbind(comp_id=rep(names(node_name), lengths(node_name)),
          node_region_id=rownames(node_name_table)[unlist(node_name)])
  mode(node_name) <- "integer"
  nodes <- node_regions[node_name[, "node_region_id"]]
  nodes$comp_id <- node_name[, "comp_id"]
  list(node_connection=node_conn, nodes=nodes, node_regions=node_regions)
}