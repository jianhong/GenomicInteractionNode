# help function getAnchor
getAnchor <- function(gr, region, ...){
  ol <- findOverlaps(gr, region, ...)
  names(region)[subjectHits(ol[match(seq_along(gr), queryHits(ol))])]
}

#' Detect the interaction hub
#' @description Define the interaction hub from input Pairs.
#' @param interaction An object of \link[S4Vectors:Pairs-class]{Pairs} to
#' represent interactions.
#' @param pval_cutoff Cutoff P value for interaction hub by poisson distribution
#' @param ... Not used.
#' @return A list of interaction hubs with elements:
#'  hub_connection, Pairs object represent interactions interacted with hubs;
#'  hubs, GRanges object represent regions involved in hubs;
#'  hub_regions, GRanges object represent regions interacted with hubs.
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
#'                  package = "interactionHub")
#' interactions <- import(con=p, format="bedpe")
#' hubs <- detectHubs(interactions)
detectHubs <- function(interaction, pval_cutoff=0.05, ...){
  stopifnot(is(interaction, "Pairs"))
  regions <- reduce(c(first(interaction), second(interaction)))
  stopifnot("Maximal regions of interaction should be less than 10G"=
              length(regions)<=1e10)
  names(regions) <- paste0("p", seq_along(regions))
  option_scipen=options(scipen=10) ## maximal interaction number is 1e10
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
  gR <- new("graphNEL", nodes=nodes, edgeL=edgeL)## there is an issue, the nodes and edgeL have limits
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
  hub_connection <- unlist(GRangesList(components))
  hub_connection$comp_id <- components_ids
  ol <- findOverlaps(hub_connection, interaction)
  hub_connection <- data.frame(connection_id=subjectHits(ol),
                               comp_id=hub_connection$comp_id[queryHits(ol)])
  hub_connection <- unique(hub_connection)
  hub_conn <- interaction[hub_connection$connection_id]
  mcols(hub_conn)$comp_id <- hub_connection$comp_id
  ## find the hub of the interactions
  sub_group <- c(first(hub_conn), second(hub_conn))
  sub_group$comp_id <- rep(mcols(hub_conn)$comp_id, 2)
  hub_regions <- unique(sub_group)
  ol <- findOverlaps(sub_group, hub_regions, type = "equal")
  stopifnot("unexpect happend at check point 2: wrong hub number"=
              length(ol)==length(sub_group))
  hub_name <- data.frame(hub_region_id=subjectHits(ol),
                         hub_comp_id=sub_group$comp_id[queryHits(ol)])
  hub_name_table <- table(hub_name)
  hub_name <- apply(hub_name_table, 2, which.max, simplify = FALSE)
  hub_name <- cbind(comp_id=rep(names(hub_name), lengths(hub_name)),
                    hub_region_id=rownames(hub_name_table)[unlist(hub_name)])
  mode(hub_name) <- "integer"
  hubs <- hub_regions[hub_name[, "hub_region_id"]]
  hubs$comp_id <- hub_name[, "comp_id"]
  list(hub_connection=hub_conn, hubs=hubs, hub_regions=hub_regions)
}