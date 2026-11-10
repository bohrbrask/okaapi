

#' Measure network metrics on a single network
#'
#'`netmetrics()` measures a set of network metrics on an input network.
#'
#'
#'#' NOTE: To measure network metrics on trait preference networks, it is better to use the `traitnetsmetrics()` function (which calls this function).
#'
#'The `netmetrics()` function uses the `igraph` package for some of the metric calculations.
#'

#' About `metrics`:
#'
#' * 'clustering' gives the global clustering coefficient, also called transitivity (measured as the ratio of triangles to connected triples).
#' * 'w. clustering' gives the weighted clustering coefficient (as defined in Barrat et al. 2004, see References). Nodes with one link get a weighted clustering of zero.
#'
#' * 'average path length' gives the average path length, i.e. the average unweighted distance between nodes.
#' * 'average w. path length' gives the average weighted path length, i.e. the average weighted distance between nodes (with smaller distances for larger weights).
#'    - For networks with more than one component: nodes in different network components are considered to have no path between them and these dyads therefore do not contribute to path length calculations.
#'
#' * 'average degree' gives the average degree, i.e. the average number of links connected to each node.
#' * 'average w. degree' gives the average weighted degree (also called strength), i.e. the sum of the weights of the links connected to each node.
#'
#' * 'degree assortativity' gives the degree assortativity.
#' * 'w. degree assortativity' gives the weighted degree assortativity.
#'
#' * 'degree variation' gives the degree variation (variance).
#' * 'w. degree variation' gives the weighted degree variation (variance).
#'
#'
#' @param net The network (adjacency matrix) to be plotted.The matrix must be undirected and symmetric, with zeros in the diagonal.
#' @param metrics A vector of strings indicating the metrics that should be measured (see Details).
#'
#' @return A vector with a value for each of the metrics specified in the `metrics` argument. The cells are named with the metric names given in the `metrics` argument.

#' @export

#' @references
#' Alain Barrat, Marc Barthelemy, Romualdo Pastor-Satorras, Alessandro Vespignani: The architecture of complex weighted networks, Proc. Natl. Acad. Sci. USA 101, 3747 (2004)
#'
#' @author
#' Josefine Bohr Brask, \email{bohrbrask@@gmail.com}
#'
#' @seealso  \code{\link{traitnet()}}
#' @seealso  \code{\link{traitnetsmetrics()}}
#' @seealso  \code{\link{igraph}}
#'
#' @examples
#'
#' # An example that uses a network created with the traitnet function (this is not recommended, it is better to use the traitnetsmetrics function):
#'
#' traitnetoutput <- traitnet(n=100, k=10, traittypes = c('cate', 'tnorm'),
#' allncats = c(5, NA), preftypes = c('simclose', 'pophigh'), linktype = 'stow',
#' wvals = c(0.8, 0.1), onecomp = TRUE, visnet = TRUE)
#'
#' net <- traitnetoutput$net
#'
#' metrics <- c('clustering', 'w. clustering', 'average path length', 'average w. path length', 'average degree', 'average w. degree',
#' 'degree assortativity', 'w. degree assortativity', 'degree variation', 'w. degree variation')
#'
#' netmetrics(net, metrics)
#'
#'
#' # An example of calculating clustering and path length on
#' two common network types made with igraph:
#'
#'library(igraph)
#'
#' # A random network
#' randnet <-  as_adjacency_matrix(sample_gnp(100, 0.1), type = 'both',
#' names = FALSE, sparse = FALSE)
#'
#' # A scale-free network
#' scalenet <-  as_adjacency_matrix(sample_pa(100, directed = FALSE),
#' type = 'both', names = FALSE, sparse = FALSE)
#'
#' # Calculate clustering and path length for the two networks
#' metrics <- c('clustering', 'average path length')
#' randnetmetrics <- netmetrics(randnet, metrics)
#' scalenetmetrics <- netmetrics(scalenet, metrics)









netmetrics <- function(net, metrics){

  if(length(net)==0 || length(netmetrics)==0 ) {
    stop("Input must not be empty")
  }
  if(any(diag(net)!=0)) {
    stop("The input network must have only zeros in the diagonal.")
  }
  if(!isSymmetric(net)) {
    stop("The input network must be symmetric.")
  }
  if(any(is.na(net))) {
    stop("The input network must not contain NA's.")
  }
  if (!all(metrics %in% c('clustering', 'w. clustering', 'average path length', 'average w. path length', 'average degree', 'average w. degree', 'degree assortativity', 'w. degree assortativity', 'degree variation', 'w. degree variation')) ) {
    stop('metrics is incorrect.')
  }

  metricvals <- rep(NA, length(metrics))
  names(metricvals) <- metrics
  netigraph <- igraph::graph_from_adjacency_matrix(net, mode = "undirected", diag = FALSE, weighted = TRUE) # get igraph version of network
  degs <- colSums(net>0)
  weidegs <- colSums(net)

  if ('clustering' %in% metrics){
    metricvals['clustering'] <- igraph::transitivity(netigraph, type= "global") # calculate and store unweighted transitivity (global clustering = ratio of triangles and connected triples)
  }
  if ('w. clustering' %in% metrics){
    transes <- igraph::transitivity(netigraph, type= "weighted")
    transesnonan <- transes
    transesnonan[transesnonan =='NaN'] <- 0 #we set NaN values to zero (individuals with only one link get NaN but we give them zero instead)
    metricvals['w. clustering'] <- mean(transesnonan)     # weighted clustering
  }
  if ('average path length' %in% metrics){
    metricvals['average path length'] <- igraph::mean_distance(netigraph, directed = FALSE, weights = NA) # calculate unweighted mean distance (average path length)
  }
  if ('average w. path length' %in% metrics){
    nonzeroweights <- net[net>0] # get nonzero links
    invnonzeroweights <- 1/nonzeroweights # get the inverse of the nonzero edge weights
    nonzeroplaces <- which(net>0) # find places of nonzero links
    invnet <- net
    invnet[nonzeroplaces] <- invnonzeroweights #insert the inverted link strengths in place of the originals
    invnetigraph <- igraph::graph_from_adjacency_matrix(invnet, mode = "undirected", diag = FALSE, weighted = TRUE)
    metricvals['average w. path length'] <- igraph::mean_distance(invnetigraph, directed = F)
  }
  if ('average degree' %in% metrics){
    metricvals['average degree'] <- mean(degs) # calculate mean degree (unweighted)
  }
  if ('average w. degree' %in% metrics){
    metricvals['average w. degree'] <- mean(weidegs) # calculate mean weighted degree
  }
  if ('degree assortativity' %in% metrics){
    metricvals['degree assortativity'] <- igraph::assortativity_degree(netigraph, directed = F)# calculate degree assortativity
  }
  if ('w. degree assortativity' %in% metrics){
    metricvals['w. degree assortativity'] <- igraph::assortativity(netigraph, weidegs, directed = F)# calculate weighted degree assortativity
  }
  if ('degree variation' %in% metrics){
    metricvals['degree variation'] <- var(degs)   # calculate degree variation
  }
  if ('w. degree variation' %in% metrics){
    metricvals['w. degree variation'] <- var(weidegs)   # calculate weighted degree variation
  }

  return(metricvals)

}


