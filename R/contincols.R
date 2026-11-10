
#' Create colours based on values from a continuous scale
#'
#' `contincols()` creates a set of colours based on a set of values and a colour scale, such that the distances between the colours correspond to the distances between the values.
#' The colours can be used to colour the nodes of an `igraph` network plot, based on a node attribute.
#'
#' NOTE: If the colours are to be used in the plotting of trait preference networks, then it is better to use the `traitnet()` function (which calls this function).
#'
#'
#' @param cols A vector with two or more colours from which the colour scale is created by extrapolation between the colours. For valid formats see the col argument of the `col2rgb()` function.
#' @param colvals A vector with values between 0 and 1. If the colours are used to colour the nodes of a network, then the length of the vector should equal the number of network nodes.
#'
#' @return A vector with colours that can be used directly to set the node colours of an `igraph` network plot.
#' @export
#'
#' @author
#' Josefine Bohr Brask, \email{bohrbrask@@gmail.com}
#'
#' @seealso
#' \code{\link{igraph}}
#'
#'
#' @examples
#'
#' # An example of colouring an igraph network based on node attribute values:
#'
#' # set colours for colour scale:
#' cols <- c('#CCFFFF','#003366')
#'
#' # Make simulated node attribute values:
#' colvals <- c(rep(0.01,30), rep(0.3,30), rep(0.99,40))
#'
#' # create colours for the network nodes, based on their attribute values:
#' nodecols <- contincols(cols, colvals)
#'
#'# make a network and colour it with the created colours:
#' library(igraph)
#' net <- sample_gnp(100,0.1)
#' V(net)$color <- nodecols
#' plot.igraph(net)
#'
#'
#' # The same with the middle value changed (note the corresponding colour change):
#'
#' # set colours for colour scale:
#' cols <- c('#CCFFFF','#003366')
#'
#' # Make simulated node attribute values:
#' colvals <- c(rep(0.01,30), rep(0.7,30), rep(0.99,40))
#'
#' # create colours for the network nodes, based on their attribute values:
#' nodecols <- contincols(cols, colvals)
#'
#'# make a network and colour it with the created colours:
#' library(igraph)
#' net <- sample_gnp(100,0.1)
#' V(net)$color <- nodecols
#' plot.igraph(net)






contincols <- function(cols, colvals){
  if(any(is.na(cols))) {
    stop("cols must not contain NA's")
  }
  if(any(colvals>1) || any(colvals<0)) {
    stop("Values in colvals must be between 0 and 1")
  }
  if(length(cols)<2) {
    stop("cols must contain at least two values")
  }
  if(length(colvals)==0) {
    stop("colvals must not be empty")
  }
  colfunction <- colorRamp(cols)
  nodecolsraw <- colfunction(colvals)
  nodecols <- rgb(nodecolsraw[,1], nodecolsraw[,2], nodecolsraw[,3], maxColorValue = 255)
  return(nodecols)
}
