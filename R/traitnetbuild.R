
#' Build a network matrix based on social attraction values
#'
#' `traitnetbuild()` generates a network (an adjacency matrix) based on social attraction values,
#' such that pairs of individuals with higher social attraction are more likely to get a link, and (for weighted networks) get a larger link weight.
#'
#' NOTE: If the network building is to be used in the generation of trait preference networks, then it is better to use the `traitnet()` or `traitnetsmetrics()` function (which call this function).
#'
#' The `traitnetbuild()` function works in the following way:
#' The number of links is calculated based on the set average degree and number of individuals.
#' Then the dyads to be linked are found by a weighted draw, where dyads with higher social attraction have a higher chance of being picked.
#' The link weights are determined based on the set link type and added to the existing links. If there are not enough dyads with social attraction values above zero to give the input average degree,
#' then random links will be added among the remaining dyads to get the set average degree.
#'
#' About `linktype`:
#' * 'unw' gives unweighted (binary) links. Existing links are indicated with the value 1, non-existing links are indicated with 0.
#' * 'stow' gives weighted links with stochasticity, where the link weights are drawn from distributions with the social attraction values (sociatvals) as means.
#' * 'detw' gives deterministic weighted links, where the link weights equal the social attraction values (sociatvals).
#'
#' @param n An integer indicating the number of individuals (nodes) in the network. Must fit with the length of the sociatvals argument.
#' @param k An integer indicating the average degree of the network, i.e. the average number of links each individual (node) has.
#' @param linktype A string indicating the type of network links ('stow', 'detw', or 'unw'; see Details).
#' @param sociatvals A vector with a social attraction value (a number between 0 and 1) for each pair of individuals (each dyad of nodes)
#'
#' @return An network matrix (adjacency matrix). The matrix is undirected and symmetric, with zeros in the diagonal.
#' @export
#'
#'@references
#' Brask, J. B., Brent, L. J., & De Moor, D. (2026a). okaapi: an R package for generating social networks based on trait preferences. EcoEvoRxiv.
#' \url{https://ecoevorxiv.org/repository/view/9622/}
#'
#' Brask, J. B., Koher, A., Croft, D. P., & Lehmann, S. (2026b). Far-reaching consequences of trait preferences for animal social network structure and function. Behavioral Ecology, 37(1), araf132.
#' \url{https://academic.oup.com/beheco/article/37/1/araf132/8339509}
#'
#' Brask, J. B., Silk, M., & Weiss, M. N. (2025). An introduction to generative network models and their use in animal sociality research. Animal Behaviour, 123364.
#' \url{https://www.sciencedirect.com/science/article/pii/S000334722500291X}
#'
#'
#' @author
#' Josefine Bohr Brask, \email{bohrbrask@@gmail.com}
#'
#' @seealso
#' \code{\link{traitnet()}}
#'
#'
#' @examples
#'
#' # Create input:
#' n <- 100
#' k <- 10
#' linktype <- 'stow'
#'
#' # For the example we use random values as social attraction values:
#' ndyads <- (n*n-n)/2
#' sociatvals <-  runif(ndyads)
#'
#' # Make the network matrix:
#' net <- traitnetbuild(n, k, linktype, sociatvals)



traitnetbuild <- function(n, k, linktype, sociatvals){

  if(length(n)==0 || length(k)==0 || length(linktype)==0 || length(sociatvals)==0) {
    stop("Input without defaults must not be empty")
  }
  if (n <= 0 || is.na(n)) {
    stop('n must be positive.')
  }
  if (k <= 0 || is.na(k)) {
    stop('k must be positive.')
  }
  if (!(linktype %in% c('stow','detw','unw')) ) {
    stop('linktype is incorrect')
  }
  if (length(sociatvals) != ((n*n-n)/2)) {
    stop('Length of sociatvals does not fit with n')
  }

  ndyads <- (n*n-n)/2
  inames <- rep(1:(n-1),(n-1):1)
  jnames <- sequence((n-1):1, from = 2:n)
  L <- (k*n)/2
  unwnet <- matrix(0,n,n)
  linkdyadnums <- wrswoR::sample_int_crank(ndyads, L, sociatvals)
  linkinames <- inames[linkdyadnums]
  linkjnames <- jnames[linkdyadnums]
  linkinamesboth <- c(linkinames, linkjnames)
  linkjnamesboth <- c(linkjnames, linkinames)
  unwnet[cbind(linkinamesboth, linkjnamesboth)] <- 1

  if (linktype =='detw'){
    net <- unwnet
    linkweightsmod <- sociatvals[linkdyadnums]
    linkweightsmod[linkweightsmod==0] <- 0.001
    linkweights <- linkweightsmod
    net[cbind(linkinamesboth, linkjnamesboth)] <- linkweights
  } else if (linktype =='stow'){
    net <- unwnet
    linkwmin <- 0.001
    linkwvar <- 0.05
    rawlinkw <- sociatvals[linkdyadnums]
    rawlinkdraw <- rnorm(length(rawlinkw), mean = rawlinkw, sd = linkwvar)
    linkweightsmod <- rawlinkdraw
    linkweightsmod[linkweightsmod <= 0] <- linkwmin
    linkweightsmod[linkweightsmod > 1] <- 1
    linkweights <- linkweightsmod
    net[cbind(linkinamesboth, linkjnamesboth)] <- linkweights
  }
  else if (linktype =='unw'){
    net <- unwnet
  }
  return(net)
}

