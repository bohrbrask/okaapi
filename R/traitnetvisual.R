



#' Plot a trait preference network
#'
#' `traitnetvisual()` makes a plot of a trait preference network, where the sizes and colours of the nodes correspond to trait values.
#'
#' NOTE: to plot trait preference networks, it is better to use the `traitnet()` function (which calls this function).
#'
#' The `traitnetvisual()` function uses two methods for determining the sizes and colours of nodes:
#' * 'Default method (when the `vistraitnum` argument is not set): The node size is based on the popularity trait (if any) with highest importance (wval). The node colour is based on the similarity trait (if any) with highest importance (wval).
#' * 'User-specified method (when the `vistraitnum` argument is set): The user sets the traits to be used for node sizes and colouring. The user can also set one of them, in which case the other (size or colour) will be uniform rather than correspond to a trait.
#'
#' Furthermore, to make a nice visualization, the average node size is adjusted based on the number of nodes, and the link width is adjusted based on whether the network is weighted or not. The layout itself is done with `igraph` plotting.
#'
#'
#' Generally about function arguments:
#'
#' All the arguments other than `net` must be the ones that were used for the creation of the input network (`net`), otherwise the visualization will be misleading (a safer option is to use the `traitnet()` function, which calls this function).
#'
#'
#' About `traitclasses`:
#' * 'cate' indicates a categorical trait.
#' * 'cont' indicates a continuous trait.
#'
#' About `preftypes`:
#'
#' Preference types belong to two general types: 'similarity' and 'popularity' (as described in Brask et al. 2026a and 2026b, see References).
#' For similarity preferences, the preference depends on the distance (similarity) between trait values.
#' For popularity preferences, the preference depends on the magnitude of the trait values.
#' * 'simclose' indicates that individuals prefer others that are closer to themselves in terms of their trait value (i.e more similar; homophily).
#' * 'simfar' indicates that individuals prefer others that are farther from themselves in terms of their trait value (i.e. less similar<, heterophily).
#' * 'pophigh' indicates that individuals prefer others that have high trait values (i.e. high trait values are popular).
#' * 'poplow' indicates that individuals prefer others that have low trait values (i.e. low trait values are popular).
#'
#' About `linktype`:
#' * 'unw' gives unweighted (binary) links. Existing links are indicated with the value 1, non-existing links are indicated with 0.
#' * 'stow' gives weighted links with stochasticity, where the link weights are drawn from distributions with the social attraction values (sociatvals) as means.
#' * 'detw' gives deterministic weighted links, where the link weights equal the social attraction values (sociatvals).
#'
#' About `traitvals`:
#'
#' The trait values must fit with the trait classes given in the `traitclasses` argument.
#' Categorical trait values should consist of integers that each indicates a category.
#' Continuous trait values should consist of numbers between 0 and 1. Empirical trait values on other scales should be transformed to this scale before being used.
#'
#'
#' @param net The network (adjacency matrix) to be plotted.The matrix must be undirected and symmetric, with zeros in the diagonal.
#' @param traitclasses A vector of strings indicating the trait class for each trait ('cate' or 'cont'; see Details).
#' @param preftypes A vector of strings indicating the preference type to be used with each trait ('simclose', 'simfar', 'pophigh' or 'poplow'; see Details).
#' @param linktype A string indicating the type of network links ('stow', 'detw', or 'unw'; see Details).
#' @param wvals A vector with numbers between 0 and 1 (must sum to 1 or less), indicating the importance (weight) of each trait preference.
#' @param traitvals A matrix with a column of trait values for each trait, and a row for each individual. Can be a vector if there is only a single trait (see Details).
#' @param vistraitnum If not set then the default plotting method is used (see Details). A vector of two integers indicating the trait to be used for the colouring (first integer) and sizes (second integer) of nodes in the network plot, with the integers refering to the places of the traits in the traittypes vector. Either can also be NA.
#'
#' @return A network plot.
#' @export

#' @references
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
#' # Here we first make a trait preference network by using `traitvalues()`,
#' # `traitnetsociat()` and `traitnetbuild()`, and then we visualize it.
#' # This is not a recommended procedure, as you can make and visualize a trait
#' # preference network with the `traitnet()` function alone, which is safer.
#'
#' # set arguments
#' n <- 100
#' k <- 10
#' traittype <- 'cate'
#' ncats <- 5
#' traitclasses <- 'cate'
#' preftypes <- 'simclose'
#' wvals <- 0.95
#' linktype <- 'stow'
#'
#' # create the network
#' traitvals <- traitvalues(n = 100, traittype = 'cate', ncats = 5)
#' sociatvals <- traitnetsociat(traitclasses, preftypes, wvals, traitvals)
#' net <- traitnetbuild(n, k, linktype, sociatvals)
#'
#' #visualize the network
#' traitnetvisual(net, traitclasses, preftypes, linktype, wvals, traitvals)







traitnetvisual <- function(net, traitclasses, preftypes, linktype, wvals, traitvals, vistraitnum = NULL) {

  if(any(is.na(net))) {
    stop("The input network must not contain NA's.")
  }
  if(any(diag(net)!=0)) {
    stop("The input network must have only zeros in the diagonal.")
  }
  if(!isSymmetric(net)) {
    stop("The input network must be symmetric.")
  }
  if(sum(wvals) > 1) {
    stop("The sum of wvals must not exceed 1.")
  }
  if(any(wvals <0)||any(wvals >1)) {
    stop("All entries in wvals must be between 0 and 1.")
  }
  inputlengths <- c(length(traitclasses), length(preftypes), length(wvals), ncol(traitvals))
  if(any(inputlengths ==0)) {
    stop("Input without defaults must not be empty.")
  }
  if(any(inputlengths != length(traitclasses))) {
    stop("Input lengths do not fit.")
  }
  if(nrow(net) != nrow(as.matrix(traitvals))) {
    stop("net and traitvals must have the same number of rows.")
  }
  if (!(linktype %in% c('stow','detw','unw')) ) {
    stop('Linktype is incorrect.')
  }
  if (!all(traitclasses %in% c('cate','cont')) ) {
    stop('traitclasses are incorrect.')
  }
  if(any(is.na(traitvals))) {
    stop("traitvals must not contain NA's.")
  }
  if (!all(preftypes %in% c('simclose', 'simfar','pophigh', 'poplow')) ) {
    stop('preftypes are incorrect.')
  }
  tmat <- as.matrix(traitvals)
  tcate <- tmat[,traitclasses == 'cate']
  tcont <- tmat[,traitclasses == 'cont']
  tcatediv <- tcate/round(tcate)
  if (any(tcatediv != 1)) {
    stop('Categorical trait values must be integers.')
  }
  if (any(tcont<0)||any(tcont>1)) {
    stop('Continuous trait values must be in the range from 0 to 1.')
  }


  # if vistraitnum not lenght 2 stop
  # if vistraitnum not integers stop
  # if vistraitnum higher than number of traits stop



  n <- nrow(net)
  traitvals <- as.matrix(traitvals)
  slope <- -0.01
  intercept <- 15
  nodesizescale <- slope*n+intercept
  if (linktype == 'unw') {
    linksizescale <- 1
  } else if (linktype == 'detw' || linktype == 'stow' ) {
    linksizescale <- 3
  }
  igraphnet = igraph::graph_from_adjacency_matrix(net, mode = "undirected", diag = FALSE, weighted = T)
  igraph::V(igraphnet)$frame.color <- 'white'
  igraph::V(igraphnet)$label <- NA




  if (!is.null(vistraitnum)) {

    if (!is.na(vistraitnum[1])) {

    sizeplace <- vistraitnum[1]
    sizetraitclass <- traitclasses[sizeplace]
    sizetraitvals <- traitvals[,sizeplace]
    if (sizetraitclass == 'cate') {
      igraph::V(igraphnet)$size <- nodesizescale*sizetraitvals/max(sizetraitvals)
    } else if (sizetraitclass =='cont') {
      igraph::V(igraphnet)$size <- nodesizescale*sizetraitvals+3
    }
    }

    if (!is.na(vistraitnum[2])) {
    coloplace <- vistraitnum[2]
    coloclass <- traitclasses[coloplace]
    tcolo <- traitvals[,coloplace]
    if (coloclass == 'cate') {
      igraph::V(igraphnet)$color <- tcolo
      ncolocats <- length(unique(tcolo))
      igraphnet$palette <- igraph::sequential_pal(ncolocats)
    } else if (coloclass == 'cont') {
      reddishcols = c('#fff7ec','#fee8c8','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#b30000','#7f0000')
      igraph::V(igraphnet)$color <- contincols(cols = reddishcols, colvals = tcolo)
    }
  } else {
    igraph::V(igraphnet)$color <- '#b30000'
  }





  } else if (is.null(vistraitnum)) {

  if(any(preftypes %in% c('pophigh', 'poplow'))){
    popplaces <- which(preftypes %in% c('pophigh', 'poplow'))
    popwvals <- wvals[popplaces]
    maxpopplacecur <- which.max(popwvals)
    maxpopplace <- popplaces[maxpopplacecur]
  popclass <- traitclasses[maxpopplace]
   tpop <- traitvals[,maxpopplace]
   if (popclass == 'cate') {
     igraph::V(igraphnet)$size <- nodesizescale*tpop/max(tpop)
   } else if (popclass =='cont') {
     igraph::V(igraphnet)$size <- nodesizescale*tpop +3
 }
 }


  if(any(preftypes %in% c('simclose', 'simfar'))){
    simplaces <- which(preftypes %in% c('simclose', 'simfar'))
    simwvals <- wvals[simplaces]
    maxsimplacecur <- which.max(simwvals)
    maxsimplace <- simplaces[maxsimplacecur]
    simclass <- traitclasses[maxsimplace]
    tsim <- traitvals[,maxsimplace]

    if (simclass == 'cate') {
      igraph::V(igraphnet)$color <- tsim
      nsimcats <- length(unique(tsim))
      igraphnet$palette <- igraph::sequential_pal(nsimcats)
    } else if (simclass == 'cont') {
      reddishcols = c('#fff7ec','#fee8c8','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#b30000','#7f0000')
      igraph::V(igraphnet)$color <- contincols(cols = reddishcols, colvals = tsim)
    }
  } else {
    igraph::V(igraphnet)$color <- '#b30000'
  }


  }


 traitnetplot <- igraph::plot.igraph(igraphnet, edge.width=igraph::E(igraphnet)$weight*linksizescale)



  return(traitnetplot)

}
