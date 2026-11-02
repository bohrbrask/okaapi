



#' Measure network metrics on a set of networks
#'
#'`traitnetsmetrics()` measures a set of network metrics on an ensemble of networks based on the trait preference model.
#'
#'
#' For details about trait preference networks, see the documentation (help page) for the `traitnet()` function, and the literature in the References section below.
#'
#' The `traittypes`, `allncats`, `preftypes`, and `wvals` arguments must fit together in terms of order (and have the same length),
#' because the function assumes that the entries in the same place in these vectors fit together.
#' For example, the first trait type in the `traittypes` vector will be used with the first preference type in the `preftypes` vector,
#' the first importance value in the `wvals` vector, and the first category number (which may be NA) in the `allncats` vector.
#'
#' About `traittypes`:
#' * 'cate'  gives categorical trait values. The values consist of the integers from 1 to `ncats`, which each indicates a category. Each category is drawn with the same probability.
#' * 'tnorm' gives normally distributed trait values (truncated). The values are drawn from a truncated normal distribution that lie in the interval \[0,1], with mean = 0.5 and standard deviation = 0.25.
#' * 'ranks' gives trait values indicating ranks (e.g. social rank). The values are equally distributed in the interval \[0,1]. This corresponds to ranks indicated by integers from 1 to n, transformed to the interval \[0,1].
#' * 'unif'  gives uniformly distributed trait values. The values are drawn from a uniform distribution with the interval \[0,1].
#' * 'own' must be given for any trait where trait values are given as input by the user. For such traits, the class needs to be specified in the `owntraitclasses` argument.
#'
#' About `preftypes`:
#' * Preference types belong to two general types: 'similarity' and 'popularity' (as described in Brask et al. 2026a and 2026b, see References).
#' * For similarity preferences, the preference depends on the distance (similarity) between trait values.
#' * For popularity preferences, the preference depends on the magnitude of the trait values.
#' * 'simclose' indicates that individuals prefer others that are closer to themselves in terms of their trait value (i.e. they are similar).
#' * 'simfar' indicates that individuals prefer others that are farther from themselves in terms of their trait value (i.e. they are dissimilar).
#' * 'pophigh' indicates that individuals prefer others that have high trait values (i.e. high trait values are popular).
#' * 'poplow' indicates that individuals prefer others that have low trait values (i.e. low trait values are popular).
#'
#' About `linktype`:
#' * 'unw' gives unweighted (binary) links. Existing links are indicated with the value 1, non-existing links are indicated with 0.
#' * 'stow' gives weighted links with stochasticity, where the weights are drawn from distributions with the social attraction values as means.
#' * 'detw' gives deterministic weighted links, where the weight equals the social attraction of the dyad.
#' * For weighted links, the link weights are in the interval \[0,1], and existing links have a minimum weight of 0.001.
#'
#' About `wvals`:
#' * The sum of the `wval` vector, and each value in the vector, must lie in the interval \[0,1], where 0 indicates no importance and 1 indicates full importance.
#'
#' About `owntraitvals` and `owntraitclasses`:
#' * These arguments should only be given if the user wants to provide trait values for one or more of the traits (for example real-world trait values).
#'
#' About `owntraitvals`:
#' * If the user provides trait values for more than one trait, then `owntraitvals` must have one column for each of these traits, and `n` rows.
#' If the user provides trait values for a single trait then `owntraitvals` can be a vector of length `n` or a matrix with `n` rows and 1 column.
#' Categorical trait values should consist of integers that each indicates a category.
#' Continuous trait values should consist of numbers between 0 and 1, and empirical trait values on other scales should be transformed to this scale before being used.
#' The trait values must fit with the trait classes given in the `owntraitclasses` argument.
#'
#' About `owntraitclasses`:
#' * Give 'cate' for categorical traits.
#' * Give 'cont' for continuous traits.
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
#' * 'average degree' gives the average degree.
#' * 'average w. degree' gives the average weighted degree (also called strength).
#'
#' * 'degree assortativity' gives the degree assortativity.
#' * 'w. degree assortativity' gives the weighted degree assortativity.
#'
#' * 'degree variation' gives the degree variation (variance).
#' * 'w. degree variation' gives the weighted degree variation (variance).

#'
#'
#' @param n An integer indicating the number of individuals (network nodes)
#' @param k An integer indicating the average degree (average number of links connected to each node)
#' @param traittypes A vector of strings indicating the trait type for each trait ('cate', 'tnorm', 'ranks', 'unif' or 'own'; see Details)
#' @param allncats A vector of integers indicating the number of categories for each trait (numbers for non-categorical traits are not used and may be set to NA, but must be included in the vector)
#' @param preftypes A vector of strings indicating the preference type to be used with each trait ('simclose', 'simfar', 'pophigh' or 'poplow'; see Details)
#' @param linktype A string indicating the type of link ('stow', 'detw', or 'unw'; see Details)
#' @param wvals  A vector with the importance (weight) of each trait-preference combination
#' @param onecomp Logical. If TRUE then only networks where all the nodes are in a single component (i.e. at least indirectly connected) are used. If FALSE then the networks may consist of more than one component
#' @param nrepls An integer indicating the number of replicates (number of networks for which metrics will be measured)
#' @param metrics A vector of strings indicating the metrics that should be measured (see Details)
#' @param showrepl Logical. If TRUE (default) then the current replicate (network number) is printed to the console
#' @param owntraitclasses A vector of strings indicating the trait class for each trait for which trait values are provided by the user ('cate' or 'cont'; see Details)
#' @param owntraitvals A matrix with `n` rows, and a column for each trait for which trait values are provided by the user (if any). Each column must contain a trait value for each each individual. Can be a vector if there is only a single trait
#'
#' @return a matrix with one row for each network and one column for each metric. The columns are named with the metric names given in the `metrics` argument.
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
#'
#' @seealso \code{\link{netmetrics()}}
#' @seealso \code{\link{traitnet()}}
#'
#' @examples
#'
#' # An example where all arguments are defined directly in the call:
#'
#' traitnetsmetrics(n=100, k=10, traittypes = c('cate', 'tnorm'),
#' allncats = c(5, NA), preftypes = c('simclose', 'pophigh'), linktype = 'stow',
#' wvals = c(0.9, 0.1), onecomp = TRUE, nrepls = 10,
#' metrics = c('clustering', 'average path length', 'average degree'))
#'
#'
#' # An example where some arguments are defined before the call,
#' # and the metric measurements are stored:
#'
#' metrics <- c('average w. path length', 'w. clustering', 'degree assortativity', 'degree variation')
#' traittypes <- c('tnorm', 'cate', 'ranks')
#' preftypes <- c('simfar', 'simclose', 'pophigh')
#' allncats <- c(5, 2, NA)
#' wvals <- c(0.5, 0.2, 0.2)
#' metricresults <- traitnetsmetrics(n=100, k=10, traittypes = traittypes,
#' allncats = allncats, preftypes = preftypes, linktype = 'stow', wvals = wvals,
#' onecomp = TRUE, nrepls = 100, metrics = metrics)
#'








traitnetsmetrics <-  function(n, k, traittypes, allncats, preftypes, linktype, wvals, onecomp, nrepls, metrics, showrepl = TRUE, owntraitclasses = NULL, owntraitvals = NULL){

  if(length(n)==0 || length(k)==0 || length(traittypes)==0 || length(allncats)==0 || length(preftypes)==0 || length(linktype)==0 || length(wvals)==0 || length(onecomp)==0 || length(nrepls)==0 || length(metrics)==0) {
    stop("Input without defaults must not be empty")
  }
  if (n <= 0 || is.na(n)) {
    stop('n must be positive.')
  }
  if (k <= 0 || is.na(k)) {
    stop('k must be positive.')
  }
  if(sum(wvals) > 1) {
    stop("The sum of wvals must not exceed 1.")
  }
  if(any(wvals <0)||any(wvals >1)) {
    stop("All entries in wvals must be between 0 and 1.")
  }
  inputlengths <- c(length(traittypes), length(allncats), length(preftypes), length(wvals))
  if(any(inputlengths != inputlengths[1])) {
    stop("Lengths of input do not fit.")
  }
  if (!(linktype %in% c('stow','detw','unw')) ) {
    stop('linktype is incorrect.')
  }
  if (!all(preftypes %in% c('simclose', 'simfar','pophigh', 'poplow')) ) {
    stop('preftypes are incorrect.')
  }
  if (!all(traittypes %in% c('cate','tnorm','ranks', 'unif', 'own')) ) {
    stop('traittypes is incorrect.')
  }
  if(any(is.na(owntraitvals))) {
    stop("owntraitvals must not contain NA's.")
  }
  if (!is.null(owntraitclasses)) {
    if (!all(owntraitclasses %in% c('cate','cont')) ) {
      stop('owntraitclasses is incorrect.')
    }
    inputlengthsown <- c(length(owntraitclasses), ncol(as.matrix(owntraitvals)), length(which(traittypes == 'own')) )
    if(any(inputlengthsown != inputlengthsown[1])) {
      stop("Input lengths do not fit.")
    }
  }
  if (!all(metrics %in% c('clustering', 'w. clustering', 'average path length', 'average w. path length', 'average degree', 'average w. degree', 'degree assortativity', 'w. degree assortativity', 'degree variation', 'w. degree variation')) ) {
    stop('metrics is incorrect.')
  }
  if (!is.null(owntraitvals)) {
    owntmat <- as.matrix(owntraitvals)
    owntcate <- owntmat[,owntraitclasses == 'cate']
    owntcont <- owntmat[,owntraitclasses == 'cont']
    owntcatediv <- owntcate/round(owntcate)
    if (any(owntcatediv != 1)) {
      stop('Categorical trait values must be integers.')
    }
    if (any(owntcont<0)||any(owntcont>1)) {
      stop('Continuous trait values must be in the range from 0 to 1.')
    }
  }

  metricvals <- matrix(NA, nrow = nrepls, ncol = length(metrics))
  colnames(metricvals) <- metrics

  for (currepl in 1:nrepls){
    traitnetoutput <- traitnet(n = n, k = k, traittypes = traittypes, allncats = allncats, preftypes = preftypes, linktype = linktype, wvals = wvals, onecomp = onecomp, visnet = FALSE, owntraitclasses = owntraitclasses, owntraitvals = owntraitvals)
    curnet <- traitnetoutput$net
    curmetricvals <- netmetrics(curnet, metrics)
    metricvals[currepl, ] <- curmetricvals
    if (showrepl==TRUE){
      print(currepl)
    }
  }

  return(metricvals)

}

