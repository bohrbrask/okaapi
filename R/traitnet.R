
#' Generate and visualize a network based on trait preferences
#'
#'`traitnet()` makes a network based on trait preferences, and if desired a plot of the network.
#'
#' This function can be used both as a tool to generate networks based on trait preferences (or more broadly node attribute effects), and as a tool to generate networks with specific structural properties (potentially disregarding the trait preference interpretation).
#'
#' Details about trait preferences and generation of trait preference networks (including case studies) are given in the papers listed in References.
#'
#'
#' Key function features:
#'
#' The function can create networks going from randomly structured to strongly influenced by trait preferences. Different trait types and preference types can be used. multiple trait preferences can influence the network structure simultaneously, and can have different importance. The trait values of the individuals can be simulated by the function, or given as input by the user.
#'
#' Network structures:
#'
#' The function can generate networks with a range of global structural patterns, including random, modular, assorted, disassorted, centralized, and mixes of these in any ratio.
#'
#' Network generation:
#'
#' The networks are generated using the trait preference model (described in References). The model calculates a social attraction value for each pair of individuals, based on their trait values and the preference type and importance set for each trait. From this, a network is constructed where pairs with higher social attraction have a higher chance of getting a link (and a higher link weight, in case of weighted networks).
#'
#' Network plotting:
#'
#' The colour and size of nodes in the network are based on the trait values of the individuals. The function uses two methods for determining the sizes and colours of nodes:
#' * 'Default method (when the `vistraitnum` argument is not set): The node size is based on the popularity trait (if any) with highest importance (wval). The node colour is based on the similarity trait (if any) with highest importance (wval).
#' * 'User-specified method (when the `vistraitnum` argument is set): The user sets the traits to be used for node sizes and colouring. The user can also set one of them, in which case the other (size or colour) will be uniform rather than correspond to a trait.
#'
#' Furthermore, to make a nice visualization, the average node size is adjusted based on the number of nodes, and the link width is adjusted based on whether the network is weighted or not. The layout itself is done with `igraph` plotting.
#'
#'
#' Generally about function arguments:
#'
#' The `traittypes`, `allncats`, `preftypes`, and `wvals` arguments must fit together in terms of the order of their contents (and have the same length),
#' because the function assumes that the entries in the same place in these vectors fit together.
#' For example, the first trait type in the `traittypes` vector will be used with the first preference type in the `preftypes` vector,
#' the first importance value in the `wvals` vector, and the first category number (which may be NA) in the `allncats` vector.
#'
#' About `traittypes`:
#'
#' Different trait types correspond to different distributions. A given trait type may be used to model any specific trait that fits with the distribution.
#' * 'cate'  gives categorical trait values. The values consist of the integers from 1 to the number of categories, with each integer ndicating a category. Each category is drawn with the same probability.
#' * 'tnorm' gives normally distributed trait values (truncated). The values are drawn from a truncated normal distribution that lie in the interval \[0,1], with mean = 0.5 and standard deviation = 0.25.
#' * 'ranks' gives trait values indicating ranks (e.g. social rank). The values are equally distributed in the interval \[0,1]. This corresponds to ranks indicated by integers from 1 to n, transformed to the interval \[0,1].
#' * 'unif'  gives uniformly distributed trait values. The values are drawn from a uniform distribution with the interval \[0,1].
#' * 'own' must be given for any trait where trait values are provided by the user.
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
#' * 'stow' gives weighted links with stochasticity, where the link weights are drawn from distributions with the social attraction values (sociatvals, see Value) as means.
#' * 'detw' gives deterministic weighted links, where the link weights equal the social attraction values (sociatvals, see Value).
#'
#' For weighted links, the link weights are in the interval \[0,1], and existing links have a minimum weight of 0.001.
#'
#' About `owntraitvals` and `owntraitclasses`:
#'
#' These arguments should only be given if the user wants to provide trait values for one or more of the traits (for example real-world trait values).
#'
#' About `owntraitvals`:
#'
#' If the user provides trait values for more than one trait, then `owntraitvals` must have one column for each of these traits, and `n` rows.
#' If the user provides trait values for a single trait, then `owntraitvals` can be a vector of length `n` or a matrix with one column and `n` rows.
#' Categorical trait values must consist of integers that each indicates a category.
#' Continuous trait values must consist of numbers between 0 and 1 (the user must transform trait values on other scales to this scale before giving them as input).
#'
#' About `owntraitclasses`:
#' * Give 'cate' for categorical traits.
#' * Give 'cont' for continuous traits.
#'
#' About `sociatvals`:
#'
#' For each pair of individuals (nodes), a social attraction value is calculated based on their trait values and the set trait preferences.
#' The social attraction values are used in the network construction, with pairs with a larger social attraction value having a higher chance of getting a link (and a larger link weight in case of weighted networks).


#' @param traittypes A vector of strings indicating the trait type for each trait ('cate', 'tnorm', 'ranks', 'unif' or 'own'; see Details).
#' @param allncats A vector of integers indicating the number of categories for each trait (must be set to NA for non-categorical traits).
#' @param preftypes A vector of strings indicating the preference type to be used with each trait ('simclose', 'simfar', 'pophigh' or 'poplow'; see Details).
#' @param wvals A vector with numbers between 0 and 1 (must sum to 1 or less), indicating the importance (weight) of each trait preference.
#' @param n An integer indicating the number of individuals (nodes) in the network.
#' @param k An integer indicating the average degree of the network, i.e. the average number of links each individual (node) has.
#' @param linktype A string indicating the type of network links ('stow', 'detw', or 'unw'; see Details).
#' @param onecomp Logical. If TRUE then the network will consist of a single component, if FALSE then the network may consist of more than one component.
#' @param visnet Logical. If TRUE then a plot of the network is made, if false then no plot is made (see Details and Value).
#' @param vistraitnum If not set then the default plotting method is used (see Details). A vector of two integers indicating the trait to be used for the colouring (first integer) and sizes (second integer) of nodes in the network plot, with the integers refering to the places of the traits in the traittypes vector. Either can also be NA.
#' @param owntraitclasses Only set if any traits are input by the user (see Details). A vector of strings indicating the trait class for each user-provided trait ('cate' or 'cont'; see Details).
#' @param owntraitvals Only set if any traits are input by the user (see Details). A matrix with trait values for user provided traits (one column for each user-provided trait). Can be a vector if there is only one user-provided trait; see Details.
#'
#' @return A plot of the network (only if `visnet` = TRUE) and a list containing the following:
#' * `net`: The network matrix.
#' * `traitvals`: A matrix with all the trait values (one for for each individual (node) for each trait).
#' * `sociatvals`: A vector with all the social attraction values (one for each pair of individuals (each pair of nodes); see Details).


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
#' \code{\link{traitnetsmetrics()}}
#'
#'
#' @examples
#'
#' # More extensive examples can be found in the papers given in References.
#'
#'
#' # A network based on two trait preferences, which have the same importance:
#'
#' traitnet(traittypes = c('cate', 'tnorm'), allncats = c(5, NA),
#' preftypes = c('simclose', 'pophigh'), wvals = c(0.5, 0.5), n=100, k=10, linktype = 'stow',
#' onecomp = TRUE, visnet = TRUE)
#'
#'
#' # The same as the above example, but where one trait preference is of
#' # high importance and the other trait preference is of no importance:
#'
#' traitnet(traittypes = c('cate', 'tnorm'), allncats = c(5, NA),
#' preftypes = c('simclose', 'pophigh'), wvals = c(0.97, 0.00), n=100, k=10, linktype = 'stow',
#' onecomp = TRUE, visnet = TRUE)
#'
#'
#' # A network based on three trait preferences, where we specify arguments before the call,
#' # and we store the network matrix afterwards:
#'
#' traittypes <- c('cate', 'tnorm', 'ranks')
#' allncats <- c(5, NA, NA)
#' preftypes <- c('simclose', 'pophigh', 'simfar')
#' wvals <- c(0.5, 0.1, 0.3)
#' n <- 100
#' k <- 10
#' linktype <- 'stow'
#'
#' traitnetoutput <- traitnet(traittypes = traittypes, allncats = allncats,
#' preftypes = preftypes, wvals = wvals, n=n, k=k, linktype = linktype, onecomp = TRUE,
#' visnet = TRUE)
#'
#' netmat <- traitnetoutput$net
#'
#'
#'# The same as the above example, but where the traits that determine the
#'# node colours and node sizes in the visualization are set by the user:
#'
#' traittypes <- c('cate', 'tnorm', 'ranks')
#' allncats <- c(5, NA, NA)
#' preftypes <- c('simclose', 'pophigh', 'simfar')
#' wvals <- c(0.5, 0.1, 0.3)
#' n <- 100
#' k <- 10
#' linktype <- 'stow'
#' vistraitnum <- c(3,1)
#'
#' traitnetoutput <- traitnet(traittypes = traittypes, allncats = allncats,
#' preftypes = preftypes, wvals = wvals, n=n, k=k, linktype = linktype, onecomp = TRUE,
#' visnet = TRUE, vistraitnum = vistraitnum)
#'
#' netmat <- traitnetoutput$net
#'
#'
#' # An example with two traits, where the trait values are provided by the user
#' # (the user-provided trait values here are simulated,
#' # but empirical trait data can be input in the same way):
#'
#' valuesfirsttrait <- sample(100)/100
#' valuessecondtrait <- sample(1:2,100, replace = T)
#'
#' traittypes <- c('own', 'own')
#' allncats <- c(NA, 2)
#' preftypes <- c('pophigh', 'simclose')
#' wvals <- c(0.1, 0.85)
#' n <- 100
#' k <- 10
#' linktype <- 'stow'
#' owntraitvals <- cbind(valuesfirsttrait, valuessecondtrait)
#' owntraitclasses <- c('cont', 'cate')
#'
#' traitnet(traittypes = traittypes, allncats = allncats,
#' preftypes = preftypes, wvals = wvals, n=100, k=10, linktype = linktype, onecomp = TRUE,
#' visnet = TRUE, owntraitclasses = owntraitclasses, owntraitvals = owntraitvals)








traitnet <- function(traittypes, allncats, preftypes, wvals, n, k, linktype, onecomp, visnet, vistraitnum = NULL, owntraitclasses = NULL, owntraitvals = NULL, fact = FALSE){

  if(length(n)==0 || length(k)==0 || length(traittypes)==0 || length(allncats)==0 || length(preftypes)==0 || length(linktype)==0 || length(wvals)==0 || length(onecomp)==0 || length(visnet)==0) {
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
    stop("Input lengths do not fit.")
  }
  if (!(linktype %in% c('stow','detw','unw')) ) {
    stop('linktype is incorrect.')
  }
  if (!all(traittypes %in% c('cate','tnorm','ranks','unif','own')) ) {
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
  if (!all(preftypes %in% c('simclose', 'simfar', 'pophigh', 'poplow')) ) {
    stop('preftypes incorrect.')
  }
  if (any(wvals[traittypes =='cate' & preftypes == 'simclose'] > 0.99) & onecomp == TRUE){
  warning(prompt = 'Using wval>0.99 is not recommended for simclose combined with a categorical trait and a single component, as it may be slow or impossible to generate the network.', immediate. = TRUE)
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

  traitvals <-  matrix(NA, nrow = n, ncol = length(traittypes))
  if (any(traittypes != 'own')){
    nonownplaces <- which(traittypes != 'own')
    for (curtraitnum in seq_along(nonownplaces)){
      curtraitplace <- nonownplaces[curtraitnum]
      traitvals[, curtraitplace] <- traitvalues(n=n, traittype = traittypes[curtraitplace], ncats = allncats[curtraitplace]) # make the traitvals for this trait and add them to the traitvals matrix
    }
  }
  if (any(traittypes == 'own')){
    ownplaces <- which(traittypes == 'own')
    traitvals[,ownplaces] <- owntraitvals
  }

  traitclassesraw <- traittypes
  traitclassesraw[traitclassesraw %in% c('tnorm', 'ranks', 'unif')] <- 'cont'  # set the relevant types from the types vector to cont
  if (any(traittypes == 'own')){ # add the classes from the classvector for traits that are inputted from the user (if any)
    ownplaces <- which(traittypes == 'own')
    traitclassesraw[ownplaces] <- owntraitclasses
  }
  traitclasses <- traitclassesraw

  sociatvals <- traitnetsociat(traitclasses = traitclasses, preftypes = preftypes, wvals = wvals, traitvals = traitvals)

  if (onecomp == FALSE){
    net <- traitnetbuild(n = n, k = k, linktype = linktype, sociatvals = sociatvals)
  } else if (onecomp == TRUE){ # if set to only one component, recreate networks until getting one with a single component
    componentnum <- 0
    while (round(componentnum) != 1) {  #to make sure the network is in one component
      net <- traitnetbuild(n = n, k = k, sociatvals = sociatvals, linktype = linktype)
      componentnum <- sna::components(net) # number of components of the current network. Note: igraph also has a function with this name.
    }
  }
  sociatmatraw <- matrix(0, n, n)
  inames <- rep(1:(n-1),(n-1):1)
  jnames <- sequence((n-1):1, from = 2:n)
  inamesboth <- c(inames, jnames)
  jnamesboth <- c(jnames, inames)
  bothboth <- cbind(inamesboth, jnamesboth)
  sociatmatraw[bothboth] <- sociatvals
  sociatmat <- sociatmatraw

  if (visnet == TRUE){
    traitnetvisual(net = net, traitclasses = traitclasses, preftypes = preftypes, linktype = linktype, wvals = wvals, traitvals = traitvals, vistraitnum = vistraitnum)
  }

  if (fact == TRUE){
    egg()
  }

  return(list(net = net, traitvals = traitvals, sociatmat = sociatmat))

}



