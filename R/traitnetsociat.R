

#' Calculate social attraction values
#'
#' `traitnetsociat()` calculates social attraction values between individuals, based on social preferences for traits (attributes).
#' The social attraction can be based on one or multiple traits, which may be used with different preference types.
#'
#'
#'
#'
#'
#' The `traitnetsociat()` function calculates the social attraction for each pair of individuals, based on the two individuals'
#' trait values (attributes), the preference type set for each trait, and importance (weight) set for each preference-trait combination.
#' The calculation follows the trait preference model (see References).
#'
#' To construct networks based on the trait preference model, it is recommended to use the `traitnet()` and `traitnetsmetrics()` functions.
#' Those functions call this one, and there is thus no need to use this one directly.
#'
#' The function assumes that each input is ordered in the same way. For example, the first trait class given in the
#' `traitclasses` vector will be matched with the first preference type given in the `preftype` vector, the first
#' importance value given in the `wvals` vector, and the first column of trait values in the `traitvals` matrix. This also
#' means that the trait values for a given trait must be of the class given for that trait.
#'
#' About `traitclasses`:
#' * 'cate' indicates a categorical trait.
#' * 'cont' indicates a continuous trait.
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
#' About `wvals`:
#' * The sum of the `wval` vector, and each value in the vector, must lie in the interval \[0,1], where 0 indicates no importance and 1 indicates full importance.
#'
#' About `traitvals`:
#' * The trait values must fit with the trait classes given in the `traitclasses` argument.
#' Categorical trait values should consist of integers that each indicates a category.
#' Continuous trait values should consist of numbers between 0 and 1. Empirical trait values on other scales should be transformed to this scale before being used.
#'
#'
#'
#'
#' @param traitclasses A vector of strings indicating the trait class for each trait ('cate' or 'cont'; see Details)
#' @param preftypes A vector of strings indicating the preference type to be used with each trait ('simclose', 'simfar', 'pophigh' or 'poplow'; see Details)
#' @param wvals A vector with the importance (weight) of each trait-preference combination
#' @param traitvals A matrix with a column of trait values for each trait, and a row for each individual. Can be a vector if there is only a single trait
#'
#' @return A vector of social attraction values (one value for each pair of individuals). The values are ordered in the following way (with integers indicating names of individuals): 1&2, 1&3, ..., 1&n, 2&3, 2&4, ..., 2&n, ... n-1 and n.
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

#' @examples
#'
#' # An example of social attraction values that are based on the
#' # combined preferences for three traits of different type:
#'
#' # Create trait values based on some trait types:
#' traittypes <- c('cate', 'tnorm', 'ranks')
#' n <- 100
#' tval1 <- traitvalues(n = n, traittype = traittypes[1], ncats = 2)
#' tval2 <- traitvalues(n = n, traittype = traittypes[2])
#' tval3 <- traitvalues(n = n, traittype = traittypes[3])
#' traitvals <- cbind(tval1, tval2, tval3)
#'
#' # Set the class, preference type and importance for each trait
#' (The trait classes must fit with the trait values):
#' traitclasses <- c('cate', 'cont', 'cont')
#' preftypes <- c('simclose', 'pophigh', 'simfar')
#' wvals <- c(0.7, 0.1, 0.1)
#'
#' # Generate the social attraction values:
#' sociatvals <- traitnetsociat(traitclasses, preftypes, wvals, traitvals)
#'
#'
#'




traitnetsociat <- function(traitclasses, preftypes, wvals, traitvals){

  if(length(traitclasses)==0 || length(preftypes)==0 || length(wvals)==0 || length(traitvals)==0) {
    stop("Input without defaults must not be empty")
  }
  if(sum(wvals) > 1) {
    stop("The sum of wvals must not exceed 1.")
  }
  if(any(wvals <0)||any(wvals >1)) {
    stop("All entries in wvals must be between 0 and 1.")
  }
  if (!all(traitclasses %in% c('cate','cont')) ) {
    stop('traitclasses are incorrect.')
  }
  inputlengths <- c(length(traitclasses), length(preftypes), length(wvals), ncol(traitvals))
  if(any(inputlengths != length(traitclasses))) {
    stop("Inputs must have the same lengths.")
  }
  if (!all(preftypes %in% c('simclose', 'simfar','pophigh', 'poplow')) ) {
    stop('preftypes are incorrect.')
  }
  if(any(is.na(traitvals))) {
    stop("traitvals must not contain NA's.")
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

  simcloseprefs <- function(simtype, tsim){
    n <- length(tsim)
    inames <- rep(1:(n-1),(n-1):1)
    jnames <- sequence((n-1):1, from = 2:n)
    tsimalli <- tsim[inames]
    tsimallj <- tsim[jnames]
    if (simtype == 'cate') {
      psimclose  <- ifelse(tsimalli-tsimallj ==0, 1,0)
    } else if (simtype == 'cont'){
      psimclose <- 1 - abs(tsimalli-tsimallj)
    }
  }

  simfarprefs <- function(simtype, tsim){
    n <- length(tsim)
    inames <- rep(1:(n-1),(n-1):1)
    jnames <- sequence((n-1):1, from = 2:n)
    tsimalli <- tsim[inames]
    tsimallj <- tsim[jnames]
    if (simtype == 'cate') {
      psimfar  <- ifelse(tsimalli-tsimallj ==0, 0,1)
    } else if (simtype == 'cont'){
      psimfar <- abs(tsimalli-tsimallj)
    }
  }

  pophighprefs <- function(poptype, tpop){
    n <- length(tpop)
    inames <- rep(1:(n-1),(n-1):1)
    jnames <- sequence((n-1):1, from = 2:n)
    tpopalli <- tpop[inames]
    tpopallj <- tpop[jnames]
    if (poptype == 'cate') {
      npopcats <- length(unique(tpop))
      popcatvals <- seq(0,1,length.out = npopcats)
      ppophigh <- (popcatvals[tpopalli]+popcatvals[tpopallj])/2
    } else if (poptype == 'cont'){
      ppophigh <- (tpopalli+tpopallj)/2
    }
  }


  poplowprefs <- function(poptype, tpop){
    n <- length(tpop)
    inames <- rep(1:(n-1),(n-1):1)
    jnames <- sequence((n-1):1, from = 2:n)
    tpopalli <- tpop[inames]
    tpopallj <- tpop[jnames]
    if (poptype == 'cate') {
      npopcats <- length(unique(tpop))
      popcatvals <- seq(0,1,length.out = npopcats)
      ppoplow <- (2-(popcatvals[tpopalli]+popcatvals[tpopallj]))/2
    } else if (poptype == 'cont'){
      ppoplow <- (2-(tpopalli+tpopallj))/2
    }
  }

  # Linking probabilities for each dyad for each trait

  traitvals <- as.matrix(traitvals)                                                                  # All the trait values in matrix format
  n <- nrow(traitvals)                                                                               # number of individuals
  ndyads <- (n*n-n)/2                                                                                # number of dyads
  ptraits <- matrix(NA, nrow = ndyads, ncol = length(traitclasses))                                  #
  for (traitnum in seq_along(traitclasses)) {                                                        # for each trait
    if (preftypes[traitnum]=='simclose'){                                                                 # if the trait is of traittype 'simclose'
      ptraits[,traitnum] <- simcloseprefs(simtype = traitclasses[traitnum], tsim = traitvals[,traitnum])  # calculate linking probabilities (note: not yet weighted by importance) for the trait with the simcloseprefs functon
    } else if (preftypes[traitnum]=='simfar'){                                                                 # if the trait is of traittype 'simfar'
        ptraits[,traitnum] <- simfarprefs(simtype = traitclasses[traitnum], tsim = traitvals[,traitnum])  # calculate linking probabilities (note: not yet weighted by importance) for the trait with the simfarprefs functon
    } else if (preftypes[traitnum]=='pophigh'){                                                          # if the trait is of traittype 'pophigh'
      ptraits[,traitnum] <- pophighprefs(poptype = traitclasses[traitnum], tpop = traitvals[,traitnum])  # calculate linking probabilities for the trait with the pophighprefs functon
    } else if (preftypes[traitnum]=='poplow'){                                                          # if the trait is of traittype 'poplow'
      ptraits[,traitnum] <- poplowprefs(poptype = traitclasses[traitnum], tpop = traitvals[,traitnum])  # calculate linking probabilities (social attractions?) for the trait with the poplowprefs functon
      }
  }

  # Weigh the trait linking probabilities by importance
  ptws <- t(t(ptraits)*wvals)     # IS THIS THE SAME ORDER OF DYADS AS ABOVE trait linking probabilities weighted by importances

  # Calculate linking probabilities for randomness weighted by importance
  prand <- runif(ndyads)          # random number for each dyad
  wrand <- 1-sum(wvals)           # importance w of randomness
  prandtws <- wrand*prand         # random probability weighted by importance

  # Calculate final social attractions
  ptwsums <- rowSums(ptws)        # weighted trait linking probabilities summed for each dyad
  sociatvals <- ptwsums+prandtws  # social attraction values for each dyad (the above with random prob added)

  return(sociatvals)

}

