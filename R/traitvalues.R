

#' Generate trait values
#'
#' `traitvalues()` generates a vector of trait values based on a given trait type.
#'
#'
#' To generate trait values, the function draws values from a distribution that corresponds to the given trait type.
#'
#' To make networks based on trait preferences, it is recommended to use the `traitnet()` or `traitnetsmetrics()` functions.
#' Those functions call this function and there is no need to use this function directly.
#'
#' About `traittype`:
#' * 'cate'  gives categorical trait values. The values consist of the integers from 1 to `ncats`, which each indicates a category. Each category is drawn with the same probability.
#' * 'tnorm' gives normally distributed trait values (truncated). The values are drawn from a truncated normal distribution that lie in the interval \[0,1], with mean = 0.5 and standard deviation = 0.25.
#' * 'ranks' gives trait values indicating ranks (e.g. social rank). The values are homogeneously distributed in the interval \[0,1]. This corresponds to ranks indicated by integers from 1 to `n`, transformed to the interval \[0,1].
#' * 'unif'  gives uniformly distributed trait values. The values are drawn from a uniform distribution with the interval \[0,1].

#' @param n An integer indicating the number of trait values to be made (i.e. number of individuals for which the trait values are made)
#' @param traittype A string indicating the trait type ('cate', 'tnorm', 'ranks', or 'unif'; see Details)
#' @param ncats An integer indicating the number of categories (only used if `traittype` is set to categorical)
#'
#' @return A vector of trait values. The vector has length `n`. The values are in random order. Non-categorical values are in the interval \[0,1]; categorical values are given by integers.

#' @export
#'
#' @references
#' Brask, J. B., Brent, L. J., & De Moor, D. (2025). okaapi: an R package for generating social networks based on trait preferences. EcoEvoRxiv.
#' \url{https://ecoevorxiv.org/repository/view/9622/}
#'
#' Brask, J. B., Koher, A., Croft, D. P., & Lehmann, S. (2026). Far-reaching consequences of trait preferences for animal social network structure and function. Behavioral Ecology, 37(1), araf132.
#' \url{https://academic.oup.com/beheco/article/37/1/araf132/8339509}
#'
#' Brask, J. B., Silk, M., & Weiss, M. N. (2025). An introduction to generative network models and their use in animal sociality research. Animal Behaviour, 123364.
#' \url{https://www.sciencedirect.com/science/article/pii/S000334722500291X}
#'
#'
#' @author
#' Josefine Bohr Brask, \email{bohrbrask@@gmail.com}
#'
#' @seealso \code{\link{traitnet()}}
#' @seealso \code{\link{traitnetsmetrics()}}
#'
#' @examples
#' # Categorical trait values with two categories, for 100 individuals:
#' traitvalues(n = 100, traittype = 'cate', ncats = 2)
#'
#' @examples
#' # Continuous trait values from a (truncated) normal distribution, for 100 individuals:
#' traitvalues(n = 100, traittype = 'tnorm')
#'
#' @examples
#' # Rank trait values, for 50 individuals:
#' traitvalues(n = 50, traittype = 'ranks')
#'
#' @examples
#' # Continuous uniform trait values, for 50 individuals:
#' traitvalues(n = 50, traittype = 'unif')


traitvalues <- function(n, traittype, ncats) {
  if (n <= 0 || is.na(n)) {
    stop('n must be positive.')
  }
  if (traittype == 'cate') {
    if (ncats <= 0 || length(ncats)==0 || is.na(ncats)) {
      stop('ncats must be a positive integer.')
    }
    traitvals <- sample(ncats,n, replace = TRUE)
  } else if (traittype == 'tnorm') {
    traitvals <- truncnorm::rtruncnorm(n, a=0,b=1,mean=0.5, 0.25)
  } else if (traittype == 'ranks') {
    traitvals <- sample(n,n, replace = FALSE)/n
  } else if (traittype == 'unif') {
    traitvals <- runif(n)
  }
  return(traitvals)
}

