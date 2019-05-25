#' Random Value Generation: Uniform Distribution Scale-Location Reparametrization
#' 
#' A long title for a small function. Helper function for \code{\link{drift_simuls}} which reparametrizes the arguments of
#' the \code{\link{runif}} function as scale and location parameters.
#' 
#' @export
#' @param n Amount of values to be generated.
#' @param loc Location parameter.
#' @param scale Scale parameter.
#' @return A vector of doubles of length "n".
runif2 <- function(n, loc, scale){ # wrap runif in location scale variant
    out <- loc + scale * runif(n)
    as.vector(out)
}