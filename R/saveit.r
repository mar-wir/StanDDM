#' Save R Objects With External Name
#' 
#' A wrapper function around \code{\link[base]{save}}. Assign a string as filename for a R object, 
#' instead the environment name. Will save file to disk.
#' @export
#' @param file Any R object in the environment of any type or class.
#' @param string A string used as filename for the saved R object
saveit <- function(..., string, file) {
    
    x <- list(...)
    names(x) <- string
    save(list=names(x), file=file, envir=list2env(x))
}