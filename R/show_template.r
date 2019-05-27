#' Model Template
#' 
#' Displays a template for creating custom models. They can be added to the "include_models"
#' argument in the \code{\link{StanDDM}} function (one at the time). They need first to be loaded in
#' memory, then, the name has to be included in the argument.
#' @export
#' @return Displays a model template in text form, to be copied and modified.
show_template <- function(){
    file.show('model_function_template.r')
}