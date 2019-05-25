#' Root Mean Squared Error (RMSE)
#' 
#' Calculates the RMSE between simulated and experimental data, stored in a data frame.
#' This is a helper function for \code{\link{fit_quality}}.
#' @param data_frame A data frame that contains simulated and experimental data, one in 
#' each column.
#' @return Returns the RMSE value as double.
rmse <- function(df){ 
    a <- df$Data
    b <- df$Sim
    sqrt(sum((a-b)**2)/length(a))
} 