#' Randomly sample cells (columns)
#'
#' This function takes a matrix or data frame as
#' input and randomly samples columns
#'
#' @param df Data frame or matrix.
#' @param n Number of cells to sample.
#'
#' @return Matrix or data frame of n randomly sampled rows
#'
#' @export
#'
#' @examples \dontrun{
#' df <- matrix(1:1000, 500, 20)
#' sub_df <- random_cells(
#'     df = df,
#'     n = 20
#')
#'}

random_cells <- function(df, n){
  return(df[, sample(ncol(df), n)])
}
