
#' plot co clustered data
#'
#' This function plots data from a lbvem object. Comparing side by side raw data and co clustered data
#'
#' @param X Data matrix
#' @param model lbvem object
#'
#' @return nothing
#'
#' @import graphics

#' @examples
#' #library(blockcluster)
#' #data(gaussiandata)
#' #model <- lbvem(gaussiandata, g=3, m=4, niter=5)
#' #plotcoclust(gaussiandata, model)
#' @export
plotcoclust <- function(X, model){
  order_row <- order( apply(model$mu, 1, sum) )
  order_col <- order( apply(model$mu, 2, sum) )

  row_posterior <- apply(model$Z_tild, 1, which.max)
  col_posterior <- apply(model$W_tild, 1, which.max)

  ind_row <- order( factor(row_posterior, levels=order_row) )
  ind_col <- order( factor(col_posterior, levels=order_col) )
  X_bloc <- X[ind_row, ind_col]
  par(mfrow=c(1,2))
  image(X)
  title("raw data")
  image(t(X_bloc))
  abline()
  title("co-clustered data")
}
