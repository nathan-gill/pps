#' Get list of paths between two nodes and their PPS
#'
#' For a given precision matrix or partial correlation matrix and a given
#' pair of nodes, computes all paths between those nodes up to a certain length
#' and their PPS
#'
#'
#' @param i,j Node numbers. Should be positive integers.
#' @param P A precision matrix or partial correlation matrix.
#' @param K The maximum path length to search up to. Defaults to 5.
#' @param prec If \code{TRUE}, indicates that \code{P} is a precision matrix.
#'             If \code{FALSE}, \code{P} is taken to be a partial correlation matrix.
#'             Defaults to TRUE.
#' @param use.names Whether or not to ouput paths using node names. Defaults to \code{TRUE}.
#'                  If the precision/partial correlation matrix does not have column names,
#'                  indices will be used.
#'
#' @return A list with the following elements:
#'   \item{path}{A list of paths between nodes \code{i} and \code{j} in order of descending PPS.}
#'   \item{pps}{A list of the PPS values for the paths in \code{path}.}
#'   \item{gamma}{The individual contributions of each path in \code{path} on the correlation scale.}
#'
#' @export


pps <- function(P, i, j, K = 5, prec = TRUE, use.names = TRUE) {

  #check if P is a valid precision or partial correlation matrix

  if (!is.matrix(P)) {
    stop("Invalid precision or partial correlation matrix.")
  } else if (dim(P)[1] != dim(P)[2]) {
    stop("P must be a square matrix.")
  } else if (!prec) {
    if (sum(diag(P)) != dim(P)[1]) {
      stop("Partial correlation matrix must have 1s on the diagonal.")
    }
  }

  names <- colnames(P)

  #check that i and j are valid
  if (is.character(i)) {
    if (i %in% names) {
      i <- which(names == i)
    } else {
      stop("Invalid node name.")
    }
  }

  if (is.character(j)) {
    if (j %in% names) {
      j <- which(names == j)
    } else {
      stop("Invalid node name.")
    }
  }

  if ((i > dim(P)[1]) | (j > dim(P)[1])) {
    stop("Index exceeds total number of nodes.")
  }

  #convert to partial correlation matrix
  if (prec == TRUE) {
    P <- flip_off_diag(cov2cor(P))
  }

  a <- get_path_contribution_k(i,j, P, K)

  #check to see if there are no paths between the nodes
  if (is.double(a)) {
    return(NULL)
  }

  #order by decreasing gamma magnitude
  ind <- order(abs(unlist(a[[1]])), decreasing = TRUE)

  #get list of paths and pps
  path.list <- list()
  if ((use.names == TRUE) && (!is.null(colnames(P)))) {
    for (k in 1:length(a[[2]])) {
      path.list <- append(path.list, list(names[unlist(a[[2]][[ind[k]]])]))
    }
  } else {
      for (k in 1:length(a[[2]])) {
        path.list <- append(path.list, list(unlist(a[[2]][[ind[k]]])))
      }
  }

  pps.list <- purrr::map(a[[1]][ind], ~ifelse(is.na(.x), 0, abs(.x))/sum(abs(unlist(a[[1]]))))

  return(list(path = path.list, pps = pps.list, gamma = a[[1]][ind]))
}

