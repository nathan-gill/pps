pps_no_order <- function (P, i, j, K = 5, prec = TRUE, use.names = TRUE)
{
  if (!is.matrix(P)) {
    stop("Invalid precision or partial correlation matrix.")
  }
  else if (dim(P)[1] != dim(P)[2]) {
    stop("P must be a square matrix.")
  }
  else if (!prec) {
    if (sum(diag(P)) != dim(P)[1]) {
      stop("Partial correlation matrix must have 1s on the diagonal.")
    }
  }
  names <- colnames(P)
  if (is.character(i)) {
    if (i %in% names) {
      i <- which(names == i)
    }
    else {
      stop("Invalid node name.")
    }
  }
  if (is.character(j)) {
    if (j %in% names) {
      j <- which(names == j)
    }
    else {
      stop("Invalid node name.")
    }
  }
  if ((i > dim(P)[1]) | (j > dim(P)[1])) {
    stop("Index exceeds total number of nodes.")
  }
  if (prec == TRUE) {
    P <- flip_off_diag(cov2cor(P))
  }
  a <- get_path_contribution_k(i, j, P, K)
  if (is.double(a)) {
    return(NULL)
  }
  #ind <- order(abs(unlist(a[[1]])), decreasing = TRUE)
  ind <- 1:length(a[[1]])
  path.list <- list()
  if ((use.names == TRUE) && (!is.null(colnames(P)))) {
    for (k in 1:length(a[[2]])) {
      path.list <- append(path.list, list(names[unlist(a[[2]][[ind[k]]])]))
    }
  }
  else {
    for (k in 1:length(a[[2]])) {
      path.list <- append(path.list, list(unlist(a[[2]][[ind[k]]])))
    }
  }
  pps.list <- purrr::map(a[[1]][ind], ~ifelse(is.na(.x), 0,
                                              abs(.x))/sum(abs(unlist(a[[1]]))))
  return(list(path = path.list, pps = pps.list, gamma = a[[1]][ind]))
}


#Function to score two results
is.equal <- function(p1, p2) {
  if (!(length(p1) == length(p2))) {
    return(0)
  } else {
    s <- 0
    for (i in 1:length(p1)) {
      s <- s + (p1[i] == p2[i])
    }
    if (s == length(p1)) {
      return(1)
    } else {
      return(0)
    }
  }
}


get_score <- function(res1, res2, k) {
  #get top 5
  res1 <- res1$path[1:k]
  res2 <- res2$path[1:k]
  #for each in res1, check to see if in res2
  score <- 0
  for (i in 1:k) {
    p <- res1[[i]] %>% unlist()
    for (j in 1:k) {
     if (is.equal(p, unlist(res2[[j]])) == 1) {
       score <- score + 1
       if (i == j) {
         score <- score + 1
       }
     }
    }
  }
  return(score)
}

#extract length of acylcarnitine

get_length <- function(ac_string) {

}


get_edges <- function(res, names) {
  edge_matrix <- diag(rep(1, length(names)))
  if (!is.null(res)) {
  path <- res$path[[1]]
  n <- length(path)
  path_inds <- rep(0, n)
  for (i in 1:n) {
    path_inds[i] <- which(names == path[i])
  }
  for (i in 1:(n-1)) {
    edge_matrix[path_inds[i], path_inds[i+1]] <- 1
    edge_matrix[path_inds[i+1], path_inds[i]] <- 1
  }
  return(edge_matrix)
  } else {
    return(NULL)
  }
}

