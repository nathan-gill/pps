
find_duplicate <- function(p) {
  if (length(unique(p)) != length(p)) {
    return(1)
  } else {
    return(0)
  }
}

#function to flip the off-diagonals of a matrix

flip_off_diag <- function(A) {
  n <- dim(A)[1]
  for (k in 1:n) {
    for (l in 1:n) {
      if (k != l) {
        A[k,l] <- -A[k,l]
      }
    }
  }
  return(A)
}



#function to compute paths up to a certain length
#between two nodes in an igraph object

#K is number of edges

paths_up_to_k <- function(gr, i, j, K) {
  paths <- purrr::imap(as.vector(igraph::neighbors(gr, i)), ~c(i, .x))

  valid_paths <- list()
  #remove paths that have duplicates
  inds <- which(purrr::map_dbl(paths, find_duplicate) == 1)
  if (length(inds) != 0) {
    paths <- paths[-inds]
  }

  #add paths ending in j to valid_paths
  inds <- which(purrr::map_dbl(paths, ~ifelse(dplyr::last(.x) == j, 1, 0)) == 1)
  valid_paths <- append(valid_paths, paths[inds])
  if (length(inds) != 0) {
    paths <- paths[-inds]
  }



  for (l in 1:(K-1)) {

    if (length(paths) == 0) {
      return(valid_paths)
    }
    #add all paths one node from existing
    temp <- list()
    for (m in 1:length(paths)) {
      temp <- append(temp, purrr::map(as.vector(igraph::neighbors(gr, dplyr::last(paths[[m]]))), ~c(paths[[m]], .x)))
    }
    paths <- temp


    #remove paths that have duplicates
    inds <- which(purrr::map_dbl(paths, find_duplicate) == 1)
    if (length(inds) != 0) {
      paths <- paths[-inds]
    }

    #add paths ending in j to valid_paths
    inds <- which(purrr::map_dbl(paths, ~ifelse(dplyr::last(.x) == j, 1, 0)) == 1)
    valid_paths <- append(valid_paths, paths[inds])
    if (length(inds) != 0) {
      paths <- paths[-inds]
    }
  }
  return(valid_paths)
}


#function to compute gamma for every path between two
#nodes in an igraph object

get_path_contribution_k <- function(i, j, pcor, K) {
  if (i == j) {
    return(1)
  }
  n <- dim(pcor)[1]
  #adjust signs of partial correlation matrix

  A <- pcor
  for (k in 1:n) {
    for (l in 1:n) {
      if (k != l) {
        A[k,l] <- -A[k,l]
      }
    }
  }

  #convert A to an adjacency matrix

  G <- ifelse(A != 0, 1, 0) - diag(rep(1, n))

  #get igraph object

  gr <- igraph::graph_from_adjacency_matrix(G, mode = "undirected")

  #get all paths

  paths <- paths_up_to_k(gr, i, j, K)

  n_paths <- length(paths)
  if (n_paths == 0) {
    return(0)
  }
  cor_val <- 0
  path_list <- list()
  cont_list <- list()
  weight_list <- list()
  for (k in 1:length(paths)) {
    path <- as.vector(paths[[k]])
    len_path <- length(path)

    #get sign of path
    sign <- ifelse(length(path) %% 2 == 0, -1, 1)

    #get product of partial correlations along the path
    prod <- 1
    for (l in 1:(len_path-1)) {
      prod <- prod*A[[path[l], path[l+1]]]
    }

    #get determinant of the relevant submatrix
    if (len_path == n-1) {
      d <- A[-path, -path]
    } else if (len_path == n) {
      d <- 1
    } else {
      d <- det(A[-path, -path])
    }

    cor_val <- cor_val + sign*prod*d
    cont_list <- append(cont_list, sign*prod*d*(det(A[-i, -i])*det(A[-j, -j]))^(-1/2))
    path_list <- append(path_list, list(path))
    weight_list <- append(weight_list, d)#*sign*(det(A[-i, -i])*det(A[-j, -j]))^(-1/2))
  }

  #normalize

  cor_val <- cor_val*(det(A[-i, -i])*det(A[-j, -j]))^(-1/2)
  return(list(cont_list, path_list, weight_list))
}
