affinity_matrix <- function(D, k, alpha=1/6, beta=1/6)
{
  n = nrow(D)
  # make it symmetric in case it is not
  D = (D + t(D)) / 2
  diag(D) = 0
  D[D<0] = 0
  finiteMean = function(x){
    return(mean(x[is.finite(x)]))
  }
  # <Tricky> k nearest neighbors should exclude itself: 1:k+1
  # <?> here we divide by k+1, maybe k is better
  d = apply(apply(D, 2, sort)[1:k+1,], 2, finiteMean)
  sigma = alpha * (outer(d, d, '+')) + beta * D + .Machine$double.eps
  return(dnorm(D, 0, sigma))
}

kNN_graph = function(W, K) {
  n = nrow(W)
  # k-nearest neighbors does not include itself
  # set similarity = 0 for n-K-1 farest neighbors
  if (n - K - 1 < 1) {
    warning("K is too big: n-K-1 < 1. Assuming n > 2, 
          set K = floor(n/2)")
    K = floor(n / 2)
  }
  idx = t(apply(W, 1, order)[1:(n - K - 1),])
  for (i in seq_len(n)) {
    W[i, idx[i,]] = 0
  }
  # row normalization --> transition matrix
  S = W / rowSums(W)
  return(S)
}
