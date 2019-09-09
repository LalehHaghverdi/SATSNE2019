cosine.norm <- function(X)
  # Computes the cosine norm, with some protection from zero-length norms.
  # Cell rows and Gene columns
{
  cellnorm <- pmax(1e-8, sqrt(colSums(X^2)))
  X/matrix(cellnorm, nrow(X), ncol(X), byrow=TRUE)
}