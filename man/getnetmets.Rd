\name{getnetmets}
\alias{getnetmets}
\title{Calculate network metrics on matrices.}
\usage{
  getnetmets(balanced, unbalanced, netmets)
}
\arguments{
  \item{balanced}{A list of matrices made from balanced
  trees.}

  \item{unbalanced}{A list of matrices made from unbalanced
  trees.}

  \item{netmets}{Network structure metrics to calculate -
  only use those that calculate single values for each
  matrix.}
}
\description{
  This is used in calculating network metrics on pairs of
  matrices made from balanced and unbalanced trees.
}
\examples{
\dontrun{
# Let's pretend these are balanced and unbalanced matrices
m <- 10
n <- 5
netmets <- c("connectance", "links per species", "nodf2")
mats_rand_bal <- replicate(20, matrix(rbinom(m * n, 1, .5), ncol = m, nrow = n), FALSE)
mats_rand_unbal <- replicate(20, matrix(rbinom(m * n, 1, .5), ncol = m, nrow = n), FALSE)
getnetmets(balanced = mats_rand_bal, unbalanced=mats_rand_unbal, netmets)
}
}

