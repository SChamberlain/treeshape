\name{traitsig}
\alias{traitsig}
\title{Measure aspects of traits on trees.}
\usage{
  traitsig(x, y)
}
\arguments{
  \item{x}{A named trait vector, the names should match the
  tip names in the tree.}

  \item{y}{A phylogeny, has to be a phylo object.}
}
\description{
  Measure aspects of traits on trees.
}
\examples{
\dontrun{
trees <- rmtree(N=10, n=10)
traitvecs <- lapply(trees, fastBM)
traitsig(traitvecs, trees)
}
}

