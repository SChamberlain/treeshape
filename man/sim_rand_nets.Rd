\name{sim_rand_nets}
\alias{sim_rand_nets}
\title{Measure aspects of traits on trees.}
\usage{
  sim_rand_nets(listoftrees)
}
\arguments{
  \item{listoftrees}{}
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

