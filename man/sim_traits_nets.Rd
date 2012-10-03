\name{sim_traits_nets}
\alias{sim_traits_nets}
\title{Simulate networks, with interactions propoprtional to trait matching}
\usage{
  sim_traits_nets(listoftraitvecs,
    method = c("ratio", "complementarity", "barrier"),
    value)
}
\arguments{
  \item{listoftraitvecs}{}

  \item{method}{The model to be used to construct
  interaction matrices. One of
  "ratio","complementarity","barrier".}

  \item{value}{Value at which to determine if species
  interact or not - this value depends on the model you are
  using.}
}
\value{
  A data.frame of network structure metrics for balanced
  and unbalanced trees.
}
\description{
  Simulate networks, with interactions propoprtional to
  trait matching
}
\examples{
\dontrun{
trees <- rmtree(N=10, n=10)
trees2 <- rmtree(N=10, n=10)
traitvecs <- lapply(trees, fastBM)
traitvecs2 <- lapply(trees, fastBM)
alltraits <- list(traitvecs, traitvecs2)
sim_traits_nets(alltraits, "r", 1.5) # where r = ratio, you can abbreviate
}
}

