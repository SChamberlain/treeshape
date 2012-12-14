\name{sim_traits_nets2}
\alias{sim_traits_nets2}
\title{Simulate networks, with interactions propoprtional to trait matching}
\usage{
  sim_traits_nets2(listoftraitvecs,
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
trees <- rmtree(N=20, n=5)
trees2 <- rmtree(N=20, n=10)
traitvecs <- lapply(trees, fastBM, bounds=c(0,Inf))
traitvecs2 <- lapply(trees2, fastBM, bounds=c(0,Inf))
alltraits <- list(traitvecs, traitvecs2)
sim_traits_nets2(listoftraitvecs=alltraits, method="c", value=1)
}
}

