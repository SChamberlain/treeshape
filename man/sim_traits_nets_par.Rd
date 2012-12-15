\name{sim_traits_nets_par}
\alias{sim_traits_nets_par}
\title{Simulate networks, with interactions propoprtional to trait matching}
\usage{
  sim_traits_nets_par(listoftraitvecs, cores = 2,
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
trees <- rmtree(N=1000, n=60)
trees2 <- rmtree(N=1000, n=120)
registerDoMC(4)
traitvecs <- llply(trees, fastBM, bounds=c(0,Inf), .parallel=TRUE)
traitvecs2 <- llply(trees2, fastBM, bounds=c(0,Inf), .parallel=TRUE)
alltraits <- list(traitvecs, traitvecs2)
out<-sim_traits_nets_par(listoftraitvecs=alltraits, method="c", value=0.5, cores=4)
length(out)
}
}

