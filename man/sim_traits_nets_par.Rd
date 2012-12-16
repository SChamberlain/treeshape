\name{sim_traits_nets_par}
\alias{sim_traits_nets_par}
\title{Simulate networks, with interactions propoprtional to trait matching}
\usage{
  sim_traits_nets_par(listoftraitvecs, type = NULL,
    traitm = NULL, matdir = "~",
    method = c("ratio", "complementarity", "barrier"),
    value)
}
\arguments{
  \item{listoftraitvecs}{}

  \item{type}{}

  \item{traitm}{}

  \item{matdir}{}

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
trees <- rmtree(N=100, n=15)
trees2 <- rmtree(N=100, n=30)
registerDoMC(4)
traitvecs <- llply(trees, fastBM, bounds=c(0,Inf), .parallel=TRUE)
traitvecs2 <- llply(trees2, fastBM, bounds=c(0,Inf), .parallel=TRUE)
alltraits <- list(traitvecs, traitvecs2)
sim_traits_nets_par(listoftraitvecs=alltraits, method="c", value=0.2, type="bal", traitm="bm", matdir="~/newfiles2")
length(out)
}
}

