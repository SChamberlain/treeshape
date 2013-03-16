\name{sim_traits_nets_par}
\alias{sim_traits_nets_par}
\title{Simulate networks, with interactions propoprtional to trait matching}
\usage{
  sim_traits_nets_par(listoftraitvecs, type = NULL,
    traitm = NULL, matdir = "~",
    method = c("complementarity", "barrier", "twomethods"),
    includetraitvar = FALSE, value = NULL,
    output = "matrix")
}
\arguments{
  \item{listoftraitvecs}{Nested list of trait vectors.}

  \item{type}{Balanced ("bal") or unbalanced ("unbal").}

  \item{traitm}{Evolutionary model, used to create file
  name, one of "bm", "ou", or "eb".}

  \item{matdir}{directory to put matrices in}

  \item{method}{The model to be used to construct
  interaction matrices. One of "complementarity",
  "barrier", or "twomethods" (for combined
  complementarity/barrier).}

  \item{includetraitvar}{Include trait variation in the
  simulation or not.  Defaults to FALSE.}

  \item{value}{Value at which to determine if species
  interact or not - this value depends on the model you are
  using.}

  \item{output}{One of matrix or edgelist.}
}
\value{
  A data.frame of network structure metrics for balanced
  and unbalanced trees (if output="matrix"), or an edgelist
  (if output="edgelist").
}
\description{
  Simulate networks, with interactions propoprtional to
  trait matching
}
\examples{
\dontrun{
trees <- rmtree(N=100, n=15)
trees2 <- rmtree(N=100, n=30)
library(doMC)
registerDoMC(4)
traitvecs <- llply(trees, fastBM, bounds=c(0,Inf), .parallel=TRUE)
traitvecs2 <- llply(trees2, fastBM, bounds=c(0,Inf), .parallel=TRUE)
alltraits <- list(traitvecs, traitvecs2)
sim_traits_nets_par(listoftraitvecs=alltraits, method="c", value=0.5, type="bal", traitm="bm", matdir="~/testtest", includetraitvar = TRUE)
sim_traits_nets_par(listoftraitvecs=alltraits, method="b", type="bal", traitm="bm", matdir="~/newfiles2")
length(out)
}
}

