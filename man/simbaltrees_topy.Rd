\name{simbaltrees_topy}
\alias{simbaltrees_topy}
\title{Simulate a set of balanced and unbalanced trees.}
\usage{
  simbaltrees_topy(networksize = 10, metric, numtrees,
    cutlow, cuthigh, a, bounds, alpha, sigma, theta,
    alpha_eb, sigma_eb, rval = NULL, cval = NULL,
    asymm = 2.47, cores = 2, dumpmatrices = TRUE,
    matdir = "~", modeltorun = "complementarity",
    output = "matrix", includetraitvar = FALSE)
}
\arguments{
  \item{networksize}{Number of total species to simulate in
  each tree - will be same for all trees. The number of
  plant and animal species will be calculated based on the
  value of the asymm parameter.}

  \item{metric}{Methods to use to generate trees, one of
  "colless", "beta", or gamma (see details). Defaults to
  "colless".}

  \item{numtrees}{Number of trees to produce. Defaults to
  10 trees.}

  \item{cutlow}{Value at which to filter trees on the low
  (e.g., unbalanced) side of the metric.}

  \item{cuthigh}{Value at which to filter trees on the high
  (e.g., balanced) side of the metric.}

  \item{a}{If model = "bm", a value for ancestral state at
  the root node.}

  \item{bounds}{If model = "bm", a vector with the lower
  and upper bounds (respectively) for bounded Brownian
  simulation - by default simulation is unbounded.}

  \item{alpha}{If model = "ou", a value for the
  Ornstein-Uhlenbeck model of trait evolution.  From ape
  documentation: "a numeric vector giving the strength of
  the selective constraint for each branch (can be a single
  value)."}

  \item{sigma}{If model = "ou", is the single value of the
  standard-deviation of the random component for each
  branch (can be a single value).}

  \item{theta}{If model = "ou", a numeric vector giving the
  optimum for each branch (can be a single value)}

  \item{alpha_eb}{If model = "eb", is the exponent of the
  relationship between rate and time in the
  exponentialchange model.}

  \item{sigma_eb}{If model = "eb", the single value of the
  standard-deviation of the random component for each
  branch (can be a single value).  This sigma means the
  same as for the OU model, but this allows to specify it
  separately.}

  \item{rval}{Value for the ratio model for constructing
  interaction network matrices.}

  \item{cval}{Value for the complimentarity model for
  constructing interation network matrices.}

  \item{asymm}{An asymmetry value in terms of ratio of
  animals:plants instead of the traditional definition like
  (Na-Np/Na+Np), where Na is number of animal species, and
  Np is number of plant species. For example, asymm=2 would
  mean twice as many animal as plant species.}

  \item{cores}{Number of cores to use in mcmapply in
  simulating networks form traits.}

  \item{dumpmatrices}{If FALSE (default) matrices are not
  spit out to file, but if TRUE, then matrices are written
  to file in separate folders for each (logical)}

  \item{matdir}{Directory to output matrices to. Ignored if
  dumpmatrices=FALSE.}

  \item{netmets}{Network structure metrics to calculate -
  only use those that calculate single values for each
  matrix.}

  \item{output}{One of matrix or edgelist.}
}
\value{
  A data.frame of network structure metrics for balanced
  and unbalanced trees.
}
\description{
  Could use this in asking questions aobut how phylogenetic
  tree balance influences ____ (madlib it).
}
\examples{
\dontrun{
temp <- simbaltrees_topy(networksize=30, metric="colless", numtrees=5, cutlow=-0.7, cuthigh=0.7, a=10, bounds=c(0,100), alpha=10, sigma=3, alpha_eb=-1.1, sigma_eb=3, cval=0.5, asymm=2, matdir="~/testtest", modeltorun="complementarity", output="edgelist", includetraitvar = TRUE)
head(temp) # traits data.frame
}
}

