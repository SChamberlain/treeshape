\name{simbaltrees}
\alias{simbaltrees}
\title{Simulate a set of balanced and unbalanced trees.}
\usage{
  simbaltrees(tips_p = 10, metric, numtrees, cutlow,
    cuthigh, a, bounds, alpha, sigma, theta, alpha_eb,
    sigma_eb, rval = NULL, cval = NULL, asymm = 1,
    dumpmatrices = FALSE, matdir = "~",
    netmets = c("connectance", "links per species", "nestedness", "web asymmetry"))
}
\arguments{
  \item{tips_p}{Number of plant species to simulate in each
  tree - will be same for all trees. The number of animal
  species will be calculated based on the value of the
  asymm parameter.}

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

  \item{dumpmatrices}{If FALSE (default) matrices are not
  spit out to file, but if TRUE, then matrices are written
  to file in separate folders for each (logical)}

  \item{matdir}{Directory to output matrices to. Ignored if
  dumpmatrices=FALSE.}

  \item{netmets}{Network structure metrics to calculate -
  only use those that calculate single values for each
  matrix.}
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
netmets <- c("connectance", "nestedness","nodf2")
temp <- simbaltrees(tips_p=10, metric="colless", numtrees=5, cutlow=-0.5, cuthigh=0.5, a=10, bounds=c(0,100), alpha=1, sigma=1, alpha_eb=-0.8, sigma_eb=3, cval=0.8, asymm=2, dumpmatrices=TRUE, matdir="~/newfiles2", netmets=netmets)
head(temp[[1]]) # traits data.frame
head(temp[[2]]) # network structure output data.frame
}
}

