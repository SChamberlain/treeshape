\name{simbaltrees}
\alias{simbaltrees}
\title{Simulate a set of balanced and unbalanced trees.}
\usage{
  simbaltrees(tips = 10, metric, numtrees, cutlow, cuthigh,
    alpha, sigma,
    netmets = c("connectance", "links per species", "nestedness", "web asymmetry"))
}
\arguments{
  \item{tips}{Number of species to simulate in each tree -
  will be same for all trees.}

  \item{metric}{Methods to use to generate trees, one of
  "colless", "beta", or gamma (see details). Defaults to
  "colless".}

  \item{numtrees}{Number of trees to produce. Defaults to
  10 trees.}

  \item{cutlow}{Value at which to filter trees on the low
  (e.g., unbalanced) side of the metric.}

  \item{cuthigh}{Value at which to filter trees on the high
  (e.g., balanced) side of the metric.}

  \item{alpha}{Alpha value for the Ornstein-Uhlenbeck model
  of trait evolution.  From ape documentation: "a numeric
  vector giving the strength of the selective constraint
  for each branch (can be a single value)."}

  \item{sigma}{Sigma is the single value of the
  standard-deviation of the random component for each
  branch (can be a single value).}

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
netmets <- c("connectance", "links per species", "nestedness", "web asymmetry")
temp <- simbaltrees(tips=10, metric="colless", numtrees=10, cutlow=-0.5, cuthigh=0.5, alpha=0.1, sigma=1, netmets=netmets)
unique(temp$model)
}
}

