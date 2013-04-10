\name{getnetmets2}
\alias{getnetmets2}
\title{Calculate network metrics on matrices.}
\usage{
  getnetmets2(filepath,
    netmets = c("connectance", "NODF2"), directory)
}
\arguments{
  \item{filepath}{Path to an edgelist}

  \item{netmets}{Network structure metrics to calculate -
  only use those that calculate single values for each
  matrix.}
}
\description{
  This is used in calculating network metrics on pairs of
  matrices made from balanced and unbalanced trees.
}
\examples{
\dontrun{
getnetmets2(filepath = "/Users/scottmac2/testtest/sp30_bal_eb_twom_4.web")

filepaths <- c("/Users/scottmac2/testtest/sp30_bal_eb_twom_4.web",
								"/Users/scottmac2/testtest/sp30_bal_bm_twom_3.web",
								"/Users/scottmac2/testtest/sp30_bal_bm_twom_2.web",
								"/Users/scottmac2/testtest/sp30_bal_bm_twom_1.web")
l_ply(filepaths, getnetmets2)
}
}

