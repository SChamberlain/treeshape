\name{simbal_plot}
\alias{simbal_plot}
\title{Plotting output from simbaltrees}
\usage{
  simbal_plot(output, netval, justdf = FALSE)
}
\arguments{
  \item{output}{Output from simbaltrees.}

  \item{netval}{A network metric used in your call to
  simbaltrees.}

  \item{justdf}{Return just the data.frame or not (default
  to FALSE).}
}
\description{
  Plotting output from simbaltrees
}
\examples{
\dontrun{
netmets <- c("connectance", "nestedness","nodf2")
outout <- suppressMessages(llply(c(10, 15), function(x)
		simbaltrees(tips_p=10, metric="colless", numtrees=5, cutlow=-0.5, cuthigh=0.5,
		a=10, bounds=c(0,100), alpha=1, sigma=1, alpha_eb=-0.8, sigma_eb=3, cval=0.8,
		asymm=2, dumpmatrices=TRUE, matdir="~/newfiles2", netmets=netmets)
simbal_plot(outout, "connectance")
}
}

