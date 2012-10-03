\name{randtips}
\alias{randtips}
\title{Random traits - randomize traits across tips - shuffles names, but doesn't
		move position in tree of course}
\usage{
  randtips(w)
}
\arguments{
  \item{w}{A named vector of traits. Each element in vector
  must be named.}
}
\description{
  This is used in randomizing traits across the tips of a
  phylogeny so as to randomize any phylogenetic signal.
}
\examples{
\dontrun{
tree <- pbtree(n=10)
traitvec <- fastBM(tree)
randtips(traitvec)
}
}

