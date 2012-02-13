\name{is.leaf.phylo}
\alias{is.leaf.phylo}
\title{Determine if the given node is a leaf or an internal node. Alias of \code{\link{tree.is.leaf}}.}
\usage{
  \method{is.leaf}{phylo} (phylo, node)
}
\arguments{
  \item{phylo,}{input phylo object}

  \item{node,}{integer index of the desired node to test}
}
\value{
  boolean, TRUE if the given node is a leaf, FALSE if it is
  an internal node
}
\description{
  Determine if the given node is a leaf or an internal
  node. Alias of \code{\link{tree.is.leaf}}.
}

