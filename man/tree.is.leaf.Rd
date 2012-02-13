\name{tree.is.leaf}
\alias{tree.is.leaf}
\title{Determine if the given node is a leaf or an internal node.}
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
  node.
}

