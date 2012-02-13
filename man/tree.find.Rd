\name{tree.find}
\alias{tree.find}
\title{Returns the index of the node with a given label. Alias for \code{\link{tree.node.with.label}}.}
\usage{
  tree.find(...)
}
\arguments{
  \item{phylo}{input phylo object}

  \item{label}{character, the label to search for in the
  tree.}
}
\value{
  integer vector corresponding to the indices of all nodes
  with the given label. Returns a zero-length vector if no
  nodes matched.
}
\description{
  Returns the index of the node with a given label. Alias
  for \code{\link{tree.node.with.label}}.
}

