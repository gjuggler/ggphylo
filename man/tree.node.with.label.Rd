\name{tree.node.with.label}
\alias{tree.node.with.label}
\title{Returns the index of the node with a given label.}
\usage{
  tree.node.with.label(phylo, label, return.one = FALSE)
}
\arguments{
  \item{phylo}{input phylo object}

  \item{label}{character, the label to search for in the
  tree.}

  \item{return.one}{boolean, TRUE if a maximum of one
  matching node should be returned, FALSE if all matching
  nodes should be returned. Defaults to TRUE.}
}
\value{
  integer vector corresponding to the indices of all nodes
  with the given label. Returns a zero-length vector if no
  nodes matched.
}
\description{
  Returns the index of the node with a given label.
}

