\name{tree.node.with.label}
\alias{tree.node.with.label}
\title{Returns the index of the node with a given label.}
\usage{
  tree.node.with.label(phylo, label)
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
  Returns the index of the node with a given label.
}

