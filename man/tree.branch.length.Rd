\name{tree.branch.length}
\alias{tree.branch.length}
\title{Extracts the length of the branch above the given node. Returns 0 if the node is root.}
\usage{
  tree.branch.length(phylo, node)
}
\arguments{
  \item{phylo}{input phylo object}

  \item{node}{integer, the node's index}
}
\value{
  numeric, the branch length of the edge leading to the
  given node. May be NA.
}
\description{
  Extracts the length of the branch above the given node.
  Returns 0 if the node is root.
}

