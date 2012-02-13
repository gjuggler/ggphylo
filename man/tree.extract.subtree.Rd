\name{tree.extract.subtree}
\alias{tree.extract.subtree}
\title{Returns the minimum sub-tree connecting the given leaf nodes}
\usage{
  tree.extract.subtree(phylo, leaves)
}
\arguments{
  \item{phylo,}{input phylo object}

  \item{leaves,}{integer vector containing the node indices
  of the leaves from which to extract a sub-tree}
}
\value{
  phylo, a phylo object containing the minimum spanning
  sub-tree defined by the given leaf nodes
}
\description{
  Returns the minimum sub-tree connecting the given leaf
  nodes
}

