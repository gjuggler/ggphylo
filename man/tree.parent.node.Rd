\name{tree.parent.node}
\alias{tree.parent.node}
\title{Extracts the parent node ID for the given node. Returns -1 if the node is root.
Return the index of the node directly above the given node. Returns -1 if the given node is root.}
\usage{
  tree.parent.node(phylo, node)
}
\arguments{
  \item{phylo,}{input phylo object}

  \item{node,}{integer index of the node whose parent is
  desired}
}
\value{
  integer, the index of the parent node or -1 if the given
  node is root.
}
\description{
  Extracts the parent node ID for the given node. Returns
  -1 if the node is root. Return the index of the node
  directly above the given node. Returns -1 if the given
  node is root.
}

