\name{tree.depth.to.root}
\alias{tree.depth.to.root}
\title{Returns the length from the tree root to the given node. Note: the node
can be given either as a node ID or a tip label.}
\usage{
  tree.depth.to.root(phylo, node)
}
\arguments{
  \item{phylo}{input phylo object}

  \item{node,}{either integer or character. When integer,
  the node index of the desired node to query. When
  character, the label of the desired node to query.}
}
\description{
  Returns the length from the tree root to the given node.
  Note: the node can be given either as a node ID or a tip
  label.
}

