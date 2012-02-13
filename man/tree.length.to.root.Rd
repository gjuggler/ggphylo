\name{tree.length.to.root}
\alias{tree.length.to.root}
\title{Returns the length from the tree root to the given node. The input
 node can either be input as a node index or a node label.}
\usage{
  tree.length.to.root(phylo, node)
}
\arguments{
  \item{phylo}{input phylo object}

  \item{node}{integer or character. When integer, the node
  index; when character, the node label}
}
\value{
  numeric, the total branch length separating the tree root
  and the given node.
}
\description{
  Returns the length from the tree root to the given node.
  The input node can either be input as a node index or a
  node label.
}

