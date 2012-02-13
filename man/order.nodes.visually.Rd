\name{order.nodes.visually}
\alias{order.nodes.visually}
\title{Returns the indices of all nodes in the tree, arranged in the order
in which they would be displayed in a 2-dimensional phylogram plot. In other words,
the parent node of two leaf nodes will be placed directly in between its two children,
since the parent node is usually placed at the midpoint between two leaf nodes in a plot.}
\usage{
  order.nodes.visually(phylo)
}
\arguments{
  \item{phylo,}{input phylo object}
}
\value{
  integer vector, the indices of all nodes in the tree,
  arranged in "visual" order.
}
\description{
  Returns the indices of all nodes in the tree, arranged in
  the order in which they would be displayed in a
  2-dimensional phylogram plot. In other words, the parent
  node of two leaf nodes will be placed directly in between
  its two children, since the parent node is usually placed
  at the midpoint between two leaf nodes in a plot.
}

