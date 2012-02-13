\name{tree.as.data.frame}
\alias{tree.as.data.frame}
\title{Create a data frame from a phylo object where each row corresponds
to a node in the tree.}
\usage{
  tree.as.data.frame(tree, order.visually = T)
}
\arguments{
  \item{tree}{input phylo object}

  \item{order.visually}{boolean, indicating whether to sort
  the returned data frame according to how the nodes of the
  tree would be arranged along the y-axis in a 2D phylogram
  plot. (e.g., ((a,b)1,c)2 would produce an order:
  [a,1,b,2,c]}
}
\value{
  data frame with at least the following columns:
  \item{label}{string, the node or tip label, or NA}
  \item{depth}{integer, the number of nodes between this
  node and the furthest leaf. Leaves have a depth of 1, and
  the root node has the highest depth.}
  \item{is_leaf}{boolean, indicates whether this row is a
  leaf (TRUE) or an internal node (FALSE)}
  \item{parent}{integer, the node index of the parent node
  (or NA for the root node)} \item{id}{integer, the node
  index of the current node} \item{[tags]}{If the phylo
  object has any attached tags (see \code{\link{get.tags}}
  and \code{\link{get.tag}}), each unique tag will be
  incorporated as an additional column. The
  \code{\link{rbind.fill}} function from \code{\link{plyr}}
  is used to combine the tags from different nodes into one
  overall data frame.}
}
\description{
  Create a data frame from a phylo object where each row
  corresponds to a node in the tree.
}

