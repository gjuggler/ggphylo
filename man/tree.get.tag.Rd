\name{tree.get.tag}
\alias{tree.get.tag}
\title{Retrieves a named tag from the given node. Tags are usually loaded
via NHX annotations when reading a tree file or string (see
\code{\link{tree.read}}) or by loading data from an external data frame
(see \code{\link{tree.load.tags}}).}
\usage{
  tree.get.tag(phylo, node, tag)
}
\arguments{
  \item{phylo}{input phylo object}

  \item{node}{the node index from which to retrieve the
  tag}

  \item{tag}{the name of the tag to retrieve}
}
\value{
  character string corresponding to the tag value, if
  exists; NA otherwise.
}
\description{
  Retrieves a named tag from the given node. Tags are
  usually loaded via NHX annotations when reading a tree
  file or string (see \code{\link{tree.read}}) or by
  loading data from an external data frame (see
  \code{\link{tree.load.tags}}).
}
\examples{
tree <- tree.read('((a,b[&&NHX:foo=bar]),c);')
tree.get.tag(tree, tree.node.with.label(tree, 'b'), 'foo') # gives "bar"
}

