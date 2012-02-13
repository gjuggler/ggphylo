\name{tree.foreach}
\alias{tree.foreach}
\title{Loops over nodes in the tree, applying a function using each node
as an argument. Optionally, either the leaves or internal nodes can
be dropped from the iteration.}
\usage{
  tree.foreach(phylo, fn, leaves = T, nodes = T)
}
\arguments{
  \item{phylo}{input phylo object}

  \item{fn}{the function to call. For each node, fn will be
  called with arguments [phylo, node] where node is the
  node's index within the phylo object.}
}
\description{
  Loops over nodes in the tree, applying a function using
  each node as an argument. Optionally, either the leaves
  or internal nodes can be dropped from the iteration.
}
\examples{
set.seed(1)
tree <- rcoal(3)
str <- ''
tree.foreach(tree, function(x, i) {
  str <<- paste(str, ' ', tree.get.label(x, i), '=', tree.depth.to.root(x, i), sep='')
})
print(str)
}

