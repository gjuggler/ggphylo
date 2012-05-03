\name{tree.set.tag}
\alias{tree.set.tag}
\title{Sets the tag for a given node to the specified value. Returns the
same tree structure with the new tag in place.}
\usage{
  tree.set.tag(phylo, node, tag, value)
}
\arguments{
  \item{phylo}{input phylo object}

  \item{node}{integer, the index of the node whose tag will
  be set}

  \item{tag}{string, the tag to set}

  \item{value}{string, the value to set}
}
\value{
  phylo, the phylo object with the new tag
}
\description{
  Sets the tag for a given node to the specified value.
  Returns the same tree structure with the new tag in
  place.
}
\examples{
tree <- tree.read('((a,b),c)d;')
tree.set.tag(tree, tree.find(tree, 'b'), 'foo', 'bar')
}

