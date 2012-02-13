\name{tree.read.nhx}
\alias{tree.read.nhx}
\title{Reads a tree annotated with NHX-formatted annotations. A phylo object is returned, with
annotations stored as tag-value lists. Several convenience functions are provided for accessing
these annotations; see 'See Also' below.}
\usage{
  tree.read.nhx(str, coming.from.tree.read = F)
}
\arguments{
  \item{str,}{a string containing the NHX-formatted tree}

  \item{coming.from.tree.read,}{a boolean indicating
  whether this function is being called from within the
  \code{\link{tree.read}} function. Since this function may
  call \code{\link{tree.read}} if no NHX annotations are
  detected, this is required to avoid an infinite loop.}
}
\value{
  phylo object containing the tree and associated
  annotations.
}
\description{
  Reads a tree annotated with NHX-formatted annotations. A
  phylo object is returned, with annotations stored as
  tag-value lists. Several convenience functions are
  provided for accessing these annotations; see 'See Also'
  below.
}
\examples{
str <- '((a,b[&&NHX:foo=bar]),c);'
x <- tree.read.nhx(str)
print(tree.get.tags(x, tree.find(x, 'b')))
print(as.character(x))
}

