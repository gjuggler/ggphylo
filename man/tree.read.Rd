\name{tree.read}
\alias{tree.read}
\title{Reads a tree and returns a phylo object. Note: this function was taken largely from the \code{\link{ape}} code
base, but some alterations were made in order to support NHX formatted trees and to fix some bugs. Only files and strings
with a single tree are allowed.}
\usage{
  tree.read(x, remove.whitespace = F)
}
\arguments{
  \item{x,}{string or filename containing the Newick- or
  NHX- formatted tree string, or a file containing such a
  string, to be read.}

  \item{remove.whitespace}{boolean, TRUE to remove
  whitespace from the tree, FALSE to retain whitespace.
  (defaults to FALSE)}
}
\value{
  phylo object containing the tree (and associated
  annotations, if any)
}
\description{
  Reads a tree and returns a phylo object. Note: this
  function was taken largely from the \code{\link{ape}}
  code base, but some alterations were made in order to
  support NHX formatted trees and to fix some bugs. Only
  files and strings with a single tree are allowed.
}
\examples{
str <- '((a,b),c)'
x <- tree.read(str)
print(as.character(x))
}

