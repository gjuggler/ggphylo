\name{tree.normalize.branchlengths}
\alias{tree.normalize.branchlengths}
\title{Normalizes the branch lengths of a tree for cleaner 2-dimensional plotting.}
\usage{
  tree.normalize.branchlengths(phylo, push.to.tips = F)
}
\arguments{
  \item{phylo}{input phylo object}

  \item{push.to.tips}{boolean, determines the branch length
  given to inner clades of the tree. When TRUE, inner
  clades are given short external branches and long
  internal branches (i.e. the clade is pushed to the tip);
  when FALSE, inner clades are given short internal
  branches (i.e. pushed to the root). See the examples for
  a demonstration of this difference.}
}
\value{
  a phylo object with altered branch lengths
}
\description{
  Normalizes the branch lengths of a tree for cleaner
  2-dimensional plotting.
}
\examples{
set.seed(1); x <- rcoal(15)
x1 <- tree.normalize.branchlengths(x, push.to.tips=T)
x2 <- tree.normalize.branchlengths(x, push.to.tips=F)
tree.plot(list('push.to.tips=T'=x1, 'push.to.tips=F'=x2), do.plot=T)
}

