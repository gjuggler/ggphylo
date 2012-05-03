\name{tree.load.data}
\alias{tree.load.data}
\title{Loads a data.frame or CSV file into a tree as tags. Requires a 'label' column
in the source data, which is used to merge rows in the dataset to labeled nodes
in the tree.}
\usage{
  tree.load.data(phylo, x, ...)
}
\arguments{
  \item{phylo}{input phylo object}

  \item{x}{input data.frame or string indicating the
  location of a CSV file.}
}
\value{
  the phylo object, with tags corresponding to the data
  from the input dataset.
}
\description{
  Loads a data.frame or CSV file into a tree as tags.
  Requires a 'label' column in the source data, which is
  used to merge rows in the dataset to labeled nodes in the
  tree.
}
\examples{
tree <- tree.read('((a,b)c,d)e;')
x <- data.frame(
  label=c('a', 'b', 'c', 'd', 'e'),
  value=c(1, 2, 3, 4, 5)
)
tree <- tree.load.data(tree, x)
print(as.character(tree)) # as NHX string
print(as.data.frame(tree, minimal.columns=TRUE)) # as data.frame
}

