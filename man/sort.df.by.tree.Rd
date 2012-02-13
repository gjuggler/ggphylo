\name{sort.df.by.tree}
\alias{sort.df.by.tree}
\title{Sorts a data frame according to a tree. Requires the data frame to
have a 'id' column, which corresponds to the index of each node in the tree. Not much error-checking is done here,
so be careful! Uses \code{\link{order.nodes.visually}} to sort nodes.}
\usage{
  sort.df.by.tree(tree.df, phylo)
}
\arguments{
  \item{tree.df,}{a data frame containing one row for each
  node in the tree. This data frame MUST contain a column
  labeled 'id', which contains, for each row, the index of
  the corresponding node in the tree.}

  \item{phylo,}{input phylo object}
}
\value{
  data.frame, the data frame sorted according to the
  natural visual order of the tree .
}
\description{
  Sorts a data frame according to a tree. Requires the data
  frame to have a 'id' column, which corresponds to the
  index of each node in the tree. Not much error-checking
  is done here, so be careful! Uses
  \code{\link{order.nodes.visually}} to sort nodes.
}

