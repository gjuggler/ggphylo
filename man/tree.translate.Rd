\name{tree.translate}
\alias{tree.translate}
\title{Translates all labels in the tree based on a list of label mappings.}
\usage{
  tree.translate(phylo, map)
}
\arguments{
  \item{phylo}{input phylo object}

  \item{map}{a list of key-value pairs by which to
  translate node labels}
}
\value{
  phylo, the original tree with all matching labels mapped
  to their new values.
}
\description{
  Translates all labels in the tree based on a list of
  label mappings.
}
\examples{
# Translate a tree of NCBI taxa to their common names
taxon.tree <- read.tree('((9606,9598),9593);')
taxon.map <- list('9606'='Human', '9598'='Chimpanzee', '9593'='Gorilla')
readable.tree <- tree.translate(taxon.tree, taxon.map)
print(as.character(readable.tree)) # ((Human,Chimpanzee),Gorilla)
}

