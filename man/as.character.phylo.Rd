\name{as.character.phylo}
\alias{as.character.phylo}
\title{Turns a phylo object into a Newick format string representation, including NHX annotations
if the tree contains any tags.}
\usage{
  \method{as.character}{phylo} (phylo, digits = 10,
    ignore.tags = F)
}
\arguments{
  \item{phylo}{input phylo object}

  \item{digits}{integer, the number of digits to output
  when formatting branch lengths}

  \item{ignore.tags}{boolean, whether or not to ignore tags
  when creating the Newick / NHX string. When FALSE, tags
  will be included in the string as NHX annotations.}
}
\value{
  character string representing the phylo object
}
\description{
  Turns a phylo object into a Newick format string
  representation, including NHX annotations if the tree
  contains any tags.
}

