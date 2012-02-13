## write.tree.R (2010-12-07)
##   Write Tree File in Parenthetic Format
## Copyright 2002-2010 Emmanuel Paradis, Daniel Lawson, and Klaus Schliep
## This file is/was part of the R-package `ape'.
## [Modified from the ape package version by Greg Jordan]

#' Checks a tree for a valid label, optionally removing bad characters on the way. This is currently a no-op,
#' to avoid interfering with the NHX annotation support.
#'
#' @param x, string input
#' @return string
checkLabel <- function(x, ...)
{
    x
}

#' Converts a phylo object into a Newick-format string representation. Taken mostly from the codebase of \code{\link{ape}}.
#'
#' @param phy, input phylo object
#' @param digits, the number of digits to use when formatting branch lengths
#' @return string, the phylo object in a string representation
tree.to.string <- function(phy, digits = 10)
{
    brl <- !is.null(phy$edge.length)
    nodelab <- !is.null(phy$node.label)
    phy$tip.label <- checkLabel(phy$tip.label)
    if (nodelab) phy$node.label <- checkLabel(phy$node.label)
    f.d <- paste("%.", digits, "g", sep = "")
    cp <- function(x){
        STRING[k] <<- x
        k <<- k + 1
    }
    add.internal <- function(i) {
        cp("(")
        desc <- kids[[i]]
        for (j in desc) {
            if (j > n) add.internal(j)
            else add.terminal(ind[j])
            if (j != desc[length(desc)]) cp(",")
        }
        cp(")")
        if (nodelab && i > n) cp(phy$node.label[i - n]) # fixed by Naim Matasci (2010-12-07)
        if (brl && !is.na(phy$edge.length[ind[i]])) {
            cp(":")
            cp(sprintf(f.d, phy$edge.length[ind[i]]))
        }
    }
    add.terminal <- function(i) {
        cp(phy$tip.label[phy$edge[i, 2]])
        if (brl && !is.na(phy$edge.length[i])) {
            cp(":")
            cp(sprintf(f.d, phy$edge.length[i]))
        }
    }

    n <- length(phy$tip.label)

    ## borrowed from phangorn:
    parent <- phy$edge[, 1]
    children <- phy$edge[, 2]
    kids <- vector("list", n + phy$Nnode)
    for (i in 1:length(parent))
        kids[[parent[i]]] <- c(kids[[parent[i]]], children[i])

    ind <- match(1:max(phy$edge), phy$edge[, 2])

    LS <- 4*n + 5
    if (brl) LS <- LS + 4*n
    if (nodelab)  LS <- LS + n
    STRING <- character(LS)
    k <- 1
    cp("(")
    getRoot <- function(phy)
        phy$edge[, 1][!match(phy$edge[, 1], phy$edge[, 2], 0)][1]
    root <- getRoot(phy) # replaced n+1 with root - root has not be n+1
    desc <- kids[[root]]
    for (j in desc) {
        if (j > n) add.internal(j)
        else add.terminal(ind[j])
        if (j != desc[length(desc)]) cp(",")
    }

    if (is.null(phy$root.edge)) {
        cp(")")
        if (nodelab) cp(phy$node.label[1])
        cp(";")
    }
    else {
        cp(")")
        if (nodelab) cp(phy$node.label[1])
        cp(":")
        cp(sprintf(f.d, phy$root.edge))
        cp(";")
    }
    paste(STRING, collapse = "")
}

#' Turns a phylo object into a Newick format string representation, including NHX annotations
#' if the tree contains any tags.
#' 
#' @method as.character phylo
#' @param phylo input phylo object
#' @param digits integer, the number of digits to output when formatting branch lengths
#' @param ignore.tags boolean, whether or not to ignore tags when creating the Newick / NHX string. When FALSE,
#'    tags will be included in the string as NHX annotations.
#' @return character string representing the phylo object
#' @export
as.character.phylo <- function(phylo, digits=10, ignore.tags=F) {
  if (!ignore.tags && tree.has.tags(phylo)) {
    phylo <- tree.tags.to.labels(phylo)
  }
  tree.to.string(phylo, digits=digits)
}

#' Converts all tags from a tree into NHX annotations and applies those annotations to each node's label.
#' @param phylo input phylo object
#' @return a phylo object, where tags have been concatenated to each node's label
tree.tags.to.labels <- function(phylo) {
  tree.foreach(phylo, function(x, node) {
    tags <- tree.get.tags(x, node)
    if (length(tags) > 0) {
      nhx.str <- tags.to.nhx.string(tags)
      cur.lbl <- tree.label.for.node(x, node)
      phylo <<- tree.set.label(phylo, node, paste(cur.lbl, nhx.str, sep=''))
    }
  })
  phylo
}

#' Writes the phylo object to disk in NHX format.
#' 
#' @method write.nhx
#' @param tree input phylo object
#' @param file the file to write to
write.nhx <- function(phylo, f) {
  str <- as.character(phylo)
  cat(str, file=f)
}

#' Converts a list of tags to an NHX annotation string.
#'
#' @param tags list of tags
#' @return character containing [&&NHX:key=val] style annotations for all tags
tags.to.nhx.string <- function(tags) {
  key.vals <- paste(names(tags), tags, sep='=')
  if (grepl("[\\[\\]\\:\\;\\,\\=]", paste(key.vals, sep=''))) {
    warning("Bad characters within NHX annotations == probably not a good thing!")
  }
  paste('[&&NHX:', key.vals, ']', sep='', collapse='')
}
