#' Reads a tree annotated with NHX-formatted annotations. A phylo object is returned, with 
#' annotations stored as tag-value lists. Several convenience functions are provided for accessing
#' these annotations; see 'See Also' below.
#' 
#' @param str, a string containing the NHX-formatted tree
#' @param coming.from.tree.read, a boolean indicating whether this function is being called from within
#' the \code{\link{tree.read}} function. Since this function may call \code{\link{tree.read}} if no NHX annotations
#' are detected, this is required to avoid an infinite loop.
#' @return phylo object containing the tree and associated annotations.
#' @export
#' @examples
#' str <- '((a,b[&&NHX:foo=bar]),c);'
#' x <- tree.read.nhx(str)
#' print(tree.get.tags(x, tree.find(x, 'b')))
#' print(as.character(x))
tree.read.nhx <- function(str, coming.from.tree.read=F) {
  if (file.exists(str)) {
    # We're reading a file -- load the file into a string.
    str <- readLines(str)
  }

  nhx.matches <- gregexpr("[,\\)\\( ]?([^,\\)\\(]+)?(:?\\d*\\.?\\d*)?\\[.*?\\]", str)
  matches <- nhx.matches[[1]]
  match.pos <- as.numeric(matches)
  match.len <- attr(matches, 'match.length')

  if (match.len[1] == -1 && !coming.from.tree.read) {
    return(tree.read(str))
  }

  match.pos <- match.pos + 1
  match.len <- match.len - 1

  #print(match.pos)
  nhx.strings <- substring(str, match.pos, match.pos+match.len-1)
  #print(nhx.strings)

  labels <- gsub("\\[&&NHX.*", "", nhx.strings)
  labels.no.bl <- gsub(":.*", "", labels)

  #print(labels.no.bl)

  # Go through and splice the stripped-down labels back into the string.
  for (i in 1:length(match.pos)) {
    new.label <- gsub(labels.no.bl[i], paste('zzz',i,'zzz',sep=''), labels[i])
    str <- paste(substr(str, 0, match.pos[i]-1 ), new.label, substr(str, match.pos[i] + match.len[i], nchar(str)), sep='')
    match.pos <- match.pos - match.len[i] + nchar(new.label)
  }
  
  # Parse the Phylo object from the cleaned-up string.
  #print(str)
  tree <- tree.read(str)

  #print(str(tree))

  # Create a list of NHX annotations keyed by the node ID.
  tag.values <- gsub(".*\\[&&NHX:(.*)\\]", "\\1", nhx.strings)
  tagval.list <- strsplit(tag.values, ":")
  names(tagval.list) <- labels.no.bl

  map.list <- llply(tagval.list, function(x) {
    list.out <- list()
    cur.list <- strsplit(x, "=")
    if (length(cur.list) > 0) {
      for (i in 1:length(cur.list)) {
        vec <- cur.list[[i]]
        list.out[vec[1]] = vec[2]
      }
    }
    return(list.out)
  })

  # Replace the labels with the true labels.
  tree$.tags <- list()
  for (i in nodes(tree)) {
    tree$.tags[[i]] <- NA
  }

  for (i in 1:length(match.pos)) {
    cur.node <- tree.node.with.label(tree, paste('zzz', i, 'zzz', sep=''))
    
    leaf <- tree.is.leaf(tree, cur.node)
    real.node.name <- names(map.list)[i]
   if (leaf) {
      tree$tip.label[cur.node] <- real.node.name
    } else {
      tree$node.label[cur.node-length(tree$tip.label)] <- real.node.name
    }
    tree$.tags[[cur.node]] <- map.list[[i]]
  }

  #tree$.tags <- map.list
  return(tree)
}

#' Reads a tree and returns a phylo object. Note: this function was taken largely from the \code{\link{ape}} code
#' base, but some alterations were made in order to support NHX formatted trees and to fix some bugs. Only files and strings
#' with a single tree are allowed.
#' 
#' @param x, string or filename containing the Newick- or NHX- formatted tree string, or a file containing such a string, to be read.
#' @param remove.whitespace boolean, TRUE to remove whitespace from the tree, FALSE to retain whitespace. (defaults to FALSE)
#' @return phylo object containing the tree (and associated annotations, if any)
#' @export
#' @examples
#' str <- '((a,b),c)'
#' x <- tree.read(str)
#' print(as.character(x))
tree.read <- function(x, remove.whitespace=F)
{
    unname <- function(treetext) {
        nc <- nchar(treetext)
	tstart <- 1
	while (substr(treetext, tstart, tstart) != "(" && tstart <= nc)
            tstart <- tstart + 1
	if (tstart > 1)
            return(c(substr(treetext, 1, tstart - 1),
                     substr(treetext, tstart, nc)))
	return(c("", treetext))
    }
    if (is.character(x)) {
      looks.like.tree.string <- grepl("\\(.*\\).*;?", x)
      if (looks.like.tree.string || !file.exists(x)) {
        # Parse as a string.
        tree <- x
      } else {
        tree <- scan(file = x, what = "", sep = "\n", quiet = TRUE)
      }
    } else {
          stop("argument `text' must be of mode character")
    }

    ## Suggestion from Eric Durand and Nicolas Bortolussi (added 2005-08-17):
    if (identical(tree, character(0))) {
        warning("empty character string.")
        return(NULL)
    }

    if (!grepl(';', tree)) {
      # Try adding a semicolon to the end.
      tree <- paste(tree, ';', sep='')
    }
    
    nhx.matches <- gregexpr("[,\\)\\( ]?([^,\\)\\(]+)?(:?\\d*\\.?\\d*)?\\[.*?\\]", tree)
    matches <- nhx.matches[[1]]
    match.pos <- as.numeric(matches)
    match.len <- attr(matches, 'match.length')
    if (match.len[1] != -1) {
      warning("Found NHX annotations within the tree -- consider using tree.read.nhx directly")
      return(tree.read.nhx(tree, coming.from.tree.read=T))
    }

    if (remove.whitespace) {
      tree <- gsub("[ \t]", "", tree)
    }
    tree <- unlist(strsplit(tree, NULL))
    y <- which(tree == ";")
    Ntree <- length(y)
    x <- c(1, y[-Ntree] + 1)
    ## Suggestion from Olivier François (added 2006-07-15):
    if (is.na(y[1])) return(NULL)
    STRING <- character(Ntree)
    for (i in 1:Ntree)
        STRING[i] <- paste(tree[x[i]:y[i]], sep = "", collapse = "")

    tmp <- unlist(lapply(STRING, unname))
    tmpnames <- tmp[c(TRUE, FALSE)]
    STRING <- tmp[c(FALSE, TRUE)]

    obj <- lapply(STRING, clado.build)
    #print(str(obj))
    obj <- obj[[1]]

    # Parse branch lengths out of the node labels.
    tip.lbls <- obj$tip.label
    n.tips <- length(tip.lbls)
    node.lbls <- obj$node.label
    n.nodes <- obj$Nnode
    edge.lengths <- rep(NA, length(obj$edge))
    for (i in 1:(n.nodes+n.tips)) {
      cur.lbl <- tree.get.label(obj, i)
      if (!is.null(cur.lbl)) {
      bl.match <- grep(':.*', cur.lbl)
      if (length(bl.match) > 0) {
        tokens <- unlist(strsplit(cur.lbl, ':'))
        #print(tokens)
        if (i <= n.tips) {
          tip.lbls[i] <- tokens[1]
        } else {
          node.lbls[i - n.tips] <- tokens[1]
        }
        edge.lengths[which(obj$edge[,2] == i)] <- as.numeric(tokens[2])
      }
      }
    }
    obj$tip.label <- tip.lbls
    if (!is.null(node.lbls)) {
      obj$node.label <- node.lbls
    }
    obj$edge.length <- edge.lengths

    obj
}

#' Internal function, taken directly from \code{\link{ape}}.
#' 
#' @param tp, the string to parse.
#' @return a phylo object
clado.build <- function(tp)
{
    add.internal <- function() {
        edge[j, 1] <<- current.node
        node <<- node + 1
        edge[j, 2] <<- current.node <<- node
        index[node] <<- j # set index
        j <<- j + 1
    }
    add.terminal <- function() {
        edge[j, 1] <<- current.node
        edge[j, 2] <<- tip
        index[tip] <<- j # set index
        tip.label[tip] <<- tpc[k]
        k <<- k + 1
        tip <<- tip + 1
        j <<- j + 1
    }
    go.down <- function() {
        l <- index[current.node]
        node.label[current.node - nb.tip] <<- tpc[k]
        k <<- k + 1
        current.node <<- edge[l, 1]
    }
    if (!length(grep(",", tp))) {
        obj <- list(edge = matrix(c(2, 1), 1, 2), Nnode = 1)
        tp <- unlist(strsplit(tp, "[\\(\\);]"))
        obj$tip.label <- tp[2]
        if (tp[3] != "") obj$node.label <- tp[3]
        class(obj) <- "phylo"
        return(obj)
    }
    tsp <- unlist(strsplit(tp, NULL))
    tp <- gsub(")", ")NA", tp)
    tp <- gsub(" ", "", tp)
    tpc <- unlist(strsplit(tp, "[\\(\\),;]"))
    tpc <- tpc[tpc != ""]
    skeleton <- tsp[tsp == "(" | tsp == ")" | tsp == "," | tsp == ";"]
    nsk <- length(skeleton)
    nb.node <- length(skeleton[skeleton == ")"])
    nb.tip <- length(skeleton[skeleton == ","]) + 1
    ## We will assume there is an edge at the root;
    ## if so, it will be removed and put in a vector
    nb.edge <- nb.node + nb.tip
    node.label <- character(nb.node)
    tip.label <- character(nb.tip)

    edge <- matrix(NA, nb.edge, 2)
    current.node <- node <- nb.tip + 1 # node number
    edge[nb.edge, 1] <- 0    # see comment above
    edge[nb.edge, 2] <- node #

    index <- numeric(nb.edge + 1)
    index[node] <- nb.edge
    ## j: index of the line number of edge
    ## k: index of the line number of tpc
    ## tip: tip number
    j <- k <- tip <- 1
    for (i in 2:nsk) {
        if (skeleton[i] == "(") add.internal()      # add an internal branch (on top)
        if (skeleton[i] == ",") {
            if (skeleton[i - 1] != ")") add.terminal()   # add a terminal branch
        }
        if (skeleton[i] == ")") {
            if (skeleton[i - 1] != ")") {   # add a terminal branch and go down one level
                add.terminal()
                go.down()
            }
            if (skeleton[i - 1] == ")") go.down()   # go down one level
        }
    }
    edge <- edge[-nb.edge, ]
    obj <- list(edge = edge, tip.label = tip.label,
                Nnode = nb.node, node.label = node.label)
    obj$node.label <-
        if (all(obj$node.label == "NA")) NULL
        else gsub("^NA", "", obj$node.label)
    class(obj) <- "phylo"
    obj
}
