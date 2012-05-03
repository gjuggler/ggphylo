#' Returns a vector of node indices for all nodes (leaves plus internal nodes) in the tree.
#'
#' @param phylo input phylo object
#' @return integer vector with the node indices for all nodes
#' @export
nodes <- function(phylo) {
  n.nodes <- length(phylo$tip.label) + phylo$Nnode
  1:n.nodes
}

#' Returns a vector of node indices for all the leaves in the tree.
#'
#' @param phylo input phylo object
#' @return integer vector, with the node indices for all leaves
#' @export
leaves <- function(phylo) {
  n.leaves <- length(phylo$tip.label)
  1:n.leaves
}

#' Returns labels corresponding to all or some subset of the nodes in a tree.
#' 
#' @method labels phylo
#' @S3method labels phylo
#' @param phylo input phylo object
#' @param nodes [optional] vector of node indices for which to return labels
#' @export
labels.phylo <- function(phylo, nodes=NA) {
  if (is.null(phylo$node.label)) {
    node.lbls <- rep(NA, phylo$Nnode)
  } else {
    node.lbls <- phylo$node.label
  }
  x <- c(phylo$tip.label, node.lbls)

  if (!is.na(nodes)) {
    x <- x[nodes]
  }
  x  
}

#' Returns the label for the given node. Alias for \code{\link{tree.label.for.node}}. Given
#' a nice short nickname due to its usefulness.
#' 
#' @method label phylo
#' @S3method label phylo
#' @param phylo, input phylo object
#' @param node, the node index for the desired label
#' @export
label.phylo <- function(phylo, node) {
  tree.label.for.node(phylo, node)
}

#' @export
label <- function(object, ...) {
  UseMethod("label")
}

#' Returns the index of the node with a given label.
#'
#' @param phylo input phylo object
#' @param label character, the label to search for in the tree.
#' @param return.one boolean, TRUE if a maximum of one matching node should be returned, FALSE
#'   if all matching nodes should be returned. Defaults to TRUE.
#' @return integer vector corresponding to the indices of all nodes with the given
#' label. Returns a zero-length vector if no nodes matched.
#' @export
tree.node.with.label <- function(phylo, label, return.one=FALSE) {
  all.labels <- c(phylo$tip.label, phylo$node.label)
  all.matches <- match(label, all.labels)
  if (length(all.matches) > 1 && return.one) {
    all.matches[1]
  } else {
    all.matches
  }
}

#' Returns the index of the node with a given label. Alias for \code{\link{tree.node.with.label}}.
#' @param phylo input phylo object
#' @param label character, the label to search for in the tree.
#' @param return.one boolean, TRUE if a maximum of one matching node should be returned, FALSE
#'   if all matching nodes should be returned. Defaults to TRUE.
#' @return integer vector corresponding to the indices of all nodes with the given
#' label. Returns a zero-length vector if no nodes matched.
#' @export
tree.find <- function(...) {
  tree.node.with.label(...)
}

#' Returns the index of the root node in the tree, or NA if there is no root
#'
#' @param phylo input phylo object
#' @return integer vector corresponding to the index of the root node, or NA if there is no root.
#' @export
tree.get.root <- function(phylo) {
  which.match <- which(!(nodes(phylo) %in% phylo$edge[, 2]))
  if (length(which.match) > 0) {
    which.match[1]
  } else {
    NA
  }
}

#' Translates all labels in the tree based on a list of label mappings.
#'
#' @param phylo input phylo object
#' @param map a list of key-value pairs by which to translate node labels
#' @return phylo, the original tree with all matching labels mapped to their new values.
#' @export
#' @examples
#' # Translate a tree of NCBI taxa to their common names
#' taxon.tree <- tree.read('((9606,9598),9593);')
#' taxon.map <- list('9606'='Human', '9598'='Chimpanzee', '9593'='Gorilla')
#' readable.tree <- tree.translate(taxon.tree, taxon.map)
#' print(as.character(readable.tree)) # ((Human,Chimpanzee),Gorilla)
tree.translate <- function(phylo, map) {
  for (node in nodes(phylo)) {
    lbl <- label(phylo, node)
    if (!is.null(lbl) && !is.null(map[[lbl]])) {
      phylo <- tree.set.label(phylo, node, map[[lbl]])
    }
  }
  phylo
}

#' Retrieves a named tag from the given node. Tags are usually loaded
#' via NHX annotations when reading a tree file or string (see
#' \code{\link{tree.read}}) or by loading data from an external data frame
#' (see \code{\link{tree.load.tags}}).
#'
#' @param phylo input phylo object
#' @param node the node index from which to retrieve the tag
#' @param tag the name of the tag to retrieve
#' @return character string corresponding to the tag value, if exists; NA otherwise.
#' @export
#' @examples
#' tree <- tree.read('((a,b[&&NHX:foo=bar]),c);')
#' tree.get.tag(tree, tree.node.with.label(tree, 'b'), 'foo') # gives "bar"
tree.get.tag <- function(phylo, node, tag) {
  tags <- phylo$.tags
  if (is.null(tags)) {
    return(NA)
  }
  cur.tags <- tags[[node]]
  if (is.null(cur.tags)) {
    return(NA)
  }
  tag <- cur.tags[[tag]]
  if(is.null(tag)) {
    return(NA)
  }
  return(tag)
}

#' Retrieves a list of all tags for the given node.
#'
#' @param phylo input phylo object
#' @param node the node index for the desired tags
#' @return list containing all tags associated with this node, if tags exist; empty list otherwise.
#' @export
#' @examples
#' tree <- tree.read('((a,b[&&NHX:foo=bar]),c);')
#' tree.get.tags(tree, tree.node.with.label(tree, 'b')) # foo => bar
tree.get.tags <- function(phylo, node) {
  if (!tree.has.tags(phylo)) {
    return(list())
  }
  
  tags <- phylo$.tags[[node]]
  #print(paste(label(phylo, node), tags))
  if (is.null(tags) || is.na(tags)) {
    return(list())
  } else {
    return(tags)
  }
}

#' Retrieves a list of all tags for the given node. Convenience wrapper around \code{\link{tree.get.tags}}.
#' @method tags phylo
#' @S3method tags phylo
#' @export
tags.phylo <- function(phylo, node) {
  tree.get.tags(phylo, node)
}

#' @export
tags <- function(object, ...) {
  UseMethod("tags")
}

#' Sets the tag for a given node to the specified value. Returns the
#' same tree structure with the new tag in place.
#' 
#' @param phylo input phylo object
#' @param node integer, the index of the node whose tag will be set
#' @param tag string, the tag to set
#' @param value string, the value to set
#' @return phylo, the phylo object with the new tag
#' @export
#' @examples
#' tree <- tree.read('((a,b),c)d;')
#' tree.set.tag(tree, tree.find(tree, 'b'), 'foo', 'bar')
tree.set.tag <- function(phylo, node, tag, value) {
  if (!tree.has.tags(phylo)) {
    phylo$.tags <- as.list(rep(NA, length(nodes(phylo))))
  }

  tags.l <- tags(phylo, node)
  tags.l[[tag]] <- value
  phylo$.tags[[node]] <- tags.l
  phylo
}

#' Determines whether the given phylo object contains tags or not.
#'
#' @param phylo input phylo object
#' @return boolean, indicating this phylo has tags (TRUE) or doesn't (FALSE).
#' @export
#' @examples
#' tree.has.tags(tree.read('((a,b[&&NHX:foo=bar]),c);')) # TRUE
#' tree.has.tags(tree.read('((a,b),c);')) # FALSE
tree.has.tags <- function(phylo) {
  !is.null(phylo$.tags)
}

tree.split.above <- function(phylo, node, fraction) {
  
}


#' Loops over nodes in the tree, applying a function using each node
#' as an argument. Optionally, either the leaves or internal nodes can
#' be dropped from the iteration.
#'
#' @param phylo input phylo object
#' @param fn the function to call. For each node, fn will be called
#' with arguments [phylo, node] where node is the node's index within the phylo object.
#' @export
#' @examples
#' set.seed(1)
#' tree <- rcoal(3)
#' str <- ''
#' tree.foreach(tree, function(x, i) {
#'   str <<- paste(str, ' ', tree.get.label(x, i), '=', tree.depth.to.root(x, i), sep='')
#' })
#' print(str)
tree.foreach <- function(phylo, fn, leaves=TRUE, nodes=TRUE) {
  n.leaves <- length(phylo$tip.label)
  n.internals <- phylo$Nnode

  indices <- c()
  if (leaves) {
    indices <- c(indices, 1:n.leaves)
  }
  if (nodes) {
    indices <- c(indices, (n.leaves+1):(n.leaves+n.internals))
  }
  
  for (node in indices) {
    do.call(fn, as.list(c(phylo, node)))
  }
}

#' Loads a data.frame or CSV file into a tree as tags. Requires a 'label' column
#' in the source data, which is used to merge rows in the dataset to labeled nodes
#' in the tree.
#'
#' @param phylo input phylo object
#' @param x input data.frame or string indicating the location of a CSV file.
#' @return the phylo object, with tags corresponding to the data from the input
#'   dataset.
#' @export
#' @examples
#' tree <- tree.read('((a,b)c,d)e;')
#' x <- data.frame(
#'   label=c('a', 'b', 'c', 'd', 'e'),
#'   value=c(1, 2, 3, 4, 5)
#' )
#' tree <- tree.load.data(tree, x)
#' print(as.character(tree)) # as NHX string
#' print(as.data.frame(tree, minimal.columns=TRUE)) # as data.frame
tree.load.data <- function(phylo, x, ...) {
  if (!is.data.frame(x)) {
    if (is.character(x)) {
      if (file.exists(x)) {
        x <- read.csv(x, ...)
      }
    }
    stop("tree.load.data given a bad input: not a data.frame, and not a valid file name.")
  }

  if (!('label' %in% colnames(x))) {
    stop("tree.load.data requires a 'label' column in the input data to match rows to
      nodes in the tree.")
  }

  if (!tree.has.tags(phylo)) {
    #print("Creating empty tags in phylo object..")
    phylo$.tags <- as.list(rep(NA, length(nodes(phylo))))
  }

  x$label <- as.character(x$label)

  any.non.matches <- FALSE

  tags <- setdiff(colnames(x), c('label'))
  for (i in 1:nrow(x)) {
    lbl <- x[i, 'label']
    match.node <- tree.find(phylo, lbl)
    if (length(match.node) > 0 && !is.na(match.node)) {
      cur.list <- as.list(x[i,])
      cur.list[['label']] <- NULL
      # Don't use 'unlist' -- it turns everything into a character vector!!
      #cur.list <- unlist(cur.list, recursive=FALSE)
      if (!is.list(cur.list)) {
        cur.list <- as.list(cur.list)
      }
      phylo$.tags[[match.node]] <- cur.list
    } else {
      any.non.matches <- TRUE
    }
  }
  if (any.non.matches) {
    warning("tree.load.data found rows in the data that did not match a labeled node in the tree")
  }

  phylo
}
 

#' Turns a phylo object into a data frame
#' @method as.data.frame phylo
#' @S3method as.data.frame phylo
#' @seealso \code{\link{tree.as.data.frame}} which this function wraps
as.data.frame.phylo <- function(x, ...) {
  tree.as.data.frame(x, ...)
}

#' Create a data frame from a phylo object where each row corresponds
#' to a node in the tree.
#'
#' @param tree input phylo object
#' @param minimal.columns boolean, whether to yield a limited number of columns. When TRUE, only the
#'   following columns are produced: 'label', 'node', and columns corresponding to any tags in the
#'   tree. Defaults to FALSE.
#' @param order.visually boolean, indicating whether to sort the
#'   returned data frame according to how the nodes of the tree would be
#'   arranged along the y-axis in a 2D phylogram plot. (e.g., ((a,b)1,c)2
#'   would produce an order: [a,1,b,2,c]
#' @return data frame with at least the following columns:
#'   \enumerate{
#'   \item{label}{string, the node or tip label, or NA}
#'   \item{node}{integer, the node index of the current node}
#' \item{[tags]}{If the phylo object has any attached tags (see
#'   \code{\link{get.tags}} and \code{\link{get.tag}}), each unique tag will be incorporated as an additional
#'   column. The \code{\link{rbind.fill}} function from \code{\link{plyr}} is used to combine the tags from different
#'   nodes into one overall data frame.}
#'   }
#'   If the parameter \code{minimal.columns} is set to FALSE, the following columns will be added:
#'   \enumerate{
#'   \item{depth}{integer, the number of nodes between this node and the
#'     furthest leaf. Leaves have a depth of 1, and the root node has the highest depth.}
#'   \item{is_leaf}{boolean, indicates whether this row is a leaf (TRUE) or an internal node (FALSE)}
#'   \item{parent}{integer, the node index of the parent node (or NA for the root node)}
#'   }
#' @importFrom plyr rbind.fill
tree.as.data.frame <- function(tree, minimal.columns=FALSE, order.visually=TRUE) {
  tree.df <- data.frame(stringsAsFactors=FALSE)
  tree.foreach(tree, function(phylo, node) {
    cur.tags <- tree.get.tags(phylo, node)
    cur.tags[['label']] <- tree.label.for.node(phylo, node)
    cur.tags[['node']] <- node    
    
    if (!minimal.columns) {
      cur.tags[['depth']] <- tree.leaves.beneath(phylo, node)
      cur.tags[['branch.length']] <- tree.branch.length(phylo, node)
      cur.tags[['parent']] <- tree.parent.node(phylo, node)
      cur.tags[['is.leaf']] <- tree.is.leaf(phylo, node)
    }
    
    tree.df <<- rbind.fill(tree.df, as.data.frame(cur.tags, stringsAsFactors=FALSE))
  })

  # Some fixing-up of the column ordering.
  init.cols <- c('node', 'label', 'parent', 'branch.length')
  if (minimal.columns) {
    init.cols <- c('label', 'node')
  }
  remaining <- setdiff(colnames(tree.df), init.cols)
  #print(tree.df)
  #print(init.cols)
  #print(remaining)
  tree.df <- subset(tree.df, select=c(init.cols, remaining))

  # Go through and make columns numeric if applicable.
  cols <- colnames(tree.df)
  for (col.s in cols) {
    vals <- tree.df[, col.s]
    if (!any(is.na(as.numeric(vals)))) {
      tree.df[, col.s] <- as.numeric(tree.df[, col.s])
    } else {
      #print(vals)
    }
  }
  

  if (order.visually) {
    tree.df <- tree.df[order.nodes.visually(tree),]
  }
  tree.df
}

#' Removes internal node labels from a tree
#'
#' @param phylo input phylo object
#' @export
tree.remove.node.labels <- function(phylo) {
  phylo$node.label <- NULL
  phylo
}

#' Removes branch lengths from a tree
#'
#' @param phylo input phylo object
#' @export
tree.strip.branchlengths <- function(phylo) {
  phylo$edge.length <- NULL
  phylo
}

#' Removes branch lengths from a tree
#'
#' @param phylo input phylo object
#' @export
tree.apply.branchlengths <- function(tree.df, phylo, column='branch.length') {
  tree.df$label <- as.character(tree.df$label)
  tree.df <- subset(tree.df, label %in% labels(phylo))

  node.indices <- tree.node.with.label(phylo, tree.df$label)
  edge.indices <- match(node.indices, phylo$edge[, 2])

  col.values <- tree.df[, column]
  col.values <- col.values[!is.na(edge.indices)]
  edge.indices <- edge.indices[!is.na(edge.indices)]
  phylo$edge.length[edge.indices] <- col.values
  phylo
}


#' Normalizes the branch lengths of a tree for cleaner 2-dimensional plotting.
#' 
#' @param phylo input phylo object
#' @param push.to.tips boolean, determines the branch length given to
#' inner clades of the tree. When TRUE, inner clades are given short 
#' external branches and long internal branches (i.e. the clade is 
#' pushed to the tip); when FALSE, inner clades are given short
#' internal branches (i.e. pushed to the root). See the examples for a
#' demonstration of this difference.
#' @return a phylo object with altered branch lengths
#' @export
#' @examples
#' set.seed(1); x <- rcoal(15)
#' x1 <- tree.normalize.branchlengths(x, push.to.tips=TRUE)
#' x2 <- tree.normalize.branchlengths(x, push.to.tips=FALSE)
#' ggphylo(list('push.to.tips=TRUE'=x1, 'push.to.tips=FALSE'=x2), do.plot=TRUE)
tree.normalize.branchlengths <- function(phylo, push.to.tips=FALSE) {
  n.leaves <- length(phylo$tip.label)
  n.nodes <- length(phylo$tip.label)+phylo$Nnode

  if (is.null(phylo$edge.length)) {
    phylo$edge.length <- rep(NA, n.nodes)
  }

  max.depth <- 0
  for (i in 1:n.nodes) {
    depth <- tree.depth.to.root(phylo, i)
    max.depth <- max(depth, max.depth)
  }
  max.depth <- max.depth + 1

  for (i in 1:n.nodes) {
    cur.depth <- tree.depth.to.root(phylo,i)
    parent.node <- tree.parent.node(phylo, i)
    edge.index <- which(phylo$edge[,2]==i)

    is.leaf <- tree.is.leaf(phylo, i)

    if (is.leaf) {
      cur.count <- 1
    } else {
      cur.count <- tree.leaves.beneath(phylo, i)
    }

    if (parent.node > 0) {
      parent.count <- tree.leaves.beneath(phylo, parent.node)
      if (push.to.tips) {
        count.diff <- parent.count - cur.count
        #print(paste(i, count.diff))
        phylo$edge.length[edge.index] <- count.diff
      } else {
        if (is.leaf) {
          # Branch length should equal diff. between depth to root and max depth.
          cur.d <- tree.depth.to.root(phylo, i)
          phylo$edge.length[edge.index] <- max.depth - cur.d
        } else {
          phylo$edge.length[edge.index] <- 1
        }
      }
    } else {
      phylo$edge.length[edge.index] <- 0
    }
  }

  phylo
}

#' Returns the number of leaves beneath a node in the tree (where a
#' leaf node is defined as having 1 leaf beneath)
#'
#' @param phylo input phylo object
#' @param node integer, the index of the node to count leaves beneath
tree.leaves.beneath <- function(phylo, node) {
  if (tree.is.leaf(phylo, node)) {
    return(1)
  }
  cld <- tree.extract.clade(phylo, node)
  length(cld$tip.label)
}


#' Returns the length from the tree root to the given node. Note: the node
#' can be given either as a node index or a tip label.
#'
#' @param phylo input phylo object
#' @param node, either integer or character. When integer, the node
#' index of the desired node to query. When character, the label of the
#' desired node to query.
#' @export
tree.depth.to.root <- function(phylo,node) {
  tip.index <- node
  if (is.character(node)) {
    tip.index <- which(phylo$tip.label==node)
  }
  cur.node.b <- tip.index
  p.edges <- phylo$edge

  length <- 0
  while(length(which(p.edges[,2]==cur.node.b)) > 0) {
    cur.edge.index <- which(p.edges[,2]==cur.node.b)
    cur.edge.length <- 1
    length <- length + cur.edge.length
    cur.node.a <- p.edges[cur.edge.index,1]
    cur.node.b <- cur.node.a # Move up to the next edge
  }
  return(length)
}

#' Replaces the label of a given node and returns a new phylo object
#' 
#' @param phylo input phylo object
#' @param node the index of the node whose label will be replaced
#' @param label the new label for the given node
#' @return a new phylo object, with the node's label replaced
#' @export
tree.set.label <- function(tree, node, label) {
  if (node <= length(tree$tip.label)) {
    tree$tip.label[node] <- label
  } else if (node <= (tree$Nnode + length(tree$tip.label))) {    
    node.label.index <- node - length(tree$tip.label)
    tree$node.label[node.label.index] <- label
  }
  tree
}

#' Returns the label for a given node in the tree. Alias for \code{\link{tree.label.for.node}}.
#'
#' @param phylo input phylo object
#' @param node integer, the node's index
#' @return the label for the given node
#' @export
tree.get.label <- function(tree, node) {
  tree.label.for.node(tree, node)
}

#' Returns the label for a given node in the tree.
#'
#' @param phylo input phylo object
#' @param node integer, the node's index
#' @return the label for the given node
#' @export
tree.label.for.node <- function(phylo, nodes) {
  if (is.null(phylo$node.label)) {
    node.lbls <- rep(NA, phylo$Nnode)
  } else {
    node.lbls <- phylo$node.label
  }
  x <- c(phylo$tip.label, node.lbls)

  if (!is.na(nodes)) {
    x <- x[nodes]
  }
  x
}

#' Extracts the length of the branch above the given node. Returns 0 if the node is root.
#' 
#' @param phylo input phylo object
#' @param node integer, the node's index
#' @return numeric, the branch length of the edge leading to the given node. May be NA.
#' @export
tree.branch.length <- function(phylo, node) {
  edge.index <- which(phylo$edge[,2]==node)
  if (is.null(phylo$edge.length)) {
    return(NA)
  }
  bl <- phylo$edge.length[edge.index]
  if (length(bl)==0) {
    bl <- 0
  }
  return(bl)
}

#' Scales the tree by a given multiplicative factor.
#'
#' @param phylo input phylo object
#' @param factor numeric, the factor to scale all branch lengths by
#' @return phylo, a new phylo object scaled by the desired amount.
#' @export
tree.scale.by <- function(phylo, factor) {
  phylo$edge.length <- phylo$edge.length * factor
  phylo
}

#' Scales the tree to a given total branch length.
#'
#' @param phylo input phylo object
#' @param factor numeric, the desired total branch length of the tree
#' @return phylo, a new phylo object scaled by the required amount.
#' @export
tree.scale.to <- function(phylo, total.length) {
  cur.total <- tree.total.branch.length(phylo)
  scale.factor <- total.length / cur.total
  tree.scale.by(phylo, scale.factor)
}

#' Returns the total branch length of the tree
#'
#' @param phylo input phylo object
#' @return numeric, the sum of all branch lengths in the tree.
#' @export
tree.total.branch.length <- function(phylo) {
  sum(phylo$edge.length, na.rm=TRUE)
}

#' Returns the maximum root-to-tip branch length in the tree. The
#' root-to-tip branch length is calculated for each node, and the
#' maximum value is returned.
#'
#' @param phylo input phylo object
tree.max.length.to.root <- function(phylo) {
  max.length <- max(tree.lengths.to.root(phylo))
  if (is.na(max.length)) {
    max.depth <- 0
    for (i in 1:length(phylo$tip.label)) {
      cur.depth <- tree.depth.to.root(phylo,i)
      max.depth <- max(max.depth, cur.depth)
    }
    return(max.depth)
  }
  return(max.length)
}

#' Returns the mean root-to-tip path length of the tree. This is
#' defined as the mean, across all leaves, of the total branch length
#' from the root to the leaf.
#'
#' @param phylo input phylo object
#' @return numeric, the mean root-to-tip path length.
#' @export
tree.mean.path.length <- function(phylo) {
  mean(tree.lengths.to.root(phylo))
}

#' Returns a vector containing the root-to-tip branch length for each
#' leaf in the tree.
#' 
#' @param phylo input phylo object
#' @return numeric vector, with one entry for each leaf node.
tree.lengths.to.root <- function(phylo) {
  lengths <- c()
  if (length(phylo$tip.label) == 0) {
    return(NA)
  }
  for (i in 1:length(phylo$tip.label)) {
    lengths[i] <- tree.length.to.root(phylo,i)
  }
  lengths
}

#' Returns the length from the tree root to the given node. The input
#'  node can either be input as a node index or a node label.
#' 
#' @param phylo input phylo object
#' @param node integer or character. When integer, the node index; when character, the node label
#' @return numeric, the total branch length separating the tree root and the given node.
#' @export
tree.length.to.root <- function(phylo, node) {
  tip.index <- node
  if (is.character(node)) {
    tip.index <- which(phylo$tip.label==node)
  }
  cur.node.b <- tip.index

  p.edges <- phylo$edge
  p.lengths <- phylo$edge.length

  if(is.null(p.lengths)) {
    p.lengths <- rep(1, length(p.edges[,1]))
  }

  length <- 0
  while(length(which(p.edges[,2]==cur.node.b)) > 0) {
    cur.edge.index <- which(p.edges[,2]==cur.node.b)
    cur.edge.length <- p.lengths[cur.edge.index]
    if (length(cur.edge.length) == 0 || is.na(cur.edge.length)) {
      cur.edge.length <- 0
    }
    length <- length + cur.edge.length
    cur.node.a <- p.edges[cur.edge.index,1]
    cur.node.b <- cur.node.a # Move up to the next edge
  }
  return(length)
}

#' Removes a single leaf from the tree and returns a copy.
#'
#' @param phylo, input phylo object
#' @param leaf, the node index of the desired leaf to remove.
#' @return phylo, a new phylo object with the altered tree structure.
#' @export
tree.remove.leaf <- function(phylo, leaf) {
  tree.remove.leaves(phylo, leaf)
}

tree.remove.leaves <- function(phylo, leaves) {
  cur.df <- as.data.frame(phylo, minimal.columns=TRUE)
  col.s <- setdiff(colnames(cur.df), 'node')
  cur.df <- cur.df[, col.s]
  extra.cols <- setdiff(colnames(cur.df), 'label')
  if (tree.has.tags(phylo)) {
    warning("[tree.remove.leaves] Tags found in tree -- tags will be re-applied after tree manipulation, but this requires unique node labels!")
  }

  phylo$.tags <- NULL
  p2 <- drop.tip(phylo, leaves)
  p2 <- tree.load.data(p2, cur.df)
  p2
}

tree.remove.outgroup <- function(phylo) {
  root.n <- tree.get.root(phylo)
  children <- tree.children(phylo, root.n)
  child.clade.sizes <- c()
  for (i in 1:length(children)) {
    child <- children[i]
    clade.size <- tree.leaves.beneath(phylo, child)
    child.clade.sizes[i] <- clade.size
  }
  
  min.index <- which(order(child.clade.sizes) == 1)
  min.child <- children[min.index]

  tree.remove.clade(phylo, min.child)
}

#' Returns the sub-tree beneath the given node.
#'
#' @param phylo, input phylo object
#' @param node, the node from which to extract the clade / sub-tree
#' @return phylo, a new phylo object corresponding to the sub-tree beneath the given node.
#' @export
tree.extract.clade <- function(phylo, node) {
  extract.clade(phylo, node)
}

#' Returns the minimum sub-tree connecting the given leaf nodes
#'
#' @param phylo, input phylo object
#' @param leaves, integer vector containing the node indices of the leaves from which to extract a sub-tree
#' @return phylo, a phylo object containing the minimum spanning sub-tree defined by the given leaf nodes
#' @export
tree.extract.subtree <- function(phylo, leaves) {
  not.in.set <- setdiff(leaves(phylo), leaves)
  phylo <- tree.remove.leaves(phylo, not.in.set)
  phylo
}

tree.subtree.leaves <- function(phylo, node) {
  if (tree.is.leaf(phylo, node)) {
    cur.leaves <- c(node)
  } else {
    cur.children <- tree.children(phylo, node)
    cur.leaves <- c()
    while (length(cur.children) > 0) {
      cur.child <- cur.children[1]
      cur.children <- cur.children[-1]
      if (tree.is.leaf(phylo, cur.child)) {
        cur.leaves <- c(cur.leaves, cur.child)
      } else {
        cur.children <- c(cur.children, tree.children(phylo, cur.child))
      }
    }
  }
  cur.leaves
}

tree.remove.clade <- function(phylo, node) {
  subclade.leaves <- tree.subtree.leaves(phylo, node)
  remaining.leaves <- setdiff(leaves(phylo), subclade.leaves)
  tree.extract.subtree(phylo, remaining.leaves)
}

#' Determine if the given node is a leaf or an internal node. Alias of \code{\link{tree.is.leaf}}.
#' @method is.leaf phylo
#' @S3method is.leaf phylo
#' @param phylo, input phylo object
#' @param node, integer index of the desired node to test
#' @return boolean, TRUE if the given node is a leaf, FALSE if it is an internal node
is.leaf.phylo <- function(phylo, node) {
  tree.is.leaf(phylo, node)
}

#' Determine if the given node is a leaf or an internal node.
#' 
#' @method is.leaf phylo
#' @S3method is.leaf phylo
#' @param phylo, input phylo object
#' @param node, integer index of the desired node to test
#' @return boolean, TRUE if the given node is a leaf, FALSE if it is an internal node
tree.is.leaf <- function(phylo, node) {
  return(node <= length(phylo$tip.label))
}

tree.children <- function(phylo, node) {
  edge.indices <- which(phylo$edge[,1]==node)
  nodes <- phylo$edge[edge.indices,2]
  if (length(nodes)==0) {
    nodes <- c()
  } else {
    nodes <- nodes
  }
  nodes
}

#' Return a list (not a vector!) of the node indices of the given
#' node's direct children. Returns (-1, -1) if the given node is a leaf.
#'
#' @param phylo, input phylo object
#' @param node, integer index of the node to test
#' @return list, a list containing the integer indices of the nodes directly beneath the given node.
#' @export
tree.child.nodes <- function(phylo, node) {
  edge.indices <- which(phylo$edge[,1]==node)
  nodes <- phylo$edge[edge.indices,2]
  if (length(nodes)==0) {
    nodes <- list(c(-1,-1))
  } else {
    nodes <- list(nodes)
  }
  return(list(nodes))
}

#' Extracts the parent node index for the given node. Returns -1 if the node is root.
#' Return the index of the node directly above the given node. Returns -1 if the given node is root.
#'
#' @param phylo, input phylo object
#' @param node, integer index of the node whose parent is desired
#' @return integer, the index of the parent node or -1 if the given node is root.
#' @export
tree.parent.node <- function(phylo, node) {
  edge.index <- which(phylo$edge[,2]==node)
  node <- phylo$edge[edge.index,1]
  if (length(node)==0) {
    node <- -1
  }
  return(node)
}

#' Returns the indices of all nodes in the tree, arranged in the order
#' in which they would be displayed in a 2-dimensional phylogram plot. In other words,
#' the parent node of two leaf nodes will be placed directly in between its two children,
#' since the parent node is usually placed at the midpoint between two leaf nodes in a plot.
#'
#' @param phylo, input phylo object
#' @return integer vector, the indices of all nodes in the tree, arranged in "visual" order.
#' @export
order.nodes.visually <- function(phylo) {
  phylo <- reorder(phylo, order="cladewise");
  df.list <- tree.layout(phylo,
    layout.ancestors=TRUE
  )
  
  nodes <- subset(df.list, type=='node')
  nodes[order(nodes$y), 'node']
}

#' Sorts a data frame according to a tree. Requires the data frame to
#' have a 'node' column, which corresponds to the index of each node in the tree. Not much error-checking is done here,
#' so be careful! Uses \code{\link{order.nodes.visually}} to sort nodes.
#'
#' @param tree.df, a data frame containing one row for each node in the tree. This data frame MUST contain a column
#' labeled 'node', which contains, for each row, the index of the corresponding node in the tree.
#' @param phylo, input phylo object
#' @return data.frame, the data frame sorted according to the natural visual order of the tree .
#' @export

sort.df.by.tree <- function(tree.df, phylo) {
  phylo.order <- order.nodes.visually(phylo)

  df.order <- match(phylo.order, tree.df$node)
  tree.df[df.order, ]
}

factorize.labels.by.tree <- function(tree.df, phylo, indent.depth=TRUE) {
  phylo.order <- order.nodes.visually(phylo)
  phylo.df <- as.data.frame(phylo)

  phylo.df <- ddply(phylo.df, .(node), function(x) {
    x$label.indent <- paste(x$label, paste(rep(' ', times=(x$depth-1)*2), collapse=''), collapse='')
    x
  })
  node.order <- match(phylo.order, phylo.df$node)
  phylo.df <- phylo.df[node.order,]

  label.order <- match(tree.df$label, phylo.df$label)
  new.lbls <- phylo.df$label.indent[label.order]
  factor(new.lbls, levels=phylo.df$label.indent)
}

#' Returns TRUE if the given object appears to be a valid tree object for the purposes of this package.
#' 
#' @param x, input phylo object
#' @return boolean, TRUE if it is a tree, FALSE if not
#' @export
is.tree <- function(x) {
  if (!is.null(x$edge) && x$Nnode > 0) {
    TRUE
  } else {
    FALSE
  }
}