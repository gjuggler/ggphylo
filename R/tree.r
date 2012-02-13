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
  if (is.null(phylo$.tags)) {
    return(list())
  }
  tags <- phylo$.tags[[node]]
  if (is.null(tags)) {
    return(list())
  } else {
    return(tags)
  }
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
tree.foreach <- function(phylo, fn, leaves=T, nodes=T) {
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

#' Turns a phylo object into a data frame
#' @method as.data.frame phylo
#' @seealso \code{\link{tree.as.data.frame}} which this function wraps
as.data.frame.phylo <- function(x, ...) {
  tree.as.data.frame(x, ...)
}

#' Create a data frame from a phylo object where each row corresponds
#' to a node in the tree.
#'
#' @param tree input phylo object
#' @param order.visually boolean, indicating whether to sort the
#'   returned data frame according to how the nodes of the tree would be
#'   arranged along the y-axis in a 2D phylogram plot. (e.g., ((a,b)1,c)2
#'   would produce an order: [a,1,b,2,c]
#' @return data frame with at least the following columns:
#'   \item{label}{string, the node or tip label, or NA}
#'   \item{depth}{integer, the number of nodes between this node and the
#'     furthest leaf. Leaves have a depth of 1, and the root node has the highest depth.}
#'   \item{is_leaf}{boolean, indicates whether this row is a leaf (TRUE) or an internal node (FALSE)}
#'   \item{parent}{integer, the node index of the parent node (or NA for the root node)}
#'   \item{id}{integer, the node index of the current node}
#' \item{[tags]}{If the phylo object has any attached tags (see
#'   \code{\link{get.tags}} and \code{\link{get.tag}}), each unique tag will be incorporated as an additional
#'   column. The \code{\link{rbind.fill}} function from \code{\link{plyr}} is used to combine the tags from different
#'   nodes into one overall data frame.}
#' @importFrom plyr rbind.fill
tree.as.data.frame <- function(tree, order.visually=T) {
  tree.df <- data.frame(stringsAsFactors=F)
  tree.foreach(tree, function(phylo, node) {
    cur.tags <- tree.get.tags(phylo, node)
    cur.tags[['label']] <- tree.label.for.node(phylo, node)
    cur.tags[['depth']] <- tree.leaves.beneath(phylo, node)
    cur.tags[['branch.length']] <- tree.branch.length(phylo, node)
    cur.tags[['parent']] <- tree.parent.node(phylo, node)
    cur.tags[['is.leaf']] <- tree.is.leaf(phylo, node)
    cur.tags[['id']] <- node
    tree.df <<- rbind.fill(tree.df, as.data.frame(cur.tags, stringsAsFactors=F))
  })

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
#' x1 <- tree.normalize.branchlengths(x, push.to.tips=T)
#' x2 <- tree.normalize.branchlengths(x, push.to.tips=F)
#' tree.plot(list('push.to.tips=T'=x1, 'push.to.tips=F'=x2), do.plot=T)
tree.normalize.branchlengths <- function(phylo, push.to.tips=F) {
  n.leaves <- length(phylo$tip.label)
  n.nodes <- length(phylo$tip.label)+phylo$Nnode

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

    is.leaf <- i <= n.leaves

    if (is.leaf) {
      cur.count <- 1
    } else {
      cur.count <- tree.leaves.beneath(phylo, i)
    }

    if (parent.node > -1) {
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
#' can be given either as a node ID or a tip label.
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

#' Returns the index of the node with a given label.
#'
#' @param phylo input phylo object
#' @param label character, the label to search for in the tree.
#' @return integer vector corresponding to the indices of all nodes with the given
#' label. Returns a zero-length vector if no nodes matched.
#' @export
tree.node.with.label <- function(phylo, label) {
  all.labels <- c(phylo$tip.label, phylo$node.label)
  return(which(all.labels %in% label))
}

#' Returns the index of the node with a given label. Alias for \code{\link{tree.node.with.label}}.
#' @param phylo input phylo object
#' @param label character, the label to search for in the tree.
#' @return integer vector corresponding to the indices of all nodes with the given
#' label. Returns a zero-length vector if no nodes matched.
#' @export
tree.find <- function(...) {
  tree.node.with.label(...)
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
tree.label.for.node <- function(tree, node) {
  if (node <= length(tree$tip.label)) {
    return(tree$tip.label[node])
  } else if (node <= (tree$Nnode + length(tree$tip.label))) {
    node.label.index <- node - length(tree$tip.label)
    return(tree$node.label[node.label.index])
  }
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
  return(phylo)
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
  sum(phylo$edge.length, na.rm=T)
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
  drop.tip(phylo, leaf)
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
  not.in.set <- setdiff(phylo$tip.label, leaves)
  tree <- drop.tip(tree, not.in.set)
  tree
}

#' Determine if the given node is a leaf or an internal node. Alias of \code{\link{tree.is.leaf}}.
#' @method is.leaf phylo
#' @param phylo, input phylo object
#' @param node, integer index of the desired node to test
#' @return boolean, TRUE if the given node is a leaf, FALSE if it is an internal node
is.leaf.phylo <- function(phylo, node) {
  tree.is.leaf(phylo, node)
}

#' Determine if the given node is a leaf or an internal node.
#' 
#' @method is.leaf phylo
#' @param phylo, input phylo object
#' @param node, integer index of the desired node to test
#' @return boolean, TRUE if the given node is a leaf, FALSE if it is an internal node
tree.is.leaf <- function(phylo, node) {
  return(node <= length(phylo$tip.label))
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

#' Extracts the parent node ID for the given node. Returns -1 if the node is root.
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
  df.list <- phylo.layout.df(phylo,
    layout.ancestors=T
  )
  
  nodes <- subset(df.list, type=='node')
  nodes[order(nodes$x), 'node']
}

#' Sorts a data frame according to a tree. Requires the data frame to
#' have a 'id' column, which corresponds to the index of each node in the tree. Not much error-checking is done here,
#' so be careful! Uses \code{\link{order.nodes.visually}} to sort nodes.
#'
#' @param tree.df, a data frame containing one row for each node in the tree. This data frame MUST contain a column
#' labeled 'id', which contains, for each row, the index of the corresponding node in the tree.
#' @param phylo, input phylo object
#' @return data.frame, the data frame sorted according to the natural visual order of the tree .
#' @export

sort.df.by.tree <- function(tree.df, phylo) {
  phylo.order <- order.nodes.visually(phylo)

  df.order <- match(phylo.order, tree.df$id)
  tree.df[df.order, ]
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