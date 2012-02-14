#' Plots a tree or a list of trees.
#'
#' @export
tree.plot <- function(
  x,
  # Options common to both the tree and aln plots.
  plot.ancestors = T,
  axis.text.size = 5,
  keep.aspect.ratio = F,

  radial.plot = F,

  # Options just for the tree plot.
  do.plot = T,
  show.y.gridlines = F,
  plot.legend = T,
  xlim = NA,
  xlim.expand = c(0.05, 0.3),
  xlab = "Branch Length",
  labels.in.axis = F,
  labels.in.tree = T,

  line.plot = TRUE,
  line.size = 1.5,
  line.size.by = NA,
  line.color.by = NA,

  node.plot = FALSE,
  node.size = 1,
  node.size.by = NA,
  node.color.by = NA,

  label.plot = TRUE,
  label.size = 4,
  label.size.by = NA,
  label.color.by = NA,

  y.lim = NA,
  ...
){

  if (is.list(x) && !is.tree(x)) {
    layout.df <- ldply(x, function(cur.tree) {
      out.df <- tree.layout(cur.tree,
        layout.ancestors = plot.ancestors
      )
      out.df
    })
    trees.df <- ldply(x, function(cur.tree) {
      out.df <- data.frame(
        max.length = tree.max.length.to.root(cur.tree),
        n.leaves = length(cur.tree$tip.label)
      )
      out.df
    })
    trees.df

  } else {
    if ('PhyloSim' %in% class(x)) {
      phylo <- x$.phylo
      aln <- x$.alignment
    } else {
      phylo <- x
      aln <- NA
    }

    layout.df <- tree.layout(phylo,
      layout.ancestors=plot.ancestors
    )
    trees.df <- data.frame(
      max.length = tree.max.length.to.root(phylo),
      n.leaves = length(phylo$tip.label)
    )
  }

  tree.df <- subset(layout.df, type=='line')
  nodes.df <- subset(layout.df, type=='node')

  #print(head(nodes.df))

#  names <- phylo$tip.label
#  if (is.na(y.lim)) {
#    y.lim <- c(0,length(names)+1)
#    y.breaks <- 1:length(names)
#    y.labels <- names
#  }

  max.length <- max(trees.df$max.length)
  max.n <- max(trees.df$n.leaves)

  q <- ggplot(tree.df)
    
  q <- q + opts(
    axis.text.y = theme_text(size=axis.text.size,hjust=1),
    axis.title.y = theme_blank()
  )
    
  # This creates the layer with the tree branches.
  if (line.plot) {
    aes.list <- list(
      x='x',
      y='y',
      xend='xend',
      yend='yend'
    )
    if (!is.na(line.color.by)) {
      aes.list[['colour']] <- line.color.by
    }
    if (!is.na(line.size.by)) {
      aes.list[['size']] <- line.size.by
    }
    q <- q + geom_joinedsegment(size=line.size, do.call(aes_string, aes.list))
  }

  if (node.plot) {
    aes.list <- list(
      x='x',
      y='y'
    )
    if (!is.na(node.size.by)) {
      aes.list[['size']] <- node.size.by
    }
    if (!is.na(node.color.by)) {
      aes.list[['colour']] <- node.color.by
    }      
    q <- q + geom_point(data=nodes.df, size=2, do.call(aes_string, aes.list))
  }

  # Put labels directly into the tree
  if (label.plot) {
    leaf.nodes <- subset(nodes.df, is.leaf==TRUE)
    aes.list <- list(
      x='x',
      y='y',
      label='label'
    )
    if (!is.na(label.size.by)) {
      aes.list[['size']] <- label.size.by
    }
    if (!is.na(label.color.by)) {
      aes.list[['colour']] <- label.color.by
    }      
    if (radial.plot) {
      aes.list[['angle']] <- '-360 * (x)/(diff(range(x))+1) + 90;'
    }

    text.nodes <- leaf.nodes
    text.margin <- 1 / 100
    text.nodes$y <- text.nodes$y + max.length * text.margin
    q <- q + geom_text(data=text.nodes, size=label.size, hjust=0, do.call(aes_string, aes.list))
  }

  # Tree x-limits.
  if (!is.na(xlim)) {
    tree.limits <- xlim
  } else {
    tree.limits <- c(0, max.length)
  }

  if (any(!is.na(xlim.expand))) {
    range <- tree.limits[2] - tree.limits[1]
    tree.limits[1] <- tree.limits[1] - range*xlim.expand[1]
    tree.limits[2] <- tree.limits[2] + range*xlim.expand[2]
  }

  if (radial.plot) {
    q <- q + scale_x_continuous(limits=c(0, max.n))
    q <- q + coord_polar(expand=F)
    q <- q + opts(axis.text.x=theme_blank())
  } else {
    q <- q + coord_flip()
    q <- q + scale_y_continuous(name=xlab, limits=tree.limits)
    q <- q + scale_x_continuous()
  }

#  q <- q + scale_x_continuous(limits=tree.limits, expand=c(0,0), name=xlab)

  if (nrow(trees.df) > 1) {
    q <- q + facet_wrap( ~ .id)
  }

  if (!plot.legend) {
    q <- q + opts(legend.position='none')
  }

  if (!show.y.gridlines) {
    q <- q + opts(panel.grid.major = theme_blank())
  }

  q <- q + theme_phylo_default()

  if (do.plot) {
    print(q)
  }
  return(invisible(q))
}

#' Returns a data frame defining segments to draw the phylogenetic tree.
#' 
#' @export
tree.layout <- function(
  phylo,
  layout.ancestors = FALSE,
  align.seq.names = NA
) {
  # Number of nodes and leaves.
  n.nodes <- length(phylo$tip.label)+phylo$Nnode
  n.leaves <- length(phylo$tip.label)

  t.labels <- phylo$tip.label
  n.labels <- ((n.leaves+1):n.nodes)
  if (!is.null(phylo$node.label)) {
    n.labels <- phylo$node.label
  }

  # Create the skeleton data frame.
  df <- data.frame(
                   node=c(1:n.nodes),                                            # Nodes with IDs 1 to N.
                   x=0,                                                          # These will contain the x and y coordinates after the layout procedure below.
                   y=0,
                   label=c(t.labels, n.labels),            # The first n.leaves nodes are the labeled tips.
                   is.leaf=c(rep(TRUE, n.leaves), rep(FALSE, n.nodes-n.leaves)),    # Just for convenience, store a boolean whether it's a leaf or not.
                   parent=0,                                                     # Will contain the ID of the current node's parent
                   children=0,                                                   # Will contain a list of IDs of the current node's children
                   branch.length=0                                               # Will contain the branch lengths
                   )

  # Collect the parents, children, and branch lengths for each node
  parent <- c()
  bl <- c()
  children <- list()
  event.count <- c()
  for (i in 1:nrow(df)) {
    node <- df[i,]$node
    parent <- c(parent, tree.parent.node(phylo, node))
    bl <- c(bl, tree.branch.length(phylo, node))
    children <- c(children, tree.child.nodes(phylo, node))
  }
  df$parent <- parent
  df$branch.length <- bl
  df$children <- children

  # Start the layout procedure by equally spacing the leaves in the y-dimension.
  df[df$is.leaf==TRUE,]$y = c(1:n.leaves)

  found.any.internal.node.sequences <- FALSE

  # For each leaf: travel up towards the root, laying out each internal node along the way.
  for (i in 1:n.leaves) {
    cur.node <- i
    while (length(cur.node) > 0 && cur.node != -1) {
      # We always use branch lengths: x-position is simply the length to the root.
      df[cur.node,]$x <- tree.length.to.root(phylo,cur.node)

      # The y-position for internal nodes is the mean of the y-position of all children.
      children <- unlist(df[cur.node,]$children)
      if (length(children) > 0 && children[1] != -1) {
        child.y.sum <- 0
        for (i in 1:length(children)) {
          child.index <- children[i]
          cur.y <- df[child.index,]$y
          child.y.sum <- child.y.sum + cur.y
        }
        df[cur.node, ]$y <- (child.y.sum) / length(children)
      }

      # Try to find the index of this node in the alignment names.
      if (!is.na(align.seq.names)) {
        lbl <- df[cur.node,]$label
        index.in.names <- which(align.seq.names == lbl | align.seq.names %in% c(paste('Node',lbl),paste('Root node',lbl)))
        if (length(index.in.names)>0) {
          df[cur.node,]$y <- index.in.names
          if (!df[cur.node,]$is.leaf) {
            found.any.internal.node.sequences <- TRUE
          }
        }
      }

      cur.node <- unlist(df[cur.node,]$parent)
    }
  }

  df$dir <- 'node'

  # We have a data frame with each node positioned.
  # Now we go through and make two line segments for each node (for a 'square corner' type tree plot).
  line.df <- data.frame()
  for (i in 1:nrow(df)) {
    row <- df[i,]            # Data frame row for the current node.
    if (row$parent == -1) {
      next; # Root node!
    }
    p.row <- df[row$parent,] # Data frame row for the parent node.
    if (layout.ancestors && found.any.internal.node.sequences) {
      horiz.line <- data.frame(
                               node=row$node,
                               y=row$x,
                               yend=p.row$x,
                               x=row$y,
                               xend=p.row$y,
                               label=row$label,
                               dir='horiz',
                               branch.length=row$branch.length
                               )
      line.df <- rbind(line.df,horiz.line)
    } else {
      horiz.line <- data.frame(
                               node=row$node,
                               y=row$x,
                               yend=p.row$x,
                               x=row$y,
                               xend=row$y,
                               label=row$label,   
                               dir='horiz',
                               branch.length=row$branch.length
                               )    # First a line from row.x to parent.
      line.df <- rbind(line.df, horiz.line)

      vert.line <- data.frame(
                               node=row$node,
                               y=p.row$x,
                               yend=p.row$x,
                               x=row$y,
                               xend=p.row$y,
                               label=row$label,
                               dir='vert',
                               branch.length=row$branch.length
          
      )
      line.df <- rbind(line.df, vert.line)
    }
  }

  line.df <- tags.into.df(phylo, line.df)
  df <- tags.into.df(phylo, df)

  # Remove weird list-of-lists from the nodes data frame.
  df$children <- NULL

  line.df <- line.df[rev(1:nrow(line.df)),]
  df <- df[rev(1:nrow(df)),]

  old.x <- df$x
  df$x <- df$y  
  df$y <- old.x

  line.df$type <- 'line'
  df$type <- 'node'
  all.df <- rbind.fill(line.df, df)
  all.df
}

tags.into.df <- function(phylo, df) {
  all.tags <- c()

  for (i in 1:nrow(df)) {
    row <- df[i,]
    node <- row$node
    tags <- tree.get.tags(phylo, node)
    for (j in 1:length(tags)) {
      all.tags <- c(all.tags, names(tags)[j])
    }
  }
  all.tags <- unique(all.tags)

  df[, all.tags] <- 0

  for (i in 1:nrow(df)) {
    row <- df[i,]
    node <- row$node
    tags <- tree.get.tags(phylo, node)

    # Always leave '0' value for tags in the vertical lines.
    # This stops vertical lines from becoming red if the leftmost node
    # has a tag value of 1
    if (row$dir == 'vert') {
      next()
    }

    for (j in 1:length(tags)) {
      key <- names(tags)[j]
      val <- tags[j]
      df[i, key] <- val
    }
  }
  return(df)
}

#' A mostly-blank theme for simple plots.
#' Ruthlessly taken from ggdendro, courtesy of Andrie de Vries
#' http://cran.r-project.org/web/packages/ggdendro/index.html
#' 
#' @export
theme_phylo_blank <- function(){
  stopifnot(require(ggplot2))
  theme_blank <- ggplot2::theme_blank
  ggplot2::opts(
      panel.grid.major = theme_blank(),
      panel.grid.minor = theme_blank(),
      panel.background = theme_blank(),
      axis.title.x = theme_text(colour=NA),
      axis.title.y = theme_blank(),
      axis.text.x = theme_blank(),
      axis.text.y = theme_blank(),
      axis.line = theme_blank(),
      axis.ticks = theme_blank()
  )
}

theme_phylo_black <- function() {
  stopifnot(require(ggplot2))
  theme_blank <- ggplot2::theme_blank
  ggplot2::opts(
      panel.grid.major = theme_segment(colour = "white"),
      panel.grid.minor = theme_segment(colour="white"),
      panel.background = theme_rect(fill = "black", colour = "black"),
      axis.title.x = theme_text(colour='white'),
      axis.title.y = theme_text(colour='white'),
      axis.text.x = theme_text(colour='white'),
      axis.text.y = theme_text(colour='white'),
      axis.line = theme_segment(colour='white'),
      axis.ticks = theme_segment(colour='white')
  )
}

theme_phylo_black <- function() {
  stopifnot(require(ggplot2))
  theme_blank <- ggplot2::theme_blank
  ggplot2::opts(
      panel.grid.major = theme_segment(colour = "white"),
      panel.grid.minor = theme_segment(colour="white"),
      panel.background = theme_rect(fill = "black", colour = "black"),
      axis.title.x = theme_text(colour='white'),
      axis.title.y = theme_text(colour='white'),
      axis.text.x = theme_text(colour='white'),
      axis.text.y = theme_text(colour='white'),
      axis.line = theme_segment(colour='white'),
      axis.ticks = theme_segment(colour='white')
  )
}


theme_phylo_default <- function() {
  stopifnot(require(ggplot2))
  stopifnot(require(grid))
  ggplot2::opts(
    plot.margin = unit(c(0,0,0,0),'npc'),
    axis.title.y = theme_blank(),
    plot.background = theme_blank(),
    panel.border = theme_blank()
  )
}

#' Create a geom_joinedsegment which joins together segments with 'round' corners. Looks MUCH nicer.
#' See https://groups.google.com/forum/?fromgroups#!topic/ggplot2/movv0f_MSuY
GeomJoinedSegment <- proto(ggplot2:::GeomSegment, {
 objname <- "geom_joinedsegment"
 draw <- function(., data, scales, coordinates, arrow=NULL, ...) {
   if (is.linear(coordinates)) {
     return(with(coord_transform(coordinates, data, scales),

       segmentsGrob(x, y, xend, yend, default.units="native",
       gp = gpar(col=alpha(colour, alpha), lwd=size * .pt,
         lty=linetype, lineend = "round"),
       arrow = arrow)
     ))
   }
    data$group <- 1:nrow(data)
    starts <- subset(data, select = c(-xend, -yend))
    ends <- rename(subset(data, select = c(-x, -y)), c("xend" = "x", "yend" = "y"))
    
    pieces <- rbind(starts, ends)
    pieces <- pieces[order(pieces$group),]
    
    GeomPath$draw_groups(pieces, scales, coordinates, arrow = arrow, lineend="round", ...)
}})

#' Create the function for geom_joinedsegment
geom_joinedsegment <- function(mapping = NULL, data = NULL, stat =
"identity", position = "identity", arrow = NULL, ...)  {
  GeomJoinedSegment$new(mapping = mapping, data = data, stat = stat,
        position = position, arrow = arrow, ...)
}
