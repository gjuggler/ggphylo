#' Plots a tree or a list of trees using \code{\link{ggplot}}.
#'
#' @param x input phylo object, list of phylo objects, or a Newick- or NHX-format tree string.
#'   When a list is given, all trees will be arranged using
#'   \code{\link{facet_wrap}} into a grid of plots. Trees of different
#'   shapes and sizes can be included; the plotting area will be expanded
#'   to fit the largest input tree.
#' @param extra.data an optional data.frame or string pointing to a CSV file, containing data to be
#'   linked to the input phylo object(s). The input data must contain a column named 'label', which
#'   will be used to match each row in the data to a node in the input tree(s). Defaults to NA.
#' @param layout a string indicating how the tree will be laid out. Defaults to 'defualt'.
#'   Available layouts include:
#'   \item{default} The 'default' layout, which arranges the tree in the standard 'phylogram' 
#'     arrangement with branch length along the x-axis and tree nodes evenly spaced along
#'     the y-axis.
#'   \item{unrooted} The 'unrooted' layout, which arranges the
#'     tree from root-to-tip giving equal viewing angle to each clade
#'     according to its size.
#'   \item{radial} The 'radial' layout, which is equivalent to the 'default' layout but plotted
#'     in polar coordinates.
#' @param do.plot a boolean indicating whether the plot should be printed to the current
#'   graphics device. When FALSE, the plot is not printed; the graphics grob is silently returned.
#'   Defaults to TRUE.
#' @param plot.internals boolean, whether the labels of internal nodes should be shown. Defaults
#'   to FALSE.
#' @param x.lab string, the label given to the x-axis (or, when layout='unrooted', given to both
#'   axes). Defaults to 'Branch Length'.
#' @param x.lim vector of length 2, the lower and upper limits to apply to the x-axis. Defaults
#'   to NA, in which case the plot is expanded to include the entire tree.
#' @param y.lim vector of length 2, the lower and upper limits to apply to the y-axis. Defaults
#'   to NA, in which case the plot is expanded to include the entire tree.
#' @param x.expand vector of length 2, the fraction by which to expand the x-axis limits to the
#'   left and to the right. This is useful to allow extra space for trees with long labels or 
#'   large node sizes. Defaults to c(0.05, 0.3) which extends the x-axis limits 5% to the left
#'   and 30% to the right of the default size.
#' @param [line|node|label].[size|color|alpha].by string, indicating the name of a tag value by
#'   which to modulate the given visual property. A value of 0 will be given to nodes which
#'   do not contain the given tag, and unless a 'x.y.scale' parameter is given (see below),
#'   the default \ggplot scaling and size/color/alpha rules will be applied.
#'   The combination of 3 visual elements and 3 visual properties creates 9 possible parameters
#'   to choose from:
#'   \item{line.size.by}
#'   \item{line.color.by}
#'   \item{line.alpha.by}
#'   \item{node.size.by}
#'   \item{node.color.by}
#'   \item{node.alpha.by}
#'   \item{label.size.by}
#'   \item{label.color.by}
#'   \item{label.alpha.by}
#' @param [line|node|label].[size|color|alpha].scale function, used in conjunction with a
#'   corresponding 'x.y.by' parameter (e.g. node.scale.by) to specify the type, limits and range
#'   of the visual scaling. The value is usually the result of one of the \code{\link{ggplot}}
#'   scale convenience functions, such as:
#'     \item{node.size.scale=scale_size_continuous(limits=c(0,100), range=c(1, 5))}
#'   In this example, 'limits' controls the range of tag values to be shown, and 'range'
#'   controls the range of the resulting visual scale (i.e. the node size will range from 1 to 5).
#'   See \code{\link{scale_size_continuous}} or the other ggplot scale_x_y functions for more info.
#' @param [line|node|label].[size|color|alpha] string or numeric, used to specify a constant
#'   value for a visual parameter if no data-based scaling is specified. Default values are:
#'   \item{line.size} 1
#'   \item{line.color} '#888888'
#'   \item{line.alpha} 1
#'   \item{node.size} 1.5
#'   \item{node.color} '#333333'
#'   \item{node.alpha} 1
#'   \item{label.size} 3
#'   \item{label.color} 'black'
#'   \item{label.alpha} 1
#' @return the \code{\link{ggplot}} grob list.
#' @export
ggphylo <- function(
  x,
  extra.data=NA,
  layout='default',
  do.plot = T,
  plot.internals=F,
  x.lab = "Branch Length",
  x.lim = NA,
  x.expand = c(0.05, 0.3),
  y.lim = NA,
  ...
){

  if (is.list(x) && !is.tree(x)) {
    if (is.null(names(x))) {
      names(x) <- paste("Tree #", 1:length(x), sep='')
    }

    layout.df <- ldply(x, function(cur.tree) {
      if (!is.na(extra.data)) {
        cur.tree <- tree.load.data(cur.tree, extra.data)
      }      
      out.df <- tree.layout(cur.tree,
        layout = layout
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
    } else if ('phylo' %in% class(x)) {
      phylo <- x
      aln <- NA
    } else if (is.character(x)) {
      phylo <- tree.read(x)
    }

    if (!is.na(extra.data)) {
      phylo <- tree.load.data(phylo, extra.data)
    }      

    layout.df <- tree.layout(phylo,
      layout = layout
    )
    trees.df <- data.frame(
      max.length = tree.max.length.to.root(phylo),
      n.leaves = length(phylo$tip.label)
    )
  }

  lines.df <- subset(layout.df, type=='line')
  nodes.df <- subset(layout.df, type=='node')
  labels.df <- subset(layout.df, type=='label')

  max.length <- max(trees.df$max.length)
  max.n <- max(trees.df$n.leaves)

  q <- ggplot(lines.df)

  args <- list(...)  
  for (aes.type in c('line', 'node', 'label')) {
    plot.s <- paste('plot', paste(aes.type, 's', sep=''), sep='.')
    
    if (is.null(args[[plot.s]]) || args[[plot.s]] == T) {
      #print(paste(plot.s, args[[plot.s]]))
      aes.list <- switch(aes.type, 
        'line' = list(
          x='x',
          y='y',
          xend='xend',
          yend='yend'
        ),
        'node' = list(
          x='x',
          y='y'
        ),
        'label' = list(
          x='x',
          y='y',
          label='label'
        )
      )

      if (aes.type == 'line') {
        # Sort by y-axis (e.g., branch length) to get better z-axis consistency.
        lines.df <- lines.df[order(lines.df$yend, lines.df$xend),]
      }
      if (aes.type == 'label') {
        if (!plot.internals) {
          labels.df <- subset(labels.df, is.leaf==TRUE)
        }

        # Add a bit of margin to the labels
        text.margin <- 1 / 50
        mrgn <- text.margin * max.length
        
        if (layout == 'radial') {
          labels.df$y <- labels.df$y + mrgn * sin(labels.df$angle / 180 * pi)
          labels.df$x <- labels.df$x + mrgn * cos(labels.df$angle / 180 * pi)
          aes.list[['angle']] <- 'pmax((y/(diff(range(y))+1) * -360) %% 360 + 90, (y/(diff(range(y))+1) * -360 - 180) %% 360 + 90)'
          aes.list[['hjust']] <- 'ifelse( (-360 * (y)/(diff(range(y))+1)) %% 360 + 90 > (-360 * (y)/(diff(range(y))+1) - 180) %% 360 + 90, 0, 1);'
        } else {
          labels.df$hjust <- ifelse((labels.df$angle-90) %% 360 > (labels.df$angle-90+180) %% 360, 0, 1)
          labels.df$angle <- pmax((labels.df$angle-90) %% 360, (labels.df$angle-90+180) %% 360) + 90
          labels.df$y <- labels.df$y + mrgn * sin(labels.df$angle / 180 * pi) * ifelse(labels.df$hjust == 1, -1, 1)
          labels.df$x <- labels.df$x + mrgn * cos(labels.df$angle / 180 * pi) * ifelse(labels.df$hjust == 1, -1, 1)
          aes.list[['angle']] <- 'angle'
          aes.list[['hjust']] <- 'hjust'
        }
      }

      # Add an aesthetic entry when we have the appropriate argument, e.g. node.color.by='pop.size' => color='pop.size'
      for (sub.aes in c('color', 'alpha', 'size')) {
        scale.key <- ifelse(sub.aes == 'color', 'colour', sub.aes)
        cur.s <- paste(aes.type, sub.aes, 'by', sep='.')
        if (!is.null(args[cur.s])) {
          aes.list[[scale.key]] <- args[[cur.s]]
        }
      }

      geom.list <- switch(aes.type,
        'line'=list(
          data=lines.df,
          do.call(aes_string, aes.list)
        ),
        'node'=list(
          data=nodes.df,
          do.call(aes_string, aes.list)
        ),
        'label'=list(
          data=labels.df,
          vjust=0.5,
          do.call(aes_string, aes.list)
        )
      )      

      if (is.null(aes.list[['hjust']]) && aes.type == 'label') {
        geom.list[['hjust']] <- 0
      }

      # Default values when the visual properties are constants.
      constant.defaults <- switch(aes.type,
        'line'=list(
          color='#888888',
          size=1,
          alpha=1
        ),
        'node'=list(
          color='#333333',
          size=1.5,
          alpha=1
        ),
        'label'=list(
          color='black',
          size=3,
          alpha=1
        )
      )

      # Add constant values if there is no aesthetic entry, e.g. node.color='black' => color='black'
      for (sub.aes in c('color', 'alpha', 'size')) {
        constant.key <- ifelse(sub.aes == 'color', 'colour', sub.aes)
        constant.s <- paste(aes.type, sub.aes, sep='.')
        by.s <- paste(aes.type, sub.aes, 'by', sep='.')
        # If we don't have an aesthetic
        if (is.null(args[[by.s]])) {
          # If we don't have a user-specified constant
          if (is.null(args[[constant.s]])) {
            # Then use the default value
            geom.list[[constant.key]] <- constant.defaults[[sub.aes]]
          } else {
            # We have a user-specified constant -- use it
            geom.list[[constant.key]] <- args[[constant.s]]
          }
        }
      }
      
      geom.fn <- switch(aes.type,
        line='geom_joinedsegment',
        node='geom_point',
        label='geom_text'
      )
      q <- q + do.call(geom.fn, geom.list)
      
      # Add the scale indicators, e.g. node.color.scale=scale_color_gradient(low='blue', hi='red')
      for (sub.aes in c('color', 'alpha', 'size')) {
        scale.s <- paste(aes.type, sub.aes, 'scale', sep='.')
        if (!is.null(args[scale.s])) {
          q <- q + args[[scale.s]]
        }
      }
    }
  }

  if (!is.na(x.lim)) {
    bl.limits <- x.lim
  } else {
    bl.limits <- range(lines.df$x)
  }

  if (!is.na(y.lim)) {
    leaf.limits <- y.lim
  } else {
    leaf.limits <- range(lines.df$y)
  }

  if (any(!is.na(x.expand))) {
    if (layout == 'unrooted') {    
      range <- leaf.limits[2] - leaf.limits[1]
      leaf.limits[1] <- leaf.limits[1] - range*x.expand[1]
      leaf.limits[2] <- leaf.limits[2] + range*x.expand[2]
    }
    range <- bl.limits[2] - bl.limits[1]
    bl.limits[1] <- bl.limits[1] - range*x.expand[1]
    bl.limits[2] <- bl.limits[2] + range*x.expand[2]
  }

  if (layout == 'radial') {
    q <- q + scale_y_continuous(limits=c(0, max.n), name='')
    q <- q + scale_x_continuous(name=x.lab)
    q <- q + coord_polar(expand=F, theta="y")
    q <- q + opts(axis.text.x=theme_blank())
  } else {
    #q <- q + coord_flip()
    q <- q + scale_x_continuous(name=x.lab, limits=bl.limits)
    y.lab <- ifelse(layout == 'unrooted', x.lab, '')
    q <- q + scale_y_continuous(name=y.lab, limits=leaf.limits)
  }

  if (nrow(trees.df) > 1) {
    q <- q + facet_wrap( ~ .id)
  }

#  q <- q + theme_phylo_default()

  if (do.plot) {
    print(q)
  }
  q
}

is.standard.layout <- function(x) {
  any(x %in% c('default', 'radial'))
}

#' Returns a data frame defining segments to draw the phylogenetic tree.
#' 
#' @export
tree.layout <- function(
  phylo,
  layout = 'default',
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
                   angle=0,
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

  if (is.standard.layout(layout)) {
    # For each leaf: travel up towards the root, laying out each internal node along the way.
    for (i in 1:n.leaves) {
      cur.node <- i
      while (length(cur.node) > 0 && cur.node != -1) {
        df[cur.node, 'angle'] <- 0

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
  }

  if (layout == 'unrooted') {
    # See http://code.google.com/p/phylowidget/source/browse/trunk/PhyloWidget/src/org/phylowidget/render/LayoutUnrooted.java
    # For unrooted layout, we start from the root.
    layout.f <- function(node, lo, hi) {
      cur.enclosed <- tree.leaves.beneath(phylo, node)
      cur.x <- df[node, 'x']
      cur.y <- df[node, 'y']

      children <- unlist(df[node, ]$children)
      if (length(children) > 0 && children[1] != -1) {
        cur.angle <- lo
        for (i in 1:length(children)) {
          child.node <- children[i]
          child.enclosed <- tree.leaves.beneath(phylo, child.node)
          child.ratio <- child.enclosed / cur.enclosed
          child.bl <- tree.branch.length(phylo, child.node)

          arc.size <- (hi - lo) * child.ratio
          mid.angle <- cur.angle + arc.size / 2
          
          child.x <- cur.x + cos(mid.angle) * child.bl
          child.y <- cur.y + sin(mid.angle) * child.bl

          df[child.node, 'x'] <<- child.x
          df[child.node, 'y'] <<- child.y
          df[child.node, 'angle'] <<- mid.angle / (2*pi) * 360

          layout.f(child.node, cur.angle, cur.angle+arc.size)
          cur.angle <- cur.angle + arc.size          
        }        
      }      
    }
    
    layout.f(tree.get.root(phylo), 0, 2 * pi)
  }

  df$dir <- 'none'

  # We have a data frame with each node positioned.
  # Now we go through and make two line segments for each node (for a 'square corner' type tree plot).
  line.df <- data.frame()
  for (i in 1:nrow(df)) {
    row <- df[i,]            # Data frame row for the current node.
    if (row$parent == -1) {
      next; # Root node!
    }
    p.row <- df[row$parent,] # Data frame row for the parent node.

    if (is.standard.layout(layout) && !(layout.ancestors && found.any.internal.node.sequences)) {
      horiz.line <- data.frame(
                               node=row$node,
                               y=row$y,
                               yend=row$y,
                               x=row$x,
                               xend=p.row$x,
                               label=row$label,   
                               dir='up',
                               branch.length=row$branch.length
                               )
      vert.line <- data.frame(
                               node=row$node,
                               y=row$y,
                               yend=p.row$y,
                               x=p.row$x,
                               xend=p.row$x,
                               label=row$label,
                               dir='across',
                               branch.length=row$branch.length
      )
      line.df <- rbind(line.df, horiz.line, vert.line)
    } else {
      up.line <- data.frame(
                               node=row$node,
                               y=row$y,
                               yend=p.row$y,
                               x=row$x,
                               xend=p.row$x,
                               label=row$label,
                               dir='up',
                               branch.length=row$branch.length
                               )
      line.df <- rbind(line.df, up.line)
    }
  }

  line.df <- tags.into.df(phylo, line.df)
  df <- tags.into.df(phylo, df)
  # Remove weird list-of-lists from the nodes data frame.
  df$children <- NULL

  label.df <- df
  line.df$type <- 'line'
  df$type <- 'node'
  label.df$type <- 'label'

  all.df <- rbind.fill(line.df, df, label.df)
  all.df
}

#' Given a \code{\link{phylo}} object and a data frame, transform all
#' the tags from the tree ' into columns of the data frame. Rows of the
#' data frame are linked to the tree via a required 'node' column, which
#' must contain integer indices of the associated node.
#'
#' This function is similar to the tree.as.data.frame method, but not
#' exactly the same. It is used internally by the tree.layout
#' function.
#' 
#' @param phylo, input phylo object
#' @param df, data.frame with a 'node' column corresponding to integer indices
#' of tree nodes.
#' @return df, a copy of the input data frame, with tags added as new columns
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

#' A white-on-black theme for futuristic plots.
#'
#' @export
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
