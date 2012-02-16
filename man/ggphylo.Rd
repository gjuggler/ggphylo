\docType{package}
\name{ggphylo}
\alias{ggphylo}
\alias{ggphylo-package}
\alias{package-ggphylo}
\title{Plots a tree or a list of trees using \code{\link{ggplot}}.}
\usage{
  ggphylo(x, extra.data = NA, layout = "default",
    do.plot = T, plot.internals = F,
    x.lab = "Branch Length", x.lim = NA,
    x.expand = c(0.05, 0.3), y.lim = NA, ...)
}
\arguments{
  \item{x}{input phylo object, list of phylo objects, or a
  Newick- or NHX-format tree string.  When a list is given,
  all trees will be arranged using \code{\link{facet_wrap}}
  into a grid of plots. Trees of different shapes and sizes
  can be included; the plotting area will be expanded to
  fit the largest input tree.}

  \item{extra.data}{an optional data.frame or string
  pointing to a CSV file, containing data to be linked to
  the input phylo object(s). The input data must contain a
  column named 'label', which will be used to match each
  row in the data to a node in the input tree(s). Defaults
  to NA.}

  \item{layout}{a string indicating how the tree will be
  laid out. Defaults to 'defualt'.  Available layouts
  include: \item{default} The 'default' layout, which
  arranges the tree in the standard 'phylogram' arrangement
  with branch length along the x-axis and tree nodes evenly
  spaced along the y-axis.  \item{unrooted} The 'unrooted'
  layout, which arranges the tree from root-to-tip giving
  equal viewing angle to each clade according to its size.
  \item{radial} The 'radial' layout, which is equivalent to
  the 'default' layout but plotted in polar coordinates.}

  \item{do.plot}{a boolean indicating whether the plot
  should be printed to the current graphics device. When
  FALSE, the plot is not printed; the graphics grob is
  silently returned.  Defaults to TRUE.}

  \item{plot.internals}{boolean, whether the labels of
  internal nodes should be shown. Defaults to FALSE.}

  \item{x.lab}{string, the label given to the x-axis (or,
  when layout='unrooted', given to both axes). Defaults to
  'Branch Length'.}

  \item{x.lim}{vector of length 2, the lower and upper
  limits to apply to the x-axis. Defaults to NA, in which
  case the plot is expanded to include the entire tree.}

  \item{y.lim}{vector of length 2, the lower and upper
  limits to apply to the y-axis. Defaults to NA, in which
  case the plot is expanded to include the entire tree.}

  \item{x.expand}{vector of length 2, the fraction by which
  to expand the x-axis limits to the left and to the right.
  This is useful to allow extra space for trees with long
  labels or large node sizes. Defaults to c(0.05, 0.3)
  which extends the x-axis limits 5% to the left and 30% to
  the right of the default size.}

  \item{[line|node|label].[size|color|alpha].by}{string,
  indicating the name of a tag value by which to modulate
  the given visual property. A value of 0 will be given to
  nodes which do not contain the given tag, and unless a
  'x.y.scale' parameter is given (see below), the default
  \ggplot scaling and size/color/alpha rules will be
  applied.  The combination of 3 visual elements and 3
  visual properties creates 9 possible parameters to choose
  from: \item{line.size.by} \item{line.color.by}
  \item{line.alpha.by} \item{node.size.by}
  \item{node.color.by} \item{node.alpha.by}
  \item{label.size.by} \item{label.color.by}
  \item{label.alpha.by}}

  \item{[line|node|label].[size|color|alpha].scale}{function,
  used in conjunction with a corresponding 'x.y.by'
  parameter (e.g. node.scale.by) to specify the type,
  limits and range of the visual scaling. The value is
  usually the result of one of the \code{\link{ggplot}}
  scale convenience functions, such as:
  \item{node.size.scale=scale_size_continuous(limits=c(0,100),
  range=c(1, 5))} In this example, 'limits' controls the
  range of tag values to be shown, and 'range' controls the
  range of the resulting visual scale (i.e. the node size
  will range from 1 to 5).  See
  \code{\link{scale_size_continuous}} or the other ggplot
  scale_x_y functions for more info.}

  \item{[line|node|label].[size|color|alpha]}{string or
  numeric, used to specify a constant value for a visual
  parameter if no data-based scaling is specified. Default
  values are: \item{line.size} 1 \item{line.color}
  '#888888' \item{line.alpha} 1 \item{node.size} 1.5
  \item{node.color} '#333333' \item{node.alpha} 1
  \item{label.size} 3 \item{label.color} 'black'
  \item{label.alpha} 1}
}
\value{
  the \code{\link{ggplot}} grob list.
}
\description{
  Plots a tree or a list of trees using
  \code{\link{ggplot}}.

  \code{ggphylo} provides convenient functions and tools
  for visualizing phylogenetic data in R. Many methods are
  provided to simplify the process of working with
  \code{phylo} objects (e.g., \code{\link{tree.scale.to}}
  and \code{\link{tree.translate}} , while two plotting
  methods allow for a wide range of expressive plots:
  \code{\link{tree.plot}} and \code{\link{aln.plot}}
}

