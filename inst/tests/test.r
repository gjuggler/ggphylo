# tree
set.seed(1)

context("tree to data.frame")
test_that("as.data.frame() works as expected", {
  n.tips <- 10
  tree <- rcoal(n.tips)
                           
  tree.df <- as.data.frame(tree)
  expect_that(tree.df, is_a("data.frame"))

  # One row for each tree node
  expect_equal(nrow(tree.df), (n.tips + n.tips-1))

  # Subsets for leaves and nodes
  expect_equal(nrow(subset(tree.df, is.leaf == T)), n.tips)
  expect_equal(nrow(subset(tree.df, is.leaf == F)), n.tips-1)

  expect_equal(min(tree.df$depth), 1)
  expect_equal(max(tree.df$depth), n.tips)

  # Use a simple tree to test node ordering.
  tree.str <- "((a,b)ab,c)abc;"
  tree <- tree.read(tree.str)

  # Ordering visually lays out nodes according to how they would have been plotted
  x <- as.data.frame(tree)
  expect_equal(paste(x$label, collapse=' '), 'a ab b abc c')

  # Not ordering visually gives the normal cladewise order
  x <- as.data.frame(tree, order.visually=F)
  expect_equal(paste(x$label, collapse=' '), 'a b c abc ab')
})

context("tree I/O")
test_that("tree I/O and NHX annotations work correctly", {
  # tree.read gives a warning when it detects NHX strings, and passes the input to tree.read.nhx
  tree.str <- "((a[&&NHX:foo=bar],b),c);"
  expect_warning(tree.read(tree.str))

  # tree.read.nhx expects NHX input.
  tree <- tree.read.nhx(tree.str)

  # NHX tags are turned into columns of the data frame, and set to NA for nodes that didn't have
  # that tag.
  tree.df <- as.data.frame(tree)
  expect_equal(nrow(subset(tree.df, !is.na(foo))), 1)


  # Let's try some round trip tests.
  tree.roundtrip <- function(str) {
    tree <- tree.read.nhx(str)
    expect_equal(str, as.character(tree))
  }
  tree.roundtrip('((a,b),c);')
  tree.roundtrip('((a,b,c,d,e,f,g),h);')
  tree.roundtrip('((a,b:2,c,d:3,e,f:4,g),h:5);')
  tree.roundtrip('((a[&&NHX:foo=bar],b[&&NHX:bizz=buzz]),c[&&NHX:ping=pong]);')
  tree.roundtrip('((a[&&NHX:long key with spaces=blah blah blah... blah blah!!],b),c);')

})

context("tree manipulations")
test_that("Tree manipulations work", {
  set.seed(1)
  tree <- rcoal(3)

  # Test the foreach
  str <- ''
  tree.foreach(tree, function(x, i) {
    str <<- paste(str, tree.get.label(x, i), sep='')
  })
  expect_equal(str, 't2t3t1')

  # Test normalizing branch lengths.
  str <- '((a,b)e:4,(c:2,d:1)f)g;'
  tree <- tree.read(str)

  a <- tree.normalize.branchlengths(tree, push.to.tips=F)
  b <- tree.normalize.branchlengths(tree, push.to.tips=T)
  c <- tree.strip.branchlengths(tree)

  # Pushing clades to the tips should create longer branchlength for internal clades.
  fn <- function(x) {
    tree.branch.length(x, tree.find(x, 'e'))
  }
  expect_true(fn(tree) > fn(b))
  expect_true(fn(b) > fn(a))
  expect_true(is.na(fn(c)))

  # Try translating leaf and node labels
  str <- '((a,b),c)d;'
  tree <- tree.read(str)
  map <- list(
    a = 'AAAAAAmazing!',
    b = 'Boo-ya',
    c = 'Cowabunga!',
    d = 'Duuude...'
  )
  tree <- tree.translate(tree, map)
  expect_true(length(tree.find(tree, 'Boo-ya')) > 0)
  expect_true(length(tree.find(tree, 'Cowabunga!')) > 0)
  expect_true(length(tree.find(tree, 'Duuude...')) > 0)

  # Try adding a tag to a tree.
  str <- '((a,b),c)d;'
  tree <- tree.read(str)
  tree <- tree.set.tag(tree, tree.find(tree, 'b'), 'foo', 'bar')
  expect_true(tree.has.tags(tree))
  expect_true(tree.get.tag(tree, tree.find(tree, 'b'), 'foo') == 'bar')

})

context("tree plots")
test_that("Tree plots work", {

})