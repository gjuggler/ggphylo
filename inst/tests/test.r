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

  # Test find the root.
  expect_equal(tree.get.root(tree), 4)

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

  # Scaling.
  str <- '((a,b),c);'
  tree <- tree.read(str)
  tree <- tree.scale.by(tree, 3)

})

context("tree plots")
test_that("Tree plots work", {
  pdf(file="~/scratch/test.pdf")

  tree <- tree.read('((a:1,b:1):1,c:1);') # Read in a simple tree.
  tree.list <- list()
  for (i in 1:10) {
    tree.list[[i]] <- tree.scale.by(tree, i)
  }
  ggphylo(tree.list)

  n <- 40;
  x <- rtree(n)
  bootstrap <- runif(length(nodes(x)), min=0, max=100)
  values2 <- runif(length(nodes(x)), min=0, max=5)
  for (i in nodes(x)) {
    x <- tree.set.tag(x, i, 'bootstrap', bootstrap[i])
    x <- tree.set.tag(x, i, 'v2', values2[i])
  }

  ggphylo(x, layout='unrooted', node.size=2, x.expand=c(0.2, 0.2))
  ggphylo(x, layout='radial', plot.nodes=F, line.color.by='bootstrap')

  x <- tree.read('(a,(b,(c,(d,(e,(f,(g,(h,(i,(j,(k,(l,(m,(n,(o,(p,(q,(r,(s,(t,(u,(v,(w,(x,(y,z)))))))))))))))))))))))));')
  x <- tree.normalize.branchlengths(x)
  ggphylo(x, layout='radial')

  x <- tree.read('((((Picrophilus_torridus:4.7126427,Feroplasma_acidarmanus:4.7126427):1.605092,(Thermoplasma_volcanium:1.576355,Thermoplasma_acidophilum:1.576355):4.7413797):6.5871563,((((Methanocorpusculum_labreanum:6.937603,(Methanoculleus_marisnigri:5.763546,Methanospirillum_hungatei:5.763546):1.174057):2.0563765,(Methanosaeta_thermophila:7.8448277,((Methanosarcina_barkeri:1.5681438,(Methanosarcina_mazei:1.3957291,Methanosarcina_acetivorans:1.3957291):0.17241478):2.8708262,Methanococcoides_burtonii:4.43897):3.4058576):1.1491518):2.4838524,((Haloquadratum_walsbyi:5.5665016,Haloferax_volcanii:5.5665016):1.0125904,(Haloarcula_marismortui:5.853859,(Natronomonas_pharaonis:4.860428,Halobacterium_halobacterium:4.860428):0.9934311):0.7252331):4.89874):1.0009365,((Methanococcus_maripaludis:10.812807,(Methanosphaera_stadtmanae:9.54023,Methanobacterium_thermoautotrophicum:9.54023):1.2725773):1.4293184,(Nanoarchaeum_equitans:10.499648,((Metallosphaera_sedula:6.152162,(Sulfolobus_solfataricus:5.188835,(Sulfolobus_tokodaii:4.2528725,Sulfolobus_acidocaldarius:4.2528725):0.9359627):0.96332693):3.7331367,(((((Thermococcus_kodakaraensis:2.5396817,((Pyrococcus_furiosus:0.14778137,Pyrococcus_horikoshii:0.14778137):0.09031296,Pyrococcus_abyssi:0.23809433):2.3015873):4.303505,Methanopyrus_kandleri:6.8431864):0.24712515,Methanocaldococcus_jannaschii:7.0903115):1.0842919,Archaeoglobus_fulgidus:8.174603):0.7848158,((Thermofilum_pendens:6.4422555,(Staphylothermus_marinus:3.8095245,(Hyperthermus_butylicus:3.1198692,Aeropyrum_pernix:3.1198692):0.6896553):2.632731):0.8470168,(Pyrobaculum_calidifontis:6.1970444,(Caldivirga_maquilingensis:4.95895,(Thermoproteus_tenax:2.665572,(Pyrobaculum_arsenaticum:0.8702793,(Pyrobaculum_islandicum:0.7389145,Pyrobaculum_aerophilum:0.7389145):0.13136482):1.7952926):2.293378):1.2380943):1.0922279):1.670147):0.9258795):0.61434937):1.7424774):0.23664284):0.42612267):0.42609692,Cenarchaeum_symbiosum:13.330988);')
  ggphylo(x, layout='unrooted', x.expand=c(0.5, 0.5))

  x <- tree.read('(((BGIOSIBCA042311_ORYSA:0.11893[&&NHX:S=ORYSA:O=BGIOSIBCA042311.1:G=BGIOSIBCA042311],At3g11830_ARATH:0.103292[&&NHX:S=ARATH:O=At3g11830.1:G=At3g11830]):0.0816687[&&NHX:S=Magnoliophyta:B=99:D=N],(cct7_SCHPO:0.166585[&&NHX:S=SCHPO:O=SPBC25H2.12c:G=SPBC25H2.12c],CCT7_YEAST:0.186209[&&NHX:S=YEAST:O=YJL111W.1:G=YJL111W]):0.029676[&&NHX:S=Ascomycota:B=99:D=N]):0.00495409[&&NHX:S=Eukaryota:B=55:D=N],((((((((CCT7_HUMAN:0.0510845[&&NHX:S=HUMAN:O=ENST00000258091.4:G=ENSG00000135624],Cct7_MOUSE:0.0548045[&&NHX:S=MOUSE:O=ENSMUST00000032078.4:G=ENSMUSG00000030007]):0.0267848[&&NHX:S=Euarchontoglires:B=54:D=N],ENSMODT00000021123_MONDO:0.0837789[&&NHX:S=MONDO:O=ENSMODT00000021123.2:G=ENSMODG00000016626]):0.0255967[&&NHX:S=Theria:B=54:D=N],ENSGALT00000025922_CHICK:0.131991[&&NHX:S=CHICK:O=ENSGALT00000025922.2:G=ENSGALG00000016088]):0.0021982[&&NHX:S=Amniota:B=53:D=N],Q5D024_XENTR:0.105271[&&NHX:S=XENTR:O=ENSXETT00000025763.2:G=ENSXETG00000011778])[&&NHX:S=Tetrapoda:B=96:D=N],(cct7_BRARE:0.0999995[&&NHX:S=BRARE:O=ENSDART00000006802.3:G=ENSDARG00000007385],NEWSINFRUT00000128352_FUGRU:0.108658[&&NHX:S=FUGRU:O=NEWSINFRUT00000128352.3:G=NEWSINFRUG00000121580]):0.00780083[&&NHX:S=Clupeocephala:B=100:D=N]):0.0510091[&&NHX:S=Euteleostomi:B=97:D=N],(ENSCINT00000001715_CIOIN:0.00425529[&&NHX:S=CIOIN:O=ENSCINT00000001715.2:G=ENSCING00000000947],ENSCINT00000003265_CIOIN:0.0104703[&&NHX:S=CIOIN:O=ENSCINT00000003265.2:G=ENSCING00000001620]):0.139761[&&NHX:S=CIOIN:B=100:D=Y]):0.0066056[&&NHX:S=Chordata:B=91:D=N],(Sm00.scaff00161.0170_SCHMA:0.187113[&&NHX:S=SCHMA:O=Sm00.scaff00161.0170:G=Sm00.scaff00161.0170],(NP_503522_CAEEL:0.0698092[&&NHX:S=CAEEL:O=T10B5.5a.1:G=T10B5.5],CBG21817_CAEBR:0.0734101[&&NHX:S=CAEBR:O=CBG21817:G=CBG21817]):0.142116[&&NHX:S=Caenorhabditis:B=94:D=N])[&&NHX:S=Bilateria:B=35:D=N])[&&NHX:S=Bilateria:D=N],(NP_649835_DROME:0.0835465[&&NHX:S=DROME:O=CG8351-RA.3:G=CG8351],GA21011-PA_DROPS:0.0840653[&&NHX:S=DROPS:O=GA21011-PA:G=GA21011]):0.108617[&&NHX:S=Sophophora:B=100:D=N]):0.00495409[&&NHX:S=Bilateria:D=N])[&&NHX:S=Eukaryota:D=N];')
  x <- tree.normalize.branchlengths(x, push.to.tips=T)
  ggphylo(x, node.color.by='D', x.expand=c(0, 1))

  dev.off()

})