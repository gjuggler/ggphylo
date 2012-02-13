aln.restrict.to.tree <- function(aln, tree) {
  tip.seqs <- tree$tip.label
  aln.restrict.to.seqs(aln, tip.seqs)
}

aln.restrict.to.seqs <- function(aln, seqs) {
  aln.names <- rownames(aln)
  overlap <- intersect(aln.names, seqs)
  aln <- as.matrix(aln[which(aln.names %in% overlap),])
  aln
}

aln.remove.allgap.columns <- function(aln) {
  nongap.cols <- apply(aln, 2, function(x) {
    sum(!grepl('(---|-)', x))
  })
  
  aln[, nongap.cols > 0]
}

aln.remove.ngap.columns <- function(aln, n=2) {
  nongap.cols <- apply(aln, 2, function(x) {
    sum(!grepl('(---|-)', x))
  })
  
  aln[, nongap.cols >= n]
}

aln.remove.anygap.columns <- function(aln) {
  gap.cols <- apply(aln, 2, function(x) {
    sum(grepl('(---|-)', x))
  })
  
  aln[, gap.cols == 0]
}

aln.is.codon <- function(aln) {
  aln <- aln.tx(aln, n.char='Z')
  if (any(aln == 'X')) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

aln.tx <- function(aln, n.char='X') {
  library(phylosim)
  assign("PSIM_FAST", TRUE, envir=.GlobalEnv)
  # Translates a codon alignment into proteins.
  ca <- CodonAlphabet()

  char.length <- stringLength(aln[1,1])
  if (char.length == 1) {
    pep.aln <- matrix(nrow=nrow(aln), ncol=ncol(aln)/3)
  } else {
    pep.aln <- matrix(nrow=nrow(aln), ncol=ncol(aln))
  }
  rownames(pep.aln) <- rownames(aln)
  for (i in 1:nrow(pep.aln)) {
    codons <- paste(aln[i,],collapse='')
    for (j in 1:ncol(pep.aln)) {
      cdna.lo <- (j-1)*3+1
      cdna.hi <- (j-1)*3+3
      codon <- substr(codons, cdna.lo, cdna.hi)
      if (codon == '---') {
        pep.aln[i,j] <- '-'
      } else if (length(grep('N', codon)) > 0) {
        pep.aln[i,j] <- n.char
      } else {
        aa <- translateCodon(ca, codon)
        pep.aln[i,j] <- aa
      }
    }
  }
  return(pep.aln)
}

aln.read <- function(file) {
  read.dna(file, format='fasta')
}

s.aln.read <- function(file) {
  getTaxaNames <- function(x) {
    x <- sub("^['\" ]+", "", x)
    x <- sub("['\" ]+$", "", x)
    x
  }   
  getNucleotide <- function(x) {
    x <- gsub(" ", "", x)
    paste(x, collapse='')
  }

  X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
  start <- grep("^ {0,}>", X)
  taxa <- X[start]
  n <- length(taxa)
  obj <- vector("list", n)
  seq.names <- NULL
  if (is.null(seq.names)) {
    taxa <- sub("^ {0,}>", "", taxa)
    seq.names <- getTaxaNames(taxa)
  }   
  start <- c(start, length(X) + 1)
  for (i in 1:n) {
    obj[[i]] <- getNucleotide(X[(start[i] +1):(start[i + 1] - 1)])
  }
  names(obj) <- seq.names
  obj
}

s.aln.to.aln <- function(aln) {
  seqs <- names(aln)
  aln.len <- as.integer(nchar(aln[1]))
  aln.out <- laply(aln, function(x){unlist(strsplit(x, ''), use.names=F)})
  rownames(aln.out) <- seqs
  colnames(aln.out) <- NULL
  aln.out
}

s.aln.concat <- function(a, b) {
 a.seqs <- names(a)
 b.seqs <- names(b)
 missing.in.b <- a.seqs[!(a.seqs %in% b.seqs)]
 missing.in.a <- b.seqs[!(b.seqs %in% a.seqs)]

 # Go through and add all-gap sequences for seqs missing from 
 # one or the other alignment.

 b.length <- nchar(b[[1]])
 b.n.missing <- length(missing.in.b)
 if (b.n.missing > 0) {
   blank.str <- paste(rep('-', times=b.length), collapse='')
   for (missing.seq in missing.in.b) {
     b[missing.seq] <- blank.str
   }
 }

 a.length <- nchar(a[[1]])
 a.n.missing <- length(missing.in.a)
 if (a.n.missing > 0) {
   blank.str <- paste(rep('-', times=a.length), collapse='', sep='')
   for (missing.seq in missing.in.a) {
     a[missing.seq] <- blank.str
   }
 }

 a.seqs <- names(a)
 b.seqs <- names(b)

 aln <- list()
 for (seq in a.seqs) {
   aln[seq] <- paste(a[seq], b[seq], collapse='', sep='')
 }
 aln
}

s.aln.write <- function(aln, file) {
  # Mostly copied from ape, but I removed the space between the Fasta
  # '>' and the sequence ID.
  zz <- file(file, "w")
  on.exit(close(zz))
  colw <- 111
  N <- length(aln)
  seq.names <- names(aln)
  for (i in 1:N) {
    cur.id <- seq.names[i]
    cur.seq <- aln[[i]]
    n.chars <- nchar(cur.seq)
    cat(">", cur.id, file=zz, sep='')
    cat("\n", file=zz)
    S <- nchar(cur.seq)
    cur.seq <- unlist(strsplit(cur.seq, ''))
    cat(cur.seq, file=zz, sep='', fill=colw)
#    cat("\n", file=zz)
#    n.rows <- ceiling(S/colw)
#    for (j in 1:n.rows) {
#      lo <- 1 + (j-1) * colw
#      hi <- colw + (j-1) * colw
#      if (hi > n.chars) {
#        hi <- n.chars
#      }
#      cur.out <- substr(cur.seq, lo, hi)
#      cat(cur.out, file=zz, sep='')
#      cat("\n", file=zz)
#    }
  }
}

aln.concat <- function(a, b) {
 a.seqs <- rownames(a)
 b.seqs <- rownames(b)
 missing.in.b <- a.seqs[!(a.seqs %in% b.seqs)]
 missing.in.a <- b.seqs[!(b.seqs %in% a.seqs)]

 # Go through and add all-gap sequences for seqs missing from 
 # one or the other alignment.

 b.length <- length(b[1,])
 b.n.missing <- length(missing.in.b)
 if (b.n.missing > 0) {
   b.add <- matrix(data=rep('-', times = b.length * b.n.missing),
     nrow=b.n.missing, ncol=b.length
   )
   rownames(b.add) <- missing.in.b
   b <- rbind(b, b.add)
 }

 a.length <- length(a[1,])
 a.n.missing <- length(missing.in.a)
 if (a.n.missing > 0) {
   a.add <- matrix(data=rep('-', times = a.length * a.n.missing),
     nrow=a.n.missing, ncol=a.length
   )
   rownames(a.add) <- missing.in.a
   a <- rbind(a, a.add)
 }

 # Re-load the seq names now that we've filled in each alignment.
 a.seqs <- rownames(a)
 b.seqs <- rownames(b)

 a.match <- match(a.seqs, b.seqs)
 a.indices <- which(!is.na(a.match))
 b.indices <- a.match[!is.na(a.match)]

 aln <- cbind(
   a[a.indices,],
   b[b.indices,]
 )
 aln
}

test.aln.concat <- function() {
  a <- read.aln(text="
>a
aa
>b
bb
>d
dd
")

  b <- read.aln(text="
>c
ccc
>b
bbb
>a
aaa
>f
fff
")

  print(aln.concat(a, b))
}

aln.write <- function(x, ff) {
  # Mostly copied from ape, but I removed the space between the Fasta
  # '>' and the sequence ID.
  zz <- file(ff, "w")
  on.exit(close(zz))
  colw <- 120
  colsep <- ''
  indent <- ''
  nbcol <- 1
  N <- dim(x)
  S <- N[2]
  N <- N[1]
  xx <- vector("list", N)
  for (i in 1:N) xx[[i]] <- x[i, ]
  names(xx) <- rownames(x)
  x <- xx
  rm(xx)
  for (i in 1:N) {
    cat(">", names(x)[i], file = zz, sep='')
    cat("\n", file = zz, sep='')
    X <- paste(x[[i]], collapse = "")
    S <- length(x[[i]])
    totalcol <- ceiling(S/colw)
    if (nbcol < 0) 
      nbcol <- totalcol
    nb.lines <- ceiling(totalcol/nbcol)
    SEQ <- character(totalcol)
    for (j in 1:totalcol) {
      SEQ[j] <- substr(X, 1 + (j - 1) * colw, colw + (j - 1) * colw)
    }
    for (k in 1:nb.lines) {
      endsel <- if (k == nb.lines) 
        length(SEQ)
        else nbcol + (k - 1) * nbcol
      cat(indent, file = zz)
      cat(SEQ[(1 + (k - 1) * nbcol):endsel], sep = colsep, file = zz)
      cat("\n", file = zz)
    }
  }
}

aln.stats <- function(aln) {
  n.seqs <- length(aln[,1])
  n.cols <- length(aln[1,])

  nongap.counts <- apply(aln, 1, function(x) {
    sum(!grepl('-', x))
  })
  nz.counts <- nongap.counts != 0
  n.nonblank.seqs <- sum(nz.counts)
  mean.seq.length <- mean(nongap.counts[nz.counts])

  nuc.counts <- apply(aln, 1, function(x) {
    str <- paste(x, collapse='')
    mtchs <- unlist(gregexpr('[GACT]', str))
    sum(mtchs != -1)
  })
  total.nucs <- sum(nuc.counts)

  nongap.cols <- apply(aln, 2, function(x) {
    sum(!grepl('(---|-)', x))
  })
  n.allgaps <- sum(nongap.cols == 0)
  n.nogaps <- sum(nongap.cols == (n.nonblank.seqs))

  percent.gap.nucs <- (n.seqs * n.cols - total.nucs/3) / (n.seqs * n.cols)
  percent.allgap.columns <- n.allgaps / n.cols
  percent.nogap.columns <- n.nogaps / n.cols

  str <- sprintf('Seqs: %d
Columns: %d
All-gap columns: %d
Mean seq length: %.1f
Total nucs: %d\n',
    n.seqs,
    n.cols,
    n.allgaps,
    mean.seq.length,
    total.nucs
  )
  #cat(str)

  out.df <- data.frame(
    n_species = n.seqs,
    n_nonblank_seqs = n.nonblank.seqs,
    mean_seq_length = mean.seq.length,
    n_columns = n.cols,
    n_allgap_columns = n.allgaps,
    n_nogap_columns = n.nogaps,
    n_nucs = total.nucs,
    percent_allgap_columns = percent.allgap.columns,
    percent_nogap_columns = percent.nogap.columns,
    percent_gap_nucs = percent.gap.nucs
  )
  out.df <- format(out.df, digits=1)
  out.df
}

tree.all.labels <- function(tree) {
  c(tree$tip.label, tree$node.label)
}

aln.remove.phylosim.internals <- function(aln) {
  all.seqs <- rownames(aln)
  
  node.seqs <- grepl("^Node \\d+$", all.seqs, perl=TRUE)
  root.seqs <- grepl("^Root node \\d+$", all.seqs, perl=TRUE)

  keep.seqs <- all.seqs[!(node.seqs | root.seqs)]
  #print(keep.seqs)
  aln.restrict.to.seqs(aln, keep.seqs)
}