##
## Prepare feature matrices for use within R
##

library(wordspace)

fetch.data <- function (filename) {
  data.mat <- read.csv(filename, fileEncoding="utf8", check.names=FALSE)
  words <- data.mat[, 1]
  data.mat <- t(as.matrix(data.mat[, -1]))
  colnames(data.mat) <- words
  files <- sub("\\.txt$", "", rownames(data.mat), perl=TRUE)
  fullnames <- gsub("-", " ", sapply(strsplit(files, split="_", fixed=TRUE), function (x) x[1]))
  authors <- gsub(",.*$", "", fullnames, perl=TRUE)
  titles <- sapply(strsplit(files, split="_", fixed=TRUE), function (x) x[2])
  rowlabels <- paste(authors, titles, sep=": ")
  rownames(data.mat) <- rowlabels
  titleinfo <- data.frame(term=rowlabels, author=authors, fullname=fullnames, title=titles, stringsAsFactors=FALSE)
  res <- dsm(M=data.mat, rowinfo=titleinfo, raw.freq=TRUE)
  dsm.score(res, score="frequency", normalize=TRUE, method="manhattan") # compute relative frequencies as scores
}

FreqDE <- fetch.data("data/delta_corpus_DE.csv")
FreqEN <- fetch.data("data/delta_corpus_EN.csv")
FreqFR <- fetch.data("data/delta_corpus_FR.csv")

## plausibility check
dim(FreqDE)
table(FreqDE$rows$author)
dim(FreqFR)
table(FreqFR$rows$author)
dim(FreqEN)
table(FreqEN$rows$author)

## use re-weighting to standardize relative frequencies
zDE <- dsm.score(FreqDE, score="reweight", scale="standardize", matrix.only=TRUE)
zEN <- dsm.score(FreqEN, score="reweight", scale="standardize", matrix.only=TRUE)
zFR <- dsm.score(FreqFR, score="reweight", scale="standardize", matrix.only=TRUE)

## extract gold standard labels
goldDE <- FreqDE$rows$author
goldEN <- FreqEN$rows$author
goldFR <- FreqFR$rows$author

## save in .rda format
save(FreqDE, FreqEN, FreqFR, zDE, zEN, zFR, goldDE, goldFR, goldEN, file="data/delta_corpus.rda", compress="xz")
