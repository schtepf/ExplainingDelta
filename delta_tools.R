## ----DT_libraries, message=FALSE----------------------------------------------
library(wordspace)
library(cluster)      # for PAM clustering
library(MASS)         # for MDS visualization
library(matrixStats)  # fast row- and column-wise statistics
library(colorspace)   # for well-designed colour palettes


## ----DT_globals---------------------------------------------------------------
## standard colour palettes of Seaborn visualization library
seaborn.pal <- c("#4C72B0","#55A868","#C44E52","#8172B2","#CCB974","#64B5CD")
muted.pal <- c("#4878CF","#6ACC65","#D65F5F","#B47CC7","#C4AD66","#77BEDB")
bright.pal <- c("#003FFF","#03ED3A","#E8000B","#8A2BE2","#FFC400","#00D7FF")
grayscale.pal <- c("black", "#888888", "#555555", "#BBBBBB")


## ----DT_adjustedRandIndex-----------------------------------------------------
adjustedRandIndex <- function (x, y) {
  x <- as.vector(x)
  y <- as.vector(y)
  if (length(x) != length(y)) stop("arguments must be vectors of the same length")
  tab <- table(x, y)
  if (all(dim(tab) == c(1, 1))) return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d)) / 
         ((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  ARI
}


## ----DT_delta.dist------------------------------------------------------------
delta.dist <- function (M, n=NA, method="euclidean", p=2, normalize=NA, norm.p=p,
                        pca=NA, whiten=TRUE, prenorm=FALSE, transform=NULL) {
  if (!is.na(n) && n < ncol(M)) M <- M[, 1:n]
  if (!is.na(pca)) {
    if (prenorm) M <- normalize.rows(M, method="euclidean")
    res <- prcomp(M, center=TRUE)
    M <- if (whiten) scaleMargins(res$x, cols=1 / res$sdev) else res$x
    if (pca < ncol(M)) M <- M[, 1:pca]
  }
  if (!is.na(normalize)) M <- normalize.rows(M, method=normalize, p=norm.p)
  res <- dist.matrix(M, method=method, p=p, convert=TRUE)
  if (!is.null(transform)) res <- transform(res)
  res
}


## ----ternarize----------------------------------------------------------------
ternarize <- function(x, neutral.p=0, neutral.range=-qnorm((1 - neutral.p)/2), crossover=NULL) {
  y <- ifelse(abs(x) <= neutral.range, 0, sign(x))
  if (!is.null(crossover)) {
    if (length(dim(x)) != 2) stop("argument must be a (sparse or dense) matrix if crossover= is specified")
    n <- seq_len(ncol(x))
    q <- ifelse(n > crossover, 1, n / crossover)
    scaleMargins(y, cols=q) + scaleMargins(x, cols=1-q)
  } else {
    y
  }
}


## ----binarize-----------------------------------------------------------------
binarize <- function(x, threshold=0, underuse.val=0, crossover=NULL) {
  y <- ifelse(x <= threshold, underuse.val, 1)
  if (!is.null(crossover)) {
    if (length(dim(x)) != 2) stop("argument must be a (sparse or dense) matrix if crossover= is specified")
    n <- seq_len(ncol(x))
    q <- ifelse(n > crossover, 1, n / crossover)
    scaleMargins(y, cols=q) + scaleMargins(x, cols=1-q)
  } else {
    y
  }
}


## ----binternize---------------------------------------------------------------
binternize <- function(z, f, df=nrow(z)/2, hapax=FALSE, crossover=NULL,
                       neutral.p=1/3, neutral.range=-qnorm((1 - neutral.p)/2)) {
  if (length(dim(z)) != 2) stop("first argument must be a (sparse or dense) matrix")
  if (length(dim(f)) != 2) stop("second argument must be a (sparse or dense) matrix")
  if (!all(dim(z) == dim(f))) stop("the two arguments must be matrices of the same dimensions")
  y3 <- ternarize(z, neutral.range=neutral.range) # ternarized matrix
  y2 <- ifelse(f > 0, 1, 0)                       # binarized matrix
  df.vec <- colSums(y2)                           # document frequencies of all features
  if (hapax) y2 <- ifelse(f > 1, 1, 0)
  y <- scaleMargins(y3, cols=(df.vec >= df)) + scaleMargins(y2, cols=(df.vec < df))
  if (!is.null(crossover)) {
    n <- seq_len(ncol(y))
    q <- ifelse(n > crossover, 1, n / crossover)
    scaleMargins(y, cols=q) + scaleMargins(z, cols=1-q)
  } else {
    y
  }
}


## ----clamp--------------------------------------------------------------------
clamp <- function (x, min=-1, max=1) {
  pmax(pmin(x, max), min)  
}


## ----spike_plot---------------------------------------------------------------
spike.pal <- rainbow_hcl(6, c=60, l=50)[c(1,5,3,2,6,4)] # suitable palette with clear contrasts
spike.plot <- function (x, lwd=1, col=spike.pal, lty="solid", stride=NA,
                        diff=FALSE, legend=NULL, ylim=range(x), ...) {  
  if (is.matrix(x)) {
    k <- nrow(x)
    xs <- lapply(1:k, function (i) x[i, ])
  } else {
    k <- 1
    xs <- list(x)
  }
  n <- length(xs[[1]])
  if (diff) {
    if (k != 2) stop("x must be a matrix with two rows if diff=TRUE")
  } else {
    col <- rep(col, length.out=k) # recycle colours / line types and adjust length
    lty <- rep(lty, length.out=k)
  }
  if (is.na(stride)) stride <- if (diff) 1 else k + 1
  x.vals <- stride * ((1:n) - 1) # coordinates of gap before first spike in each group
  plot(0, 0, type="n", xaxs="i", xaxt="n", xlim=c(0, stride * n + 1), xlab="",
       ylim=ylim, ylab="", ...)
  if (diff) {
    ## this is legacy code and may not respect all parameter settings
    y1 <- xs[[1]]
    y2 <- xs[[2]]
    col.vec <- ifelse(y2 > y1, col[3], col[1])
    segments(x.vals + 1, y1, x.vals + 1, y2, col=col.vec, lwd=lwd, lend=1)
    abline(h = 0)
    if (!is.null(legend)) {
      legend("topright", inset=.02, bg="white", legend=sprintf("%s %s %s", legend[2], c(">", "<"), legend[1]), col=col[c(3, 1)], lwd=lwd + 1)
    }
  } else {
    for (i in 1:k) {
      segments(x.vals + i, rep(0, n), x.vals + i, xs[[i]], col=col[i], lty=lty[i], lwd=lwd, lend=1)
    }
    abline(h=0)
    if (!is.null(legend)) {
      legend("topright", inset=.02, bg="white", legend=legend, col=col, lty=lty, lwd=lwd + 1)
    }
  }
}
# spike.plot(zDE[1:2,1:50], diff=TRUE, lwd=3, legend=c("A", "B"), ylim=c(-2.5, 2.5))  # for testing


## ----DT_nn.classify-----------------------------------------------------------
nn.classify <- function (M, gold, ..., predicted=FALSE) {
  stopifnot(nrow(M) == length(gold))
  gold <- as.character(gold) # in case gold is a factor
  if (is.null(names(gold))) names(gold) <- rownames(M)
  DM <- delta.dist(M, ...)
  ## nearest neighbour of each vector = result of NN classifier in leave-one-out CV
  neighbours <- nearest.neighbours(DM, rownames(DM), n=1, method=method)
  nn.predict <- gold[sapply(neighbours, names)] # classifier predictions
  accuracy <- 100 * sum(nn.predict == gold) / length(gold)
  if (predicted) nn.predict else accuracy
}


## ----DT_pam.cluster-----------------------------------------------------------
pam.cluster <- function (M, gold, ..., clusters=NULL, predicted=FALSE,
                         clust.method=c("pam", "ward", "complete", "single", "average", "hclust")) {
  stopifnot(nrow(M) == length(gold))
  clust.method <- match.arg(clust.method)
  gold <- as.character(gold) # in case gold is a factor
  if (is.null(clusters)) clusters <- length(unique(gold))
  if (is.null(names(gold))) names(gold) <- rownames(M)
  DM <- delta.dist(M, ...)
  if (clust.method == "pam") {
    if (length(clusters) > 1) {
      clustering.list <- lapply(clusters, function (k) pam(DM, k, diss=TRUE, cluster.only=TRUE))
      sil.vec <- sapply(clustering.list, function (cl) summary(silhouette(cl, dmatrix=DM))$avg.width)
      idx <- which.max(sil.vec)
      clusters <- clusters[idx]
      clustering <- clustering.list[[idx]]
    } else {
      clustering <- pam(DM, clusters, diss=TRUE, cluster.only=TRUE)
    }
  } else {
    if (clust.method == "hclust") {
      hcl <- hclust(as.dist(DM), method="ward.D")
    } else {
      hcl <- agnes(DM, diss=TRUE, method=clust.method)
    }
    if (length(clusters) > 1) {
      sil.vec <- sapply(clusters, function (cl) summary(silhouette(cutree(hcl, k=cl), dmatrix=DM))$avg.width)
      clusters <- clusters[which.max(sil.vec)]
    }
    clustering <- cutree(hcl, k=clusters)
  }
  if (predicted) {
    ct <- table(gold, clustering) # gold author vs. cluster number
    labels <- rownames(ct)[apply(ct, 2, which.max)] # find majority label for each cluster number
    res <- data.frame(cluster=clustering, label=labels[clustering],
                      row.names=rownames(M), stringsAsFactors=FALSE)
  } else {
    res <- 100 * adjustedRandIndex(clustering, gold)
  }
  sil <- summary(silhouette(clustering, dmatrix=DM))$avg.width
  structure(res, n.clusters=clusters, avg.width=sil)
}


## ----DT_evaluate--------------------------------------------------------------
evaluate <- function (M, gold, n=NA, p=2, norm.p=2, clusters=NULL, label=NULL, 
                      do.nn=TRUE, do.cluster=TRUE, clust.method="pam", ...) {
  if (all(is.na(n))) n <- ncol(M)
  n.runs <- max(length(n), length(p), length(norm.p))
  n <- rep(n, length.out=n.runs)
  p <- rep(p, length.out=n.runs)
  norm.p <- rep(norm.p, length.out=n.runs)
  if (!is.null(label)) label <- sprintf("%s | n=%d", label, n)
  res <- data.frame(n=n, p=p, norm.p=norm.p, row.names=label)
  if (do.nn) {
    acc <- sapply(1:n.runs, function (i) {
      nn.classify(M, gold, n=n[i], p=p[i], norm.p=norm.p[i], ...) })
    res$accuracy <- round(acc, 2)
  }
  if (do.cluster) {
    rand <- sapply(1:n.runs, function (i) {
      pam.cluster(M, gold, n=n[i], p=p[i], norm.p=norm.p[i], 
                  clusters=clusters, clust.method=clust.method, ...) })
    res$adj.rand <- round(rand, 2)
  }
  res
}

