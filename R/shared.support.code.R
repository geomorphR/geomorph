### Functions that both RRPP and geomorph use, but should remain internal
### any alterations to these functions must be saved for both RRPP and geomorph

# get.names
# a universal subject name search
# used in lm.rrpp and its descendants
get.names <- function(Y) {
  nms <- if(is.vector(Y)) names(Y) else if(inherits(Y, "dist")) attr(Y, "Labels") else
    if(inherits(Y, "matrix")) rownames(Y) else dimnames(Y)[[3]]
  nms
}

# center
# centers a matrix faster than scale()
# used in various functions where mean-centering is required
center <- function(x){
  if(is.vector(x)) x - mean(x) else {
    x <- as.matrix(x)
    dims <- dim(x)
    fast.center(x, dims[1], dims[2])
  }
}

fast.center <- function(x, n, p){
  m <- colMeans(x)
  x - rep.int(m, rep_len(n, p))
}

fast.scale <- function(x, n, p){
  if(p > 1) {
    x <- fast.center(x, n, p)
    scale <- apply(x, 2, sd)
    x / rep.int(scale, rep_len(n, p))
  } else {
    x <- x - mean(x)
    x/sd(x)
  }
}

# csize
# calculates centroid size
# digitsurface
csize <- function(x) sqrt(sum(center(as.matrix(x))^2))

# cs.scale
# divide matrices by centroid size
# used in other functions for gpagen
cs.scale <- function(x) x/csize(x)

# center.scale
# center and divide matrices by centroid size; faster than scale()
# used in other functions for gpagen
center.scale <- function(x) {
  x <- center(x)
  cs <- sqrt(sum(x^2))
  y <- x/cs
  list(coords=y, CS=cs)
}

# apply.pPsup
# applies a partial Procrustes superimposition to matrices in a list
# used in gpagen functions
apply.pPsup<-function(M, Ya) {	# M = mean (reference); Ya all Y targets
  dims <- dim(Ya[[1]])
  k <- dims[2]; p <- dims[1]; n <- length(Ya)
  M <- cs.scale(M)
  lapply(1:n, function(j){
    y <- Ya[[j]]
    MY <- crossprod(M,y)
    sv <- La.svd(MY,k,k)
    u <- sv$u; u[,k] <- u[,k]*determinant(MY)$sign
    tcrossprod(y,u%*%sv$vt)
  })
}

# fast.ginv
# same as ginv, but without traps (faster)
# used in any function requiring a generalized inverse
fast.ginv <- function(X, tol = sqrt(.Machine$double.eps)){
  k <- ncol(X)
  if(inherits(X, "matrix")) {
    Xsvd <- La.svd(X, k, k)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    rtu <-((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
    v <-t(Xsvd$vt)[, Positive, drop = FALSE]
  } else {
    Xsvd <- svd(X, k, k)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    rtu <-((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
    v <-Xsvd$v[, Positive, drop = FALSE]
  }

  v %*% rtu
}

# fast.solve
# same as solve, but without traps (faster)
# used in any function requiring a generalized inverse
fast.solve <- function(x) { 
  res <- try(solve(x), silent = TRUE)
  if(inherits(res, "try-error")) res <- fast.ginv(x)
  return(res)
}


# pcoa
# acquires principal coordinates from distance matrices
# used in all linear model functions with data input
pcoa <- function(D){
  options(warn=-1)
  if(!inherits(D, "dist")) stop("function only works with distance matrices")
  cmd <- cmdscale(D, k=attr(D, "Size") -1, eig=TRUE)
  options(warn=0)
  d <- cmd$eig
  min.d <- min(d)
  if(min.d < 0) {
    options(warn=-1)
    cmd.c <- cmdscale(D, k=attr(D, "Size") -1, eig=TRUE, add= TRUE)
    options(warn=0)
    d <- cmd.c$eig
  } else cmd.c <- cmd
  p <- length(cmd.c$eig[zapsmall(d) > 0])
  Yp <- cmd.c$points[,1:p]
  Yp
}


# perm.index
# creates a permutation index for resampling
# used in all functions with a resampling procedure
perm.index <-function(n, iter, block = NULL, seed = NULL){
  
  if(!is.null(block)){
    block <- as.factor(block)
    if(length(block) != n)
      stop("block factor not the same length as the number of observations.\n",
           call. = FALSE)
    n <- length(block)
    indx <- 1:n
    g <- nlevels(block)
    gp <- levels(block)
    loc.list <- lapply(1:g, function(j) {
      indx[block == gp[[j]]]
    })
    rindx <- unlist(loc.list)
    
  }
  
  if(is.null(seed)) seed = iter else
    if(seed == "random") seed = sample(1:iter,1) else
      if(!is.numeric(seed)) seed = iter
      set.seed(seed)
      get.samp <- function(n){
        if(!is.null(block)){
          out <- array(NA, n)
          r <- unlist(lapply(loc.list, sample, replace = FALSE))
          for(i in indx) out[rindx[i]] <- r[i]
        } else  out <- sample.int(n, n)
        out
      }
      ind <- c(list(1:n), (Map(function(x) get.samp(n), 1:iter)))
      rm(.Random.seed, envir=globalenv())
      attr(ind, "seed") <- seed
      names(ind) <- c("obs", paste("iter", 1:(length(ind) - 1), sep = "."))
      ind
      
}


# boot.index
# creates a bootstrap index for resampling
# used in lm.rrpp for intercept models
boot.index <-function(n, iter, block = NULL, seed = NULL){
  
  if(!is.null(block)){
    block <- as.factor(block)
    if(length(block) != n)
      stop("block factor not the same length as the number of observations.\n",
           call. = FALSE)
    n <- length(block)
    indx <- 1:n
    g <- nlevels(block)
    gp <- levels(block)
    loc.list <- lapply(1:g, function(j) {
      indx[block == gp[[j]]]
    })
    rindx <- unlist(loc.list)
  }
  
  if(is.null(seed)) seed = iter else
    if(seed == "random") seed = sample(1:iter,1) else
      if(!is.numeric(seed)) seed = iter
      set.seed(seed)
      get.samp <- function(n){
        if(!is.null(block)){
          out <- array(NA, n)
          r <- unlist(lapply(loc.list, sample, replace = TRUE))
          for(i in indx) out[rindx[i]] <- r[i]
        } else  out <- sample.int(n, n)
        out
      }
      ind <- c(list(1:n), (Map(function(x) get.samp(n), 
                              1:iter)))
      rm(.Random.seed, envir=globalenv())
      attr(ind, "seed") <- seed
      ind
}

# fastFit
# calculates fitted values for a linear model, after decomoposition of X to get U
# used in SS.iter
fastFit <- function(U,y,n,p){
  if(p > n) tcrossprod(U) %*% y else
    U %*% crossprod(U,y)
} 


# fastLM
# calculates fitted values and residuals, after fastFit
# placeholder in case needed later
fastLM<- function(U, y){
  p <- dim(y)[2]; n <- dim(y)[1]
  yh <- fastFit(U,y,n,p)
  list(fitted = yh, residuals = y-yh)
}

# pval
# P-values form random outcomes
# any analytical function
pval = function(s){# s = sampling distribution
  p = length(s)
  r = rank(s)[1]-1
  pv = 1-r/p
  pv
}

# box.cox
# Box-Cox transformation for normalizing distributions.  Similar to MASS::boxcox
# without unneeded arguments, plus faster
# Used in effect.size
box.cox.true <- function(y, eps = 0.001){
  
  if(any(y <= 0)) y = y - min(y) + 0.0001
  
  y.obs <- y[1]
  y <- y[-1]
  
  # Note the code below (to get yy) is short-cut for Jacobian adjustment 
  n <- length(y)
  yy <- y / exp(mean(log(y)))
  logy <- log(yy)
  
  lambda <- seq(-5, 5, 0.001)
  m <- length(lambda)
  
  loglik <- sapply(1:m, function(j){ # same as MASS::boxcox loglik 
    la <- lambda[j]
    yt <- if(abs(la) > eps) yt <- (yy^la - 1)/la else
      logy * (1 + (la * logy)/2 * (1 + (la * logy)/3 * (1 + (la * logy)/4)))
    
    -n/2 * log(sum(center(yt)^2))
  })
  
  lambda.opt <- lambda[which.max(loglik)][[1]]
  if(abs(lambda.opt) < eps) lambda.opt <- 0
  y <- c(y.obs, y)
  res <- if(lambda.opt == 0) log(y) else (y^lambda.opt - 1)/lambda.opt
  
  list(opt.lambda = lambda.opt, transformed = res, lambda = lambda, 
       loglik = loglik)
  
}

box.cox.spline <- function(y, eps = 0.001) {
  
  if(any(y <= 0)) y = y - min(y) + 0.0001
  
  y.obs <- y[1]
  y <- y[-1]
  
  n <- length(y)
  yy <- y / exp(mean(log(y)))
  logy <- log(yy)
  m <- 20
  lambda <- seq(-0.5, 1.5, length.out = m)
  
  loglik <- sapply(1:m, function(j){ # same as MASS::boxcox loglik 
    la <- lambda[j]
    yt <- if(abs(la) > eps) yt <- (yy^la - 1)/la else
      logy * (1 + (la * logy)/2 * (1 + (la * logy)/3 * (1 + (la * logy)/4)))
    
    -n/2 * log(sum(center(yt)^2))
  })
  
  sp <- spline(lambda, loglik, n = 300)
  lambda.opt <- sp$x[which.max(sp$y)]
  if(abs(lambda.opt) < eps) lambda.opt <- 0
  y <- c(y.obs, y)
  res <- if(lambda.opt == 0) log(y) else (y^lambda.opt - 1)/lambda.opt
  
  list(opt.lambda = lambda.opt, transformed = res, lambda = sp$x, loglik = sp$y)
  
}

box.cox.iter <- function(y, eps = 0.001) {
  bc <- box.cox.spline(y, eps = eps)
  if(bc$opt.lambda == -0.5 || bc$opt.lambda == 1.5) 
    bc <- box.cox.true(y, eps = eps)
  return(bc)
}

box.cox.fast <- function(y, eps = 0.001) {
  if(any(y <= 0)) y = y - min(y) + 0.0001
  y.obs <- y[1]
  y <- y[-1]
  
  n <- length(y)
  yy <- y / exp(mean(log(y)))
  logy <- log(yy)
  
  logLik <- function(lambda) {
    la <- lambda
    yt <- if(abs(la) > eps) yt <- (yy^la - 1)/la else
      logy * (1 + (la * logy)/2 * (1 + (la * logy)/3 * (1 + (la * logy)/4)))
    
    -n/2 * log(sum(center(yt)^2))
  }
  
  result <- optimise(logLik, lower = -5, upper = 5, maximum = TRUE)
  lambda.opt <- result$maximum
  if(abs(lambda.opt) < eps) lambda.opt <- 0
  y <- c(y.obs, y)
  res <- if(lambda.opt == 0) log(y) else (y^lambda.opt - 1)/lambda.opt
  list(opt.lambda = lambda.opt, transformed = res, lambda = NULL, loglik = NULL)
}

box.cox <- function(y, eps = 0.001, iterate = FALSE) {
  result <- if(iterate) box.cox.iter(y, eps = eps) else
    box.cox.fast(y, eps = eps)
  return(result)
}


# effect.size
# Effect sizes (standard deviates) form random outcomes
# any analytical function

effect.size <- function(x, center = TRUE) {
  if(length(unique(x)) == 1) {
    sdx <- 1
    x <- 0
  } else {
    x <- box.cox(x)$transformed
    n <- length(x)
    if(center) x <- center(x)
    sdx <- sqrt((sum(x^2)/n))
  }
  
  (x[1]- mean(x)) / sdx
}

# Pval.matrix
# P-values form random outcomes that comprise matrices
# any analytical function with results in matrices
Pval.matrix = function(M){
  P = matrix(0,dim(M)[1],dim(M)[2])
  for(i in 1:dim(M)[1]){
    for(j in 1:dim(M)[2]){
      y = M[i,j,]
      p = pval(y)
      P[i,j]=p
    }
  }
  if(dim(M)[1] > 1 && dim(M)[2] >1) diag(P)=1
  rownames(P) = dimnames(M)[[1]]
  colnames(P) = dimnames(M)[[2]]
  P
}

# Effect.size.matrix
# Effect sizes form random outcomes that comprise matrices
# any analytical function with results in matrices
Effect.size.matrix <- function(M, center=TRUE){
  Z = matrix(0,dim(M)[1],dim(M)[2])
  for(i in 1:dim(M)[1]){
    for(j in 1:dim(M)[2]){
      y = M[i,j,]
      z = effect.size(y, center=center)
      Z[i,j]=z
    }
  }
  if(dim(M)[1] > 1 && dim(M)[2] >1) diag(Z)=0
  rownames(Z) = dimnames(M)[[1]]
  colnames(Z) = dimnames(M)[[2]]
  Z
}


# Cov.proj
# generates projection matrix from covariance matrix
# used in lm.rrpp

Cov.proj <- function(Cov, id = NULL, symmetric = FALSE){
  Cov <- if(is.null(id)) Cov else Cov[id, id]
  if(inherits(Cov, "matrix")) {
    Cov.s <- Matrix(Cov, sparse = TRUE)
    if(object.size(Cov.s) < object.size(Cov)) Cov <- Cov.s
    Cov.s <- NULL
  }
  ow <- options()$warn
  options(warn = -1)
  if(!symmetric) {
    Chol <- try(chol(Cov), silent = TRUE)
    if(inherits(Chol, "try-error")) 
      symmetric <- TRUE
  }
  
 if(symmetric) {
   sym <- isSymmetric(Cov)
   eigC <- eigen(Cov, symmetric = sym)
   eigC.vect = t(eigC$vectors)
   L <- eigC.vect * sqrt(abs(eigC$values))
   P <- fast.solve(crossprod(L, eigC.vect))
   dimnames(P) <- dimnames(Cov)
 } else P <- solve(t(Chol))
  options(warn = ow)
  P
}



# ape replacement functions below --------------------------------------------

# sim.char has a similar function to this but it is called in every simulation
# and defers to C for help.  This is done once only here. 
# (Produces a projection matrix)
phy.sim.mat <- function(phy) {
  N <- length(phy$tip.label)
  n <- nrow(phy$edge)
  m <- matrix(0, N, n)
  edg <- phy$edge.length
  idx <- phy$edge[, 2]
  anc <- phy$edge[, 1]
  tips <- which(idx <= N)
  non.tips <- which(idx > N)
  for(i in 1:length(tips)) m[idx[tips[i]], tips[i]] <- sqrt(edg[tips[i]])
  for(i in 1:length(non.tips)) {
    x <- idx[non.tips[i]]
    anc.i <- which(anc == x)
    edg.i <- idx[anc.i]
    if(any(edg.i > N)) {
      edg.i <- as.list(edg.i)
      while(any(edg.i > N)) {
        edg.i <- lapply(1:length(edg.i), function(j){
          edg.i.j <- edg.i[[j]]
          if(edg.i.j > N) {
            idx[which(anc == edg.i.j)]
          } else edg.i.j <- edg.i.j
        })
        edg.i <- unlist(edg.i)
      }
    }
    m[edg.i, which(idx == x)] <- sqrt(edg[which(idx == x)])
  }
  rownames(m) <- phy$tip.label
  m
}

# fast.phy.vcv
# same as vcv.phylo but without options, in order to not use ape

fast.phy.vcv <- function (phy) tcrossprod(phy.sim.mat(phy))

# reorder.phy
# same as reorder function, but without options

reorder.phy <- function(phy){
  edge <- phy$edge
  edge.length <- phy$edge.length
  edge <- cbind(edge, edge.length)
  n <- nrow(edge)
  ind <-rank(edge[,1], ties.method = "last")
  edge <- edge[order(ind, decreasing = TRUE), ]
  edge.length <- edge[,3]
  phy$edge <- edge[,-3]
  phy$edge.length <- edge.length
  attr(phy, "order") <- "postorder"
  return(phy)
}


# anc.BM 
# via PICs
# same as ace, but multivariate
pic.prep <- function(phy, nx, px){
  phy <- reorder.phy(phy)
  ntip <- length(phy$tip.label)
  nnode <- phy$Nnode
  edge <- phy$edge
  edge1 <- edge[, 1]
  edge2 <- edge[,2]
  edge_len <- phy$edge.length
  phe <- matrix(0, ntip + nnode, px)
  contr <- matrix(0, nnode, px)
  var_contr <- rep(0, nnode)
  i.seq <- seq(1, ntip * 2 -2, 2)
  list(ntip = ntip, nnode = nnode, edge1 = edge1,
       edge2 = edge2, edge_len = edge_len, phe = phe,
       contr = contr, var_contr = var_contr, 
       tip.label = phy$tip.label,
       i.seq = i.seq)
}

ace.pics <- function(ntip, nnode, edge1, edge2, edge_len, phe, contr,
                     var_contr, tip.label, i.seq, x) {
  phe[1:ntip,] <- if (is.null(rownames(x))) x else x[tip.label,]
  N <- ntip + nnode
  for(ii in 1:nnode) {
    anc <- edge1[i.seq[ii]]
    ij <- which(edge1 == anc)
    i <- ij[1]
    j <- ij[2]
    d1 <- edge2[i]
    d2 <- edge2[j]
    sumbl <- edge_len[i] + edge_len[j]
    ic <- anc - ntip 
    ya <- (phe[d1,] - phe[d2,])/sqrt(sumbl)
    contr[ic, ] <- ya
    var_contr[ic] <- sumbl
    phe[anc,] <- (phe[d1, ] * edge_len[j] + phe[d2, ] * edge_len[i])/sumbl
    k <- which(edge2 == anc)
    edge_len[k] <- edge_len[k] + edge_len[i] * edge_len[j] / sumbl
  }
  phe
}

# anc.BM
# multivariate as opposed to fastAnc

anc.BM <- function(phy, Y){
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  Y <- as.matrix(Y[phy$tip.label,])
  phy <- reorder.phy(phy)
  n <- length(phy$tip.label)
  out <- t(sapply(1:phy$Nnode, function(j){
    phy.j <- multi2di.phylo(root.phylo(phy, node = j + n))
    phy.j <- collapse.singles(phy.j)
    preps <- pic.prep(phy.j, NROW(Y), NCOL(Y))
    preps$x <- Y
    preps$tip.label <- phy$tip.label
    out <- do.call(ace.pics, preps)
    out[n + 1,]
  }))
  
  if(NROW(out) == 1) out <- t(out)
  dimnames(out) <- list(1:phy$Nnode + length(phy$tip.label), colnames(Y))
  out
}

# getNode Depth
# replaces node.depth.edgelength

getNodeDepth <- function(phy){
  phy <- reorder.phy(phy)
  E <- phy$edge
  anc <- E[,1]
  des <- E[,2]
  ntip <- length(phy$tip.label)
  nnode <- phy$Nnode
  N <- ntip + nnode
  L <- phy$edge.length
  base.tax <- ntip + 1
  full.depth.seq <- (ntip + 1):N
  full.node.depth <- 0
  tips <- which(des <= ntip)
  non.tips <- which(des > ntip)
  
  get.edge.ind <- function(tax){
    root.t <- ntip +1
    des.i <- which(des == tax)
    edge <- numeric()
    tax.i <- tax
    while(tax.i != root.t) {
      anc.i <- anc[which(des == tax.i)]
      tax.i <- anc[des.i]
      edge <- c(edge, des.i)
      des.i <- which(des == anc.i)
    }
    edge
  }
  
  tips.taxa <- lapply(as.list(tips), function(j) des[j])
  tips.edges <- lapply(tips.taxa, get.edge.ind)
  tips.depths <- sapply(1:ntip, function(j) sum(L[tips.edges[[j]]]))
  
  nodes <- as.list((ntip + 1):(N))
  nodes.edges <- lapply(nodes, get.edge.ind)
  nodes.depths <- sapply(1:nnode, function(j) sum(L[nodes.edges[[j]]]))
  
  c(tips.depths, nodes.depths)
}
