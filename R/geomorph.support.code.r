#' @name geomorph-package
#' @docType package
#' @aliases geomorph
#' @title Geometric morphometric analyses for 2D/3D data
#' @author Dean C. Adams, Michael Collyer, & Emma Sherratt
#'
#' Functions in this package allow one to read, manipulate, and digitize landmark data; generate shape
#'  variables via Procrustes analysis for points, curves and surface data, perform statistical analyses
#'  of shape variation and covariation, and provide graphical depictions of shapes and patterns of
#'  shape variation.
#' 
#' 
#' @import ape
#' @import rgl
#' @import stats
#' @import utils
#' @import graphics
#' @import grDevices
#' @importFrom geiger sim.char
#' @importFrom jpeg readJPEG
#' @importFrom phytools fastAnc
#' @importFrom Matrix nearPD
#' @importFrom MASS cmdscale

NULL

#' Landmark data from Plethodon salamander heads
#'
#' @name plethodon
#' @docType data
#' @author Dean Adams
#' @references Adams, D. C. 2004. Character displacement via aggressive interference in Appalachian salamanders. 
#' Ecology. 85:2664-2670.
#' @references Adams, D.C. 2010. Parallel evolution of character displacement driven by competitive selection 
#' in terrestrial salamanders. BMC Evolutionary Biology. 10(72)1-10.
#' @keywords datasets
NULL

#' Head shape and food use data from Plethodon salamanders
#'
#' @name plethShapeFood
#' @docType data
#' @author Dean Adams
#' @references Adams, D. C., and F. J. Rohlf. 2000. Ecological character 
#' displacement in Plethodon: biomechanical differences found from a geometric 
#' morphometric study. Proceedings of the National Academy of Sciences, 
#' U.S.A. 97:4106-4111
#' @keywords datasets
NULL

#' Landmark data from dataset rat
#'
#' @name ratland
#' @docType data
#' @author Dean Adams
#' @references Bookstein, F. L. 1991. Morphometric tools for landmark data: Geometry and Biology. 
#'  Cambridge Univ. Press, New York.
#' @keywords datasets
NULL

#' Landmark data from hummingbird bills (includes sliding semilandmarks on curves)
#'
#' @name hummingbirds
#' @docType data
#' @author Chelsea Berns and Dean Adams
#' @references Berns, C.M., and Adams, D.C. 2010. Bill shape and sexual shape dimorphism between two species 
#' of temperate hummingbirds: Archilochus alexandri (black-chinned hummingbirds) and Archilochus colubris 
#' (ruby-throated hummingbirds). The Auk. 127:626-635.
#' @keywords datasets
NULL

#' Average head shape and phylogenetic relationships for several Plethodon salamander species
#'
#' @name plethspecies
#' @docType data
#' @author Dean Adams
#' @references Phylogeny pruned from: Wiens et al. (2006). Evol.
#' @references Data from: Adams and Rohlf (2000); Adams et al. (2007); Arif et al. (2007) Myers and Adams (2008)
#' @keywords datasets
NULL

#' Landmark data from scallop shells
#'
#' @name scallops
#' @docType data
#' @author Dean Adams and Erik Otarola-Castillo
#' @references Serb et al. (2011). "Morphological convergence of shell shape in distantly related
#' scallop species (Mollusca: Pectinidae)." Zoological Journal of the Linnean Society 163: 571-584.
#' @keywords datasets
NULL

#' 3D scan of a scallop shell from a .ply file in mesh3d format
#'
#' @name scallopPLY
#' @docType data
#' @author Emma Sherratt
#' @references Serb et al. (2011). "Morphological convergence of shell shape in distantly related
#' scallop species (Mollusca: Pectinidae)." Zoological Journal of the Linnean Society 163: 571-584.
#' @keywords datasets
NULL

#' Simulated motion paths
#'
#' @name motionpaths
#' @docType data
#' @author Dean Adams
#' @references Adams, D. C., and M. L. Collyer. 2009. A general framework for the analysis of phenotypic 
#'   trajectories in evolutionary studies. Evolution 63:1143-1154.
#' @keywords datasets
NULL

#' landmarks on mosquito wings
#'
#' @name mosquito
#' @docType data
#' @author Dean Adams
#' @keywords datasets
NULL

#' landmarks on pupfish
#'
#' @name pupfish
#' @docType data
#' @author Michael Collyer
#' @keywords datasets
#' @description Landmark data from Cyprindon pecosensis body shapes, with indication of Sex and 
#' Population from which fish were sampled (Marsh or Sinkhole).  These data were previously aligned 
#' with GPA.  Centroid size (CS) is also provided.
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic 
#' change for phenotypes described by high-dimensional data. Heredity. 113: doi:10.1038/hdy.2014.75.
NULL

#####----------------------------------------------------------------------------------------------------

# SUPPORT FUNCTIONS

# center
# centers a matrix faster than scale()
# used in other functions for gpagen
center <- function(x) x - rep(colMeans(x), rep.int(nrow(x), ncol(x)))

# cs.scale
# divide matrices by centroid size
# used in other functions for gpagen
cs.scale <- function(x) x/sqrt(sum(center(x)^2))

# center.scale
# center and divide matrices by centroid size; faster than scale()
# used in other functions for gpagen
center.scale <- function(x) {
  x <- center(x)
  cs <- sqrt(sum(x^2))
  y <- x/cs
  list(coords=y, CS=cs)
}

# orp
# projection in GPA
# used in gpagen functons
orp<-function(A){			
  if(is.array(A)) {
    n<-dim(A)[3]; k<-dim(A)[2]; p<-dim(A)[1]  
    Y <- lapply(1:n, function(j) A[,,j])
  } else
    if(is.list(A)){
      Y <- A
      n <- length(A); k <- ncol(A[[1]]); p <- nrow(A[[1]])
    } else stop("Input must be either a list or array")
  
  Y1<-as.vector(center.scale((Reduce("+", Y)/n))$coords)
  oo<-as.matrix(rep(1,n))%*%Y1
  mat <- t(matrix(unlist(Y),k*p,n))
  Xp <- (mat%*%(diag(1,p*k) - (tcrossprod(Y1)))) +oo
  lapply(1:n, function(j) matrix(Xp[j,],p,k))
}

# apply.pPsup
# applies a partial Procrustes superimposition to matrices in a list
# used in gpagen functions
apply.pPsup<-function(M, Ya) {	# MY = list of cross products of M and Yi; Y = List of Yi
  k <- ncol(Ya[[1]]); p <- nrow(Ya[[1]]); n <- length(Ya)
  M <- cs.scale(M)
  MY <- Map(function(y) .Internal(crossprod(M,y)), Ya)
  sv <- Map(function(y) .Internal(La_svd("A", y, double(k), matrix(0,k,k), matrix(0,k,k))), MY)
  sig <- Map(function(y) determinant(y)$sign, MY)
  Ut <- lapply(1:n, function(j) {u <- sv[[j]]$u; u[,k] <- u[,k]*sig[[j]]; t(u)})
  V <- Map(function(y) y$vt, sv)
  Yrot <- Map(function(u,v,y) y%*%.Internal(crossprod(v,u)), Ut,V,Ya)
  Yrot
}

# fast.ginv
# same as ginv, but without traps (faster)
# used in any function requiring a generalized inverse
fast.ginv <- function(X, tol = sqrt(.Machine$double.eps)){
  k <- ncol(X)
  Xsvd <- .Internal(La_svd("A", X, double(k), matrix(0,k,k), matrix(0,k,k)))
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  rtu <-((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
  v <-t(Xsvd$vt)[, Positive, drop = FALSE]
  v%*%rtu
}

# fast.solve
# chooses between fast.ginv or qr.solve, when det might or might not be 0
# used in any function requring a matrix inverse where the certainty of
# singular matrices is in doubt; mostly phylo. functions
fast.solve <- function(x) if(det(x) > 1e-8) qr.solve(x) else fast.ginv(x)

# tangents
# finds tangents in a matrix based on sliders
# used in all functions associated with pPga.wCurvs
tangents = function(s,x, scaled=FALSE){ # s = curves, x = landmarks
  ts <- x[s[,3],] - x[s[,1],]
  if(scaled==TRUE) {
    ts.scale = sqrt(rowSums(ts^2))
    ts <- ts/ts.scale
  }
  y <- matrix(0, nrow(x), ncol(x))
  y[s[,2],] <- ts
  y
}

# nearest
# finds nearest points on surfaces for sliding semilandmakrs
# used in all functions associated with pPga.wCurves
nearest <- function(X,m,k=4) {
  a <- X[m,]
  b <- sapply(1:nrow(X), function (j) sum((a-X[j,])^2))
  match(sort(b)[2:(k+1)],b)
}

# getU
# calculates U matrix for sliding semilandmarks
# currently not used but retained for posterity
getU <- function(y,tn, surf){
  p <- nrow(tn); k <- ncol(tn)
  Ux <- Uy <- Uz <- matrix(0,p,p)
  if(!is.null(tn)) {
    if(k == 3){ diag(Ux) <- tn[,1]; diag(Uy) <- tn[,2]; diag(Uz) <- tn[,3]} else
    { diag(Ux) <- tn[,1]; diag(Uy) <- tn[,2]; Uz <- NULL}
    U <- rbind(Ux,Uy,Uz)
  }
  if(!is.null(surf)){
    U2 <- array(0,dim=c(k*p,p))  
    U <- cbind(U,U2)
    nearpts <- sapply(1:length(surf), function(j) nearest(y,surf[j], k=4))
    nearpts <- cbind(t(nearpts), surf)
    tmp.pts <- lapply(1:nrow(nearpts), function(j) {
      k <- nearpts[j,]
      x <- y[k,]; x})
    pc.dir<-Map(function(y) .Internal(La_svd("A", var(y), double(k), matrix(0,k,k), matrix(0,k,k)))$u, tmp.pts)
    z11 <- cbind(surf,surf); z21 <- cbind(p+surf, surf); if(k==3) z31 <- cbind(2*p+surf, surf)
    z12 <- cbind(surf,p+surf); z22 <- cbind(p+surf, p+surf); if(k==3) z32 <- cbind(2*p+surf, p+surf)
    pc11 <- sapply(1:length(surf), function(j) pc.dir[[j]][1,1])
    pc12 <- sapply(1:length(surf), function(j) pc.dir[[j]][1,2])
    if(k==3) pc13 <- sapply(1:length(surf), function(j) pc.dir[[j]][1,3])
    pc21 <- sapply(1:length(surf), function(j) pc.dir[[j]][2,1])
    pc22 <- sapply(1:length(surf), function(j) pc.dir[[j]][2,2])
    if(k==3) pc23 <- sapply(1:length(surf), function(j) pc.dir[[j]][2,3])
    U[z11] <- pc11; U[z21] <- pc21; U[z21] <- pc12; U[z22] <- pc22
    if(k==3) U[z31] <- pc13; U[z32] <- pc23
    }                   
    U                  
}

# Ltemplate
# calculates inverse of bending energy matrix
# used in any function that calculates bending energy
# currently not used but retained for posterity
Ltemplate <-function(Mr, Mt=NULL){
  p <-nrow(Mr); k <- ncol(Mr)
  if(!is.null(Mt)) P <- as.matrix(dist(Mr-Mt)) else P <- as.matrix(dist(Mr))
  if(k==2) {P <-P^2*log(P); P[is.na(P)] <- 0}
  Q <- rbind(cbind(1,Mr))
  L<-rbind(cbind(P,Q), cbind(t(Q),matrix(0,k+1,k+1)))
  Linv <- -fast.ginv(L)[1:p,1:p]
  Linv
}

# Ltemplate
# calculates inverse of bending energy matrix and expand it to dimensions of landmarks
# only used if getU is used
# used in a function that calculates bending energy, if get U is used
bigLtemplate <-function(Mr, Mt=NULL){
  p <-nrow(Mr); k <- ncol(Mr)
  if(!is.null(Mt)) P <- as.matrix(dist(Mr-Mt)) else P <- as.matrix(dist(Mr))
  if(k==2) {P <-P^2*log(P); P[is.na(P)] <- 0}
  Q <- rbind(cbind(1,Mr))
  L<-rbind(cbind(P,Q), cbind(t(Q),matrix(0,k+1,k+1)))
  Linv <- -fast.ginv(L)[1:p,1:p]
  if(k==2) Linv <- rbind(cbind(Linv,array(0,dim(Linv))),cbind(array(0,dim(Linv)),Linv))
  if(k==3) Linv <- rbind(cbind(Linv,array(0,dim(Linv)), array(0,dim(Linv))),
                         cbind(array(0,dim(Linv)),Linv,array(0,dim(Linv))),
                         cbind(array(0,dim(Linv)),array(0,dim(Linv)),Linv))
  Linv
}

# pGPA
# GPA with partial Procrustes superimposition
# used in gpagen
pGpa <- function(Y, PrinAxes = FALSE, Proj = FALSE, max.iter = 5){
  n <- length(Y); p <- nrow(Y[[1]]); k <- ncol(Y[[1]])
  Yc <- Map(function(y) center.scale(y), Y)
  CS <- sapply(Yc,"[[","CS")
  Ya <- lapply(Yc,"[[","coords")
  M <- Reduce("+",Ya)/n
  Ya <- apply.pPsup(M, Ya)
  M <- Reduce("+",Ya)/n
  Q <- ss <- n*(1-sum(M^2))
  M <- cs.scale(M)
  iter <- 1
  while(Q > 0.0001){
    iter <- iter+1
    Ya <- apply.pPsup(M, Ya)
    M <- Reduce("+",Ya)/n
    ss2 <- n*(1-sum(M^2))
    Q <- abs(ss-ss2)
    ss <- ss2
    M <- cs.scale(M)
    if(iter > max.iter) break
  }
  if (PrinAxes == TRUE) {
    ref <- M
    rot <- prcomp(ref)$rotation
    for (i in 1:k) if (sign(rot[i, i]) != 1) 
      rot[1:k, i] = -rot[1:k, i]
    Ya <- Map(function(y) y%*%rot, Ya)
    M <- cs.scale(Reduce("+", Ya)/n)
  }
  list(coords= Ya, CS=CS, iter=iter, consensus=M, Q=Q, nsliders=NULL)
}

# getSurfPCs
# finds PC loadings for surface landmarks
# used in semilandmarks functions, within the larger gpagen framework
getSurfPCs <- function(y, surf){
  p <- nrow(y); k <- ncol(y)
  pc.match <- as.vector(match(1:p, surf))
  nearpts <- lapply(1:length(pc.match), function(j) {
    x <- pc.match[j]
    if(is.na(x)) 0 else 
      c(nearest(y,pc.match[j], k=k+1), x)})
  tmp.pts <- lapply(1:p, function(j) {
    k <- nearpts[[j]]
    if(sum(k) > 0) x <- y[k,] else x <- NA
    x})
  pc.dir <- lapply(1:p, function(j) {
    y <- tmp.pts[[j]]
    if(is.matrix(y)) 
      .Internal(La_svd("A", var(y), double(k), matrix(0,k,k), matrix(0,k,k)))$u else 0
  })
  p1x <- sapply(1:p, function(j) {x <- pc.dir[[j]]; if(is.matrix(x)) x[1,1] else 0})
  p1y <- sapply(1:p, function(j) {x <- pc.dir[[j]]; if(is.matrix(x)) x[1,2] else 0})
  p2x <- sapply(1:p, function(j) {x <- pc.dir[[j]]; if(is.matrix(x)) x[2,1] else 0})
  p2y <- sapply(1:p, function(j) {x <- pc.dir[[j]]; if(is.matrix(x)) x[2,2] else 0})
  if(k==3) {
    p1z <- sapply(1:p, function(j) {x <- pc.dir[[j]]; if(is.matrix(x)) x[1,3] else 0})
    p2z <- sapply(1:p, function(j) {x <- pc.dir[[j]]; if(is.matrix(x)) x[2,3] else 0})
  } else
  {p1z <- NULL; p2z <- NULL}
    
  list(p1x=p1x,p1y=p1y, p2x=p2x, p2y=p2y, p1z=p1z, p2z=p2z)
}

# semilandmarks.slide.tangents.BE
# slides landmarks along tangents of curves using bending energy
# used in pGpa.wSliders
semilandmarks.slide.tangents.BE <- function(y, tans, ref, L){
  yc <- y - ref
  p <- nrow(yc); k <-ncol(yc)
  if(k==3) {tx <- tans[,1]; ty <- tans[,2]; tz <- tans[,3 ]} else {tx <- tans[,1]; ty <- tans[,2]}
 
    if(k==3) {
      int.part <- fast.ginv(t(t(tx*L)*tx)+t(t(ty*L)*ty)+ 
                              t(t(tz*L)*tz))%*%cbind(tx*L,ty*L,tz*L)
      Ht <- rbind(tx*int.part, ty*int.part, tz*int.part) 
    } else {
      int.part <- fast.ginv(t(t(tx*L)*tx)+t(t(ty*L)*ty))%*%cbind(tx*L,ty*L)
      Ht <- rbind(tx*int.part, ty*int.part) 
    }
  y  - matrix(Ht%*%as.vector(yc), p,k)
}

# semilandmarks.slide.surf.BE
# slides landmarks in PC planes tangent to surfaces using bending energy
# used in pGpa.wSliders
semilandmarks.slide.surf.BE <- function(y, surf, ref, L){
  yc <- y - ref
  p <- nrow(yc); k <-ncol(yc)
    PC <- getSurfPCs(y, surf)
    p1x <- PC$p1x; p1y <- PC$p1y; p1z <- PC$p1z; p2x <- PC$p2x; p2y <- PC$p2y; p2z <- PC$p2z
    if(k==3) {
      int.part <- fast.ginv(t(t(p1x*L)*p1x)+t(t(p1y*L)*p1y)+ 
                              t(t(p1z*L)*p1z))%*%cbind(p1x*L,p1y*L,p1z*L)
      Hp1 <- rbind(p1x*int.part, p1y*int.part, p1z*int.part) 
    } else {
      int.part <- fast.ginv(t(t(p1x*L)*p1x)+t(t(p1y*L)*p1y))%*%cbind(p1x*L,p1y*L)
      Hp1 <- rbind(p1x*int.part, p1y*int.part) 
    }
    if(k==3) {
      int.part <- fast.ginv(t(t(p2x*L)*p2x)+t(t(p2y*L)*p2y)+ 
                              t(t(p2z*L)*p2z))%*%cbind(p2x*L,p2y*L,p2z*L)
      Hp2 <- rbind(p2x*int.part, p2y*int.part, p2z*int.part) 
    } else {
      int.part <- fast.ginv(t(t(p2x*L)*p2x)+t(t(p2y*L)*p2y))%*%cbind(p2x*L,p2y*L)
      Hp2 <- rbind(p2x*int.part, p2y*int.part) 
    }
    y  - matrix( Hp1%*%as.vector(yc) + Hp2%*%as.vector(yc), p,k) 
}

# semilandmarks.slide.tangents.surf.BE
# slides landmarks along tangents of curves and PC planes of surfaces using bending energy
# used in pGpa.wSliders
semilandmarks.slide.tangents.surf.BE <- function(y, tans, surf, ref, L){
  yc <- y - ref
  p <- nrow(yc); k <-ncol(yc)
  if(k==3) {tx <- tans[,1]; ty <- tans[,2]; tz <- tans[,3 ]} else {tx <- tans[,1]; ty <- tans[,2]}
    if(k==3) {
      int.part <- fast.ginv(t(t(tx*L)*tx)+t(t(ty*L)*ty)+ 
                              t(t(tz*L)*tz))%*%cbind(tx*L,ty*L,tz*L)
      Ht <- rbind(tx*int.part, ty*int.part, tz*int.part) 
    } else {
      int.part <- fast.ginv(t(t(tx*L)*tx)+t(t(ty*L)*ty))%*%cbind(tx*L,ty*L)
      Ht <- rbind(tx*int.part, ty*int.part) 
    }
    PC <- getSurfPCs(y, surf)
    p1x <- PC$p1x; p1y <- PC$p1y; p1z <- PC$p1z; p2x <- PC$p2x; p2y <- PC$p2y; p2z <- PC$p2z
    if(k==3) {
      int.part <- fast.ginv(t(t(p1x*L)*p1x)+t(t(p1y*L)*p1y)+ 
                              t(t(p1z*L)*p1z))%*%cbind(p1x*L,p1y*L,p1z*L)
      Hp1 <- rbind(p1x*int.part, p1y*int.part, p1z*int.part) 
    } else {
      int.part <- fast.ginv(t(t(p1x*L)*p1x)+t(t(p1y*L)*p1y))%*%cbind(p1x*L,p1y*L)
      Hp1 <- rbind(p1x*int.part, p1y*int.part) 
    }
    if(k==3) {
      int.part <- fast.ginv(t(t(p2x*L)*p2x)+t(t(p2y*L)*p2y)+ 
                              t(t(p2z*L)*p2z))%*%cbind(p2x*L,p2y*L,p2z*L)
      Hp2 <- rbind(p2x*int.part, p2y*int.part, p2z*int.part) 
    } else {
      int.part <- fast.ginv(t(t(p2x*L)*p2x)+t(t(p2y*L)*p2y))%*%cbind(p2x*L,p2y*L)
      Hp2 <- rbind(p2x*int.part, p2y*int.part) 
    }
  y  - matrix(Ht%*%as.vector(yc) + Hp1%*%as.vector(yc) + Hp2%*%as.vector(yc), p,k) 
}

# semilandmarks.slide.tangents.procD
# slides landmarks along tangents of curves using minimized ProcD
# used in pGpa.wSliders
semilandmarks.slide.tangents.procD <- function(y,tans, ref){
  yc <- y - ref
  p <- nrow(yc); k <-ncol(yc)
  if(k==3) {ycx <- yc[,1]; ycy <- yc[,2]; ycz <- yc[,3 ]} else {ycx <- yc[,1]; ycy <- yc[,2]}
  if(k==3) {tx <- tans[,1]; ty <- tans[,2]; tz <- tans[,3 ]} else {tx <- tans[,1]; ty <- tans[,2]}
  if(k==3){
    sx <-tx*tx*ycx+tx*ty*ycy+tx*tz*ycz
    sy <-ty*tx*ycx+ty*ty*ycy+ty*tz*ycz
    sz <-tz*tx*ycx+tz*ty*ycy+tz*tz*ycz
  } else
  {
    sx <-tx*tx*ycx+tx*ty*ycy
    sy <-ty*tx*ycx+ty*ty*ycy
    sz <- NULL
  }
  y - cbind(sx,sy,sz)
}

# semilandmarks.slide.surf.procD
# slides landmarks within PC planes tangent to surfaces using minimized ProcD
# used in pGpa.wSliders
semilandmarks.slide.surf.procD <- function(y,surf, ref){
  yc <- y - ref
  p <- nrow(yc); k <-ncol(yc)
  PC <- getSurfPCs(y, surf)
  p1x <- PC$p1x; p1y <- PC$p1y; p1z <- PC$p1z; p2x <- PC$p2x; p2y <- PC$p2y; p2z <- PC$p2z
  if(k==3) {ycx <- yc[,1]; ycy <- yc[,2]; ycz <- yc[,3 ]} else {ycx <- yc[,1]; ycy <- yc[,2]}
  if(k==3){
    sx1 <-p1x*p1x*ycx+p1x*p1y*ycy+p1x*p1z*ycz
    sy1 <-p1y*p1x*ycx+p1y*p1y*ycy+p1y*p1z*ycz
    sz1 <-p1z*p1x*ycx+p1z*p1y*ycy+p1z*p1z*ycz
  } else
  {
    sx1 <-p1x*p1x*ycx+p1x*p1y*ycy
    sy1 <-p1y*p1x*ycx+p1y*p1y*ycy
    sz1 <- NULL
  }
  if(k==3){
    sx2 <-p2x*p2x*ycx+p2x*p2y*ycy+p2x*p2z*ycz
    sy2 <-p2y*p2x*ycx+p2y*p2y*ycy+p2y*p2z*ycz
    sz2 <-p2z*p2x*ycx+p2z*p2y*ycy+p2z*p2z*ycz
  } else
  {
    sx2 <-p2x*p2x*ycx+p2x*p2y*ycy
    sy2 <-p2y*p2x*ycx+p2y*p2y*ycy
    sz2 <- NULL
  }
  y - (cbind(sx1,sy1,sz1)+cbind(sx2,sy2,sz2))
}

# semilandmarks.slide.tangents.surf.procD
# slides landmarks along tangents of curves and within tangent planes on surfaces using minimized ProcD
# used in pGpa.wSliders
semilandmarks.slide.tangents.surf.procD <- function(y,tans,surf, ref){
  yc <- y - ref
  p <- nrow(yc); k <-ncol(yc)
  if(k==3) {tx <- tans[,1]; ty <- tans[,2]; tz <- tans[,3 ]} else {tx <- tans[,1]; ty <- tans[,2]}
  PC <- getSurfPCs(y, surf)
  p1x <- PC$p1x; p1y <- PC$p1y; p1z <- PC$p1z; p2x <- PC$p2x; p2y <- PC$p2y; p2z <- PC$p2z
  if(k==3) {ycx <- yc[,1]; ycy <- yc[,2]; ycz <- yc[,3 ]} else {ycx <- yc[,1]; ycy <- yc[,2]}
  if(k==3){
    sxt <-tx*tx*ycx+tx*ty*ycy+tx*tz*ycz
    syt <-ty*tx*ycx+ty*ty*ycy+ty*tz*ycz
    szt <-tz*tx*ycx+tz*ty*ycy+tz*tz*ycz
  } else
  {
    sxt <-tx*tx*ycx+tx*ty*ycy
    syt <-ty*tx*ycx+ty*ty*ycy
    szt <- NULL
  }
  if(k==3){
    sx1 <-p1x*p1x*ycx+p1x*p1y*ycy+p1x*p1z*ycz
    sy1 <-p1y*p1x*ycx+p1y*p1y*ycy+p1y*p1z*ycz
    sz1 <-p1z*p1x*ycx+p1z*p1y*ycy+p1z*p1z*ycz
  } else
  {
    sx1 <-p1x*p1x*ycx+p1x*p1y*ycy
    sy1 <-p1y*p1x*ycx+p1y*p1y*ycy
    sz1 <- NULL
  }
  if(k==3){
    sx2 <-p2x*p2x*ycx+p2x*p2y*ycy+p2x*p2z*ycz
    sy2 <-p2y*p2x*ycx+p2y*p2y*ycy+p2y*p2z*ycz
    sz2 <-p2z*p2x*ycx+p2z*p2y*ycy+p2z*p2z*ycz
  } else
  {
    sx2 <-p2x*p2x*ycx+p2x*p2y*ycy
    sy2 <-p2y*p2x*ycx+p2y*p2y*ycy
    sz2 <- NULL
  }
  y - (cbind(sxt,syt,szt)+cbind(sx1,sy1,sz1)+cbind(sx2,sy2,sz2))
}

# BE.slide
# performs sliding iterations using bending energy
# used in pGpa.wSliders
BE.slide <- function(curves, surf, Ya, ref, max.iter=5){# see pGpa.wCurves for variable meaning
  n <- length(Ya); p <- nrow(Ya[[1]]); k <- ncol(Ya[[1]])
  iter <- 1 # from initial rotation of Ya
  slid0 <- Ya
  Q <-ss0 <-1
  while(Q > 0.0001){
    iter <- iter+1
    if(!is.null(curves)) tans <- Map(function(y) tangents(curves, y, scaled=TRUE), slid0)
    L <- Ltemplate(ref)
    if(is.null(surf) & !is.null(curves))
      slid <- Map(function(tn,y) semilandmarks.slide.tangents.BE(y, tn, ref, L), tans, slid0)
    if(!is.null(surf) & is.null(curves))
      slid <- Map(function(y) semilandmarks.slide.surf.BE(y, surf, ref, L), slid0)
    if(!is.null(surf) & !is.null(curves))
      slid <- Map(function(tn,y) semilandmarks.slide.tangents.surf.BE(y, tn, surf, ref, L), tans, slid0)
    M <- Reduce("+", slid)/n
    ss <-(1-sum(M^2))*n
    ref <- cs.scale(M)
    Ya <- apply.pPsup(ref, slid)
    Q <- abs(ss0-ss)
    slid0 <- Ya
    ss0 <- ss
    if(iter >= max.iter) break
  }
  gpa.final <- pGpa(Ya, PrinAxes = F, Proj=F)
  list(coords=gpa.final$coords, consensus=gpa.final$consensus, iter=iter+1, Q=Q)
}

# procD.slide
# performs sliding iterations using minimized ProcD
# used in pGpa.wSliders
procD.slide <- function(curves, surf, Ya, ref, max.iter=5){# see pGpa.wCurves for variable meaning
  n <- length(Ya); p <- nrow(Ya[[1]]); k <- ncol(Ya[[1]])
  iter <- 1 # from initial rotation of Ya
  slid0 <- Ya
  iter.s <- 0
  Q <- QQ <-1
  ss0 <- n*(1-sum(ref^2))
  while(Q > 0.0001){
    iter <- iter+1
    while(QQ > 0.0001){
      iter.s <- iter.s + 1
      if(!is.null(curves)) tans <- Map(function(y) tangents(curves, y, scaled=TRUE), slid0)
      if(is.null(surf) & !is.null(curves))
        slid <- Map(function(tn,y) semilandmarks.slide.tangents.procD(y, tn, ref), tans, slid0)
      if(!is.null(surf) & is.null(curves))
        slid <- Map(function(y) semilandmarks.slide.surf.procD(y, surf, ref), slid0)
      if(!is.null(surf) & !is.null(curves))
        slid <- Map(function(tn,y) semilandmarks.slide.tangents.surf.procD(y, tn, surf, ref), tans, slid0)
      M <- Reduce("+", slid)/n
      ss2 <-(1-sum(M^2))*n
      slid0 <- slid
      ss <- ss2
      ref <- M
      QQ <- abs(ss-ss2)
      if(iter.s >= max.iter) break
    }
    Q <- abs(ss0-ss)
    ss0 <-ss
    if(iter >=max.iter) break
  }
  gpa.final <- pGpa(slid, PrinAxes = F, Proj = F)
  list(coords=gpa.final$coords, consensus=gpa.final$consensus, iter=iter+iter.s+1, Q=Q)
}

# pGPA.wSliders
# GPA with partial Procrustes superimposition, incorporating semilandmarks
# used in gpagen
pGpa.wSliders <- function(Y, curves, surf, ProcD = TRUE, PrinAxes = FALSE, Proj = FALSE, max.iter = 5){
  n <- length(Y); p <- nrow(Y[[1]]); k <- ncol(Y[[1]])
  Yc <- Map(function(y) center.scale(y), Y)
  CS <- sapply(Yc,"[[","CS")
  Ya <- lapply(Yc,"[[","coords")
  Ya <- apply.pPsup(Ya[[1]], Ya)
  M <- Reduce("+", Ya)/n
  if(ProcD == FALSE) gpa.slide <- BE.slide(curves, surf, Ya, ref=M, max.iter=max.iter) else 
    gpa.slide <- procD.slide(curves, surf, Ya, ref=M, max.iter=max.iter)
  Ya <- gpa.slide$coords
  M <- gpa.slide$consensus
  iter <- gpa.slide$iter
  Q <- gpa.slide$Q
  if (PrinAxes == TRUE) {
    ref <- M
    rot <- prcomp(ref)$rotation
    for (i in 1:k) if (sign(rot[i, i]) != 1) 
      rot[1:k, i] = -rot[1:k, i]
    Ya <- Map(function(y) y%*%rot, Ya)
    M <- center.scale(Reduce("+", Ya)/n)$coords
  }
  list(coords= Ya, CS=CS, iter=iter, iter.s=NULL, consensus=M, Q=Q, nsliders=NULL)
}

# description still needed
#
#
tps<-function(matr, matt, n,sz=1.5, pt.bg="black",
              grid.col="black", grid.lwd=1, grid.lty=1, refpts=FALSE){		#DCA: altered from J. Claude: 2D only	
  xm<-min(matt[,1])
  ym<-min(matt[,2])
  xM<-max(matt[,1])
  yM<-max(matt[,2])
  rX<-xM-xm; rY<-yM-ym
  a<-seq(xm-1/5*rX, xM+1/5*rX, length=n)
  b<-seq(ym-1/5*rX, yM+1/5*rX,by=(xM-xm)*7/(5*(n-1)))
  m<-round(0.5+(n-1)*(2/5*rX+ yM-ym)/(2/5*rX+ xM-xm))
  M<-as.matrix(expand.grid(a,b))
  ngrid<-tps2d(M,matr,matt)
  plot(ngrid, cex=0.2,asp=1,axes=FALSE,xlab="",ylab="")
  for (i in 1:m){lines(ngrid[(1:n)+(i-1)*n,], col=grid.col,lwd=grid.lwd,lty=grid.lty)}
  for (i in 1:n){lines(ngrid[(1:m)*n-i+1,], col=grid.col,lwd=grid.lwd,lty=grid.lty)}
  if(refpts==FALSE) points(matt,pch=21,bg=pt.bg,cex=sz) else points(matr,pch=21,bg=pt.bg,cex=sz)
}

# description still needed
#
#
tps2d<-function(M, matr, matt)
{p<-dim(matr)[1]; q<-dim(M)[1]; n1<-p+3
P<-matrix(NA, p, p)
for (i in 1:p)
{for (j in 1:p){
  r2<-sum((matr[i,]-matr[j,])^2)
  P[i,j]<- r2*log(r2)}}
P[which(is.na(P))]<-0
Q<-cbind(1, matr)
L<-rbind(cbind(P,Q), cbind(t(Q),matrix(0,3,3)))
m2<-rbind(matt, matrix(0, 3, 2))
coefx<-qr.solve(qr(L))%*%m2[,1]
coefy<-qr.solve(qr(L))%*%m2[,2]
fx<-function(matr, M, coef)
{Xn<-numeric(q)
for (i in 1:q)
{Z<-apply((matr-matrix(M[i,],p,2,byrow=TRUE))^2,1,sum)
Xn[i]<-coef[p+1]+coef[p+2]*M[i,1]+coef[p+3]*M[i,2]+sum(coef[1:p]*(Z*log(Z)))}
Xn}
matg<-matrix(NA, q, 2)
matg[,1]<-fx(matr, M, coefx)
matg[,2]<-fx(matr, M, coefy)
matg}

# description still needed
#
#
tps2d3d<-function(M, matr, matt){		#DCA: altered from J. Claude 2008  
  p<-dim(matr)[1]; k<-dim(matr)[2];q<-dim(M)[1]
  Pdist<-as.matrix(dist(matr))
  ifelse(k==2,P<-Pdist^2*log(Pdist^2),P<- Pdist) 
  P[which(is.na(P))]<-0
  Q<-cbind(1, matr)
  L<-rbind(cbind(P,Q), cbind(t(Q),matrix(0,k+1,k+1)))
  m2<-rbind(matt, matrix(0, k+1, k))   
  coefx<-qr.solve(qr(L))%*%m2[,1]
  coefy<-qr.solve(qr(L))%*%m2[,2]
  if(k==3){coefz<-qr.solve(qr(L))%*%m2[,3]}
  fx<-function(matr, M, coef){
    Xn<-numeric(q)
    for (i in 1:q){
      Z<-apply((matr-matrix(M[i,],p,k,byrow=TRUE))^2,1,sum)  
      ifelse(k==2,Z1<-Z*log(Z),Z1<-sqrt(Z)); Z1[which(is.na(Z1))]<-0
      ifelse(k==2,Xn[i]<-coef[p+1]+coef[p+2]*M[i,1]+coef[p+3]*M[i,2]+sum(coef[1:p]*Z1),
             Xn[i]<-coef[p+1]+coef[p+2]*M[i,1]+coef[p+3]*M[i,2]+coef[p+4]*M[i,3]+sum(coef[1:p]*Z1))
    }    
    Xn}
  matg<-matrix(NA, q, k)
  matg[,1]<-fx(matr, M, coefx)
  matg[,2]<-fx(matr, M, coefy)
  if(k==3){matg[,3]<-fx(matr, M, coefz)}  
  matg
}

# geomorph.data.frame
# data frames for shape data and independent variables
# used in all functions with 'procD.lm"
geomorph.data.frame <- function(...) {
  dots <- list(...)
  list.check0 <- sapply(1:length(dots), function(j) any(is.geomorph.data.frame(dots[[j]])))
  dots0 <- unlist(dots[list.check0], recursive=FALSE)
  dots.updated <- dots[!list.check0]
  if(length(dots.updated) > 0) {
    list.check1 <- sapply(1:length(dots.updated), function(j) is.gpagen(dots.updated[[j]]))
    dots1 <- dots.updated[list.check1]
  } else dots1 <- NULL
  if(length(dots1) > 0){
    dots2 <- lapply(1:length(dots1), function(j){
      x <- unlist(dots1[j], recursive = FALSE)
      list(x$coords, x$Csize)
    })
    dots2 <- unlist(dots2, recursive = FALSE)
    dots2.names <- rep(c("coords", "Csize"), length(dots1))
    names(dots2) <- dots2.names
  } else dots2 <- NULL
  dots.updated <- dots.updated[!list.check1]
  if(length(dots.updated) > 0){
    list.check2 <- sapply(1:length(dots.updated), function(j) is.phylo(dots.updated[[j]]))
    dots3 <- dots.updated[list.check2]
    dots4 <- dots.updated[!list.check2]
  } else {
    dots3 <- NULL
    dots4 <- NULL
  }
  if(length(dots3) == 0) dots3 <- NULL
  if(length(dots4) == 0) dots4 <- NULL
  dots <- c(dots0,dots2, dots3,dots4)
  N <- length(dots)
  dots.ns <- array(NA,N)
  for(i in 1:N){
    if(is.array(dots[[i]])) {
      if(length(dim(dots[[i]])) == 3) dots.ns[i] <- dim(dots[[i]])[[3]]
      if(length(dim(dots[[i]])) == 2) dots.ns[i] <- dim(dots[[i]])[[2]]
      if(length(dim(dots[[i]])) == 1) dots.ns[i] <- dim(dots[[i]])[[1]]
    }
    if(is.matrix(dots[[i]])) dots.ns[i] <- dim(dots[[i]])[[1]]
    if(class(dots[[i]]) == "dist") dots.ns[i] <- attr(dots[[i]], "Size")
    if(class(dots[[i]]) == "phylo") dots.ns[i] <- length(dots[[i]]$tip.label)
    if(is.data.frame(dots[[i]])) dots.ns[i] <- dim(dots[[i]])[[2]]
    if(is.vector(dots[[i]])) dots.ns[i] <- length(dots[[i]])
    if(is.factor(dots[[i]])) dots.ns[i] <- length(dots[[i]])
    if(is.logical(dots[[i]])) dots.ns[i] <- length(dots[[i]])
  }
  if(any(is.na(dots.ns))) stop("Some input is either dimensionless or inappropriate for data frames")
  if(length(unique(dots.ns)) > 1) stop("Inputs have different numbers of observations")
  class(dots) <- c("geomorph.data.frame",class(dots))
  dots
}

# pcoa
# acquires principal coordimates from distance matrices
# used in all functions with 'procD.lm" via porcD.fit
pcoa <- function(D){
  options(warn=-1)
  if(class(D) != "dist") stop("function only works with distance matrices")
  cmd <- cmdscale(D, k=attr(D, "Size") -1, eig=TRUE)
  options(warn=0)
  p <- length(cmd$eig[zapsmall(cmd$eig) > 0])
  Yp <- cmd$points[,1:p]
  Yp
}

# description still needed
#
#
simplify2data.frame <- function(g, n) { # g = geomorph.data.frame class object
  newdf <- list(ind=rep(1,n))
  df.names <- names(g)
  for(i in 1:length(df.names)) {
    if(is.array(g[[i]])) newdf[[i+1]] <- two.d.array(g[[i]])
    if(is.matrix(g[[i]])) newdf[[i+1]] <- g[[i]]
    if(is.vector(g[[i]])) newdf[[i+1]] <- g[[i]]
    if(is.factor(g[[i]])) newdf[[i+1]] <- g[[i]]
  }
  newdf <- newdf[-1]
  names(newdf) <- df.names
  as.data.frame(newdf)
}

# procD.fit
# lm-like fit modified for Procrustes residuals 
# general workhorse for all 'procD.lm' functions
# used in all 'procD.lm' functions
procD.fit <- function(f1, keep.order=FALSE, pca=TRUE, data=NULL,...){
  form.in <- formula(f1)
  if(any(class(f1)=="lm")) {
    weights <- f1$weights
    contrasts <- f1$contrasts
    offset <-f1$offset
    data <- model.frame(f1)
    Y <- as.matrix(data[1])
  } else {
    if(!is.null(data)) {
      dat <- data 
    Y <- dat[which(!is.na(match(names(dat),as.character(form.in[[2]]))))][[1]]
    } else {
      Y <- eval(form.in[[2]], parent.frame())
      dat <- model.frame(form.in[-2])
    }
    
    if(class(Y) == "dist") Y <- pcoa(Y) else
      if(length(dim(Y)) == 3)  Y <- two.d.array(Y) else 
        Y <- as.matrix(Y)
    weights <- NULL
    contrasts <- NULL
    offset <- NULL
  }
  n <- nrow(Y)
  if(ncol(Y) > n & pca==TRUE){
    y.pca <- prcomp(Y)
    Y <- y.pca$x[,zapsmall(y.pca$sdev) > 0]
  }
  dots <- list(...)
  if(!is.null(contrasts)) {
    op.c <- options()$contrasts
    options(contrasts = unlist(contrasts))
    contrasts <- contrasts
  } else if(!is.null(dots$contrasts)){
    op.c <- options()$contrasts
    options(contrasts = unlist(dots$contrasts))
    contrasts <- dots$contrasts
  }
  if (any(is.na(Y)) == T) stop("Response data matrix (shape) contains missing values. Estimate these first (see 'estimate.missing').")
  form.terms <- form.in
  form.terms[[2]] <- NULL
  df.check <- c("matrix", "vector", "numeric", "integer", "character", "logical", "factor")
  df.check <- sapply(1:length(dat), function(j) any(class(dat[[j]]) == df.check))
  mf <- as.data.frame(dat[df.check])
  Terms <- terms(form.terms, keep.order=keep.order, data = mf)
  if(nrow(Y) != nrow(mf)) stop("Different numbers of specimens in dependent and independent variables")
  X <- model.matrix(Terms, data=mf)
  if(is.null(weights) & !is.null(dots$weights)) w <- dot$weights else 
    if(is.null(weights)) w <- rep(1,n)
  if(sum(w)==n){
    wY <- Y; wX <- X
  } else {
    wY <- Y*sqrt(w); wX <- X*sqrt(w)
  }
  X.k <- attr(X, "assign")
  k <- length(X.k)
  QRx <- qr(X)
  X <- X[, QRx$pivot, drop = FALSE]
  X <- X[, 1:QRx$rank, drop = FALSE]
  X.k <- X.k[QRx$pivot][1:QRx$rank]
  uk <- unique(X.k)
  k <- length(uk) - 1
  Xs <- lapply(1:length(uk), function(j)  Xj <- X[, X.k %in% uk[1:j]])
  QRs <- lapply(Xs, function(x) qr(x))
  fitted <- lapply(QRs, function(x) qr.fitted(x,Y))
  residuals <- lapply(QRs, function(x) qr.resid(x,Y))
  coefficients <- lapply(QRs, function(x) qr.coef(x,Y))
  if(sum(w)==n){
    wXs <- Xs; wQRs <- QRs; wFitted <- fitted; wResiduals <- residuals
    wCoefficients <- coefficients
  } else {
    wXs <- lapply(Xs, function(x) x*w)
    wQRs <- lapply(wXs, function(x) qr(x))
    wFitted <- lapply(wQRs, function(x) qr.fitted(x,Y))
    wResiduals <- lapply(wQRs, function(x) qr.resid(x,Y))
    wCoefficients <- lapply(wQRs, function(x) qr.coef(x,Y))
  }
  
  if(!is.null(offset)) {
    fitted <- fitted + offset
    wFitted <- wFitted + offset
    residuals <- residuals + offset
    wResiduals <- wResiduals + offset
  } else if (!is.null(dots$offset)) {
    fitted <- fitted + dots$offset
    wFitted <- wFitted + dots$offset
    residuals <- residuals + dots$offset
    wResiduals <- wResiduals + dots$offset
  }
  term.labels <- attr(Terms, "term.labels")
  if(length(term.labels) > 0) mf.out <- model.frame(Terms, data= mf) else
    mf.out <- data.frame(Int = rep(1,n))
  mf.out <- data.frame(Y=Y, mf.out)
  out <- list(Y=Y, wY=wY, X=X, Xs=Xs, wX=wX, wXs=wXs,
              QRs = QRs, wQRs=wQRs, fitted=fitted, wFitted=wFitted,
              residuals = residuals, wResiduals=wResiduals,
              coefficients=coefficients, wCoefficients=wCoefficients,
              weights = w, data = mf.out, 
              contrasts = contrasts,
              Terms = Terms, term.labels = term.labels)
  class(out) <- "procD.fit"
  invisible(out)
}

# perm.index
# creates a permutation index for resampling
# used in all functions with a resampling procedure
perm.index <-function(n,iter){
  set.seed(iter)
  ind <- c(list(1:n),(Map(function(x) sample(1:n), 1:iter)))
  ind
}

# boot.index
# creates a bootstrap index for resampling
# used in modularity test functions
boot.index <-function(n,iter){
  set.seed(iter)
  ind <- c(list(1:n),(Map(function(x) sample(1:n, replace = TRUE), 1:iter)))
  ind
}

# SS.iter
# calculates SS in random iterations of a resmapling procedure
# used in nearly all 'procD.lm' functions, unless pgls in used
SS.iter = function(pfit,iter, Yalt="RRPP"){
  Y <- as.matrix(pfit$Y)
  k <- length(pfit$QRs)
  n <- nrow(Y)
  Yh <- pfit$fitted
  E <- pfit$residuals
  w<- pfit$weights
  Xr <- lapply(pfit$wXs[1:(k-1)], function(x) as.matrix(x))
  Xf <- lapply(pfit$wXs[2:k], function(x) as.matrix(x))
  ind = perm.index(n,iter)
  if(Yalt=="RRPP") {
    if(sum(w)==n) {
      Yr = Map(function(x) (Map(function(y,e) e[x,]+y, Yh[1:(k-1)], E[1:(k-1)])),ind)
    } else {
      Yr = Map(function(x) (Map(function(y,e) (e[x,]+y)*sqrt(w), Yh[1:(k-1)], E[1:(k-1)])),ind) 
    }} else {
      if(sum(w)==n) {
        Yr = Map(function(x) Map(function(y) y[x,], lapply(1:(k-1),function(.) Y)),ind)
      } else {
        Yr = Map(function(x) Map(function(y) (y[x,])*sqrt(w), lapply(1:(k-1),function(.) Y)),ind)
      }
    }
  SS <- sapply(1:(iter+1), function(j){
    mapply(function(x1,x2,y) sum(.lm.fit(x1,y)$residuals^2 - .lm.fit(x2,y)$residuals^2), 
           Xr, Xf,Yr[[j]])})
  SS
}

# Fpgls.iter
# calculates F values in random iterations of a resmapling procedure, with pgls involved
# used in the 'procD.lm' functions where pgls in used

Fpgls.iter = function(pfit,Pcor,iter, Yalt="RRPP"){
  Y <- as.matrix(pfit$Y)
  k <- length(pfit$QRs)
  n <- nrow(Y)
  Yh <- pfit$fitted
  E <- pfit$residuals
  w<- pfit$weights
  wQRs <- pfit$wQRs
  dfE <- sapply(1:k, function(j) wQRs[[j]]$rank)
  df <- dfE[-1] - dfE[1:(k-1)]
  Pcor <- Pcor[rownames(Y),rownames(Y)]
  PwXs <- lapply(pfit$wXs, function(x) crossprod(Pcor,as.matrix(x)))
  Xr <- lapply(PwXs[1:(k-1)], function(x) as.matrix(x))
  Xf <- lapply(PwXs[2:k], function(x) as.matrix(x))
  ind = perm.index(n,iter)
  if(Yalt=="RRPP") {
    Yr = Map(function(x) (Map(function(y,e) crossprod(Pcor,as.matrix(e[x,]+y)*sqrt(w)), Yh[1:(k-1)], E[1:(k-1)])),ind) 
  } else {
    Yr = Map(function(x) Map(function(y) crossprod(Pcor,as.matrix((y[x,])*sqrt(w))), lapply(1:(k-1),function(.) Y)),ind)
  }
  SS <- lapply(1:(iter+1), function(j){
    mapply(function(x1,x2,y) sum(.lm.fit(x1,y)$residuals^2 - .lm.fit(x2,y)$residuals^2), 
           Xr, Xf,Yr[[j]])})
  SSEs <- Map(function(y) sum(.lm.fit(Xf[[k-1]],y[[k-1]])$residuals^2), Yr)
  Fs <- mapply(function(s1,s2) (s1/df)/(s2/(n-k)), SS, SSEs)
  list(SS=simplify2array(SS), Fs=Fs)
}

# anova.parts
# makes an ANOVA table based on SS from SS.iter
# used in nearly all 'procD.lm' functions
anova.parts <- function(pfit, SS){ # SS from SS.iter
  Y <- pfit$wY
  k <- length(pfit$term.labels)
  dfe <-sapply(pfit$wQRs, function(x) x$rank)
  df <- dfe[2:(k+1)] - dfe[1:k]
  if(k==1) SS <- SS[1] else SS <- SS[,1]
  anova.terms <- pfit$term.labels
  SSY <- sum(qr.resid(qr(pfit$wX[,1]),pfit$wY)^2)
  MS <- SS/df
  R2 <- SS/SSY
  SSE.model <- SSY - sum(SS)
  dfE <- nrow(Y)-(sum(df)+1)
  MSE <- SSE.model/dfE
  Fs <- MS/MSE
  df <- c(df,dfE,nrow(Y)-1)
  SS <- c(SS,SSE.model, SSY)
  MS <- c(MS,MSE,NA)
  R2 <- c(R2,NA,NA)
  Fs <- c(Fs,NA,NA)
  anova.tab <- data.frame(df,SS,MS,Rsq=R2,F=Fs)
  rownames(anova.tab) <- c(anova.terms, "Residuals", "Total")
  out <- list(anova.table = anova.tab, anova.terms = anova.terms,
              SS = SS, df = df, R2 = R2[1:k], F = Fs[1:k])
  out
}

# anova.parts.symmetry
# makes an ANOVA table based on SS from SS.iter
# used in symmetry analysis, which is a nested model
anova.parts.symmetry <- function(pfit, SS,object.sym){ # SS from SS.iter
  Y <- pfit$wY
  k <- length(pfit$term.labels)
  dfe <-sapply(pfit$wQRs, function(x) x$rank)
  df <- dfe[2:(k+1)] - dfe[1:k]
  if(k==1) SS <- SS[1] else SS <- SS[,1]
  anova.terms <- pfit$term.labels
  SSY <- sum(qr.resid(qr(pfit$wX[,1]),pfit$wY)^2)
  if(object.sym==TRUE){
    SS<-SS/2;SSY<-SSY/2
  }
  MS <- SS/df
  R2 <- SS/SSY
  Fs <- array(NA,k)
  Fs[1]<-MS[1]/MS[3]; Fs[2]<-MS[2]/MS[3]
  if(k==4) Fs[3]<-MS[3]/MS[4] else Fs[3]=NA  
  anova.tab <- data.frame(df,SS,MS,Rsq=R2,F=Fs)
  rownames(anova.tab) <- anova.terms
  out <- list(anova.table = anova.tab, anova.terms = anova.terms,
              SS = SS, df = df, R2 = R2[1:k], F = Fs[1:k])
  out
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

# pval
# Effect sizes (standard deviates) form random outcomes
# any analytical function
effect.size <- function(x, center = FALSE) {
  z = scale(x, center=center)
  n <- length(z)
  z[1]*sqrt((n-1)/(n))
}

# Effect.size.matrix
# Effect sizes form random outcomes that comprise matrices
# any analytical function with results in matrices
Effect.size.matrix <- function(M, center=F){
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

# single.factor
# converts factorial designs to single-factor variables
# advanced.procD.lm
single.factor <- function(pfit) {# pfit = Procrustes fit
  Terms <- pfit$Terms
  dat <- pfit$data
  datClasses <- sapply(dat, function(x) data.class(x))
  facs <- dat[which(datClasses == "factor")]
  facs <- as.data.frame(facs)
  if(ncol(facs) > 1) fac <- factor(apply(facs, 1,function(x) paste(x, collapse=":"))) else 
    fac <- as.factor(unlist(facs))
  fac
}

# cov.extract
# Extacts covariates from design matrices
# advanced.procD.lm
cov.extract <- function(pfit) {
  Terms <- pfit$Terms
  dat <- pfit$data
  datClasses <- sapply(dat, function(x) data.class(x))
  datClasses <- which(datClasses == "numeric")
  covs <- dat[datClasses]
  covs
}

# ls.means
# estimates ls.means from models with slopes and factors
# advanced.procD.lm
ls.means = function(pfit, Y=NULL, g=NULL, data=NULL) { 
  if(is.null(Y)) Y <- pfit$wY
  n <- nrow(Y)
  if(is.null(data)) dat <- pfit$data else dat <- data
  if(!is.null(g)) {
    if(ncol(as.matrix(g)) > 1){
      fac.check <- sapply(g, is.factor)
      facs <- g[,fac.check]
    } else if(is.factor(g)) facs <- g else stop("groups input neither a factor nor factors")
    if(ncol(as.matrix(facs)) > 1) fac <- factor(apply(facs, 1,function(x) paste(x, collapse=":"))) else 
      fac <- as.factor(unlist(facs))
  } else fac <- single.factor(pfit)
  covs <- cov.extract(pfit)
  if(length(covs) == 0){
    lsm <- .lm.fit(model.matrix(~fac+0),Y)$coefficients
  } else if(ncol(as.matrix(covs)) == 1){
    covs <- as.vector(covs[[1]])
    fit <- .lm.fit(model.matrix(~covs+fac+0),Y)
    X <- cbind(mean(covs),model.matrix(~fac+0))
    lsm <- .lm.fit(model.matrix(~fac+0),X%*%coef(fit))$coefficients
  } else {
    covs <- sapply(covs, function(x) matrix(x))
    fit <- .lm.fit(model.matrix(~covs+fac+0),Y)
    X <- cbind(matrix(rep(colMeans(covs), rep.int(nrow(covs), ncol(covs))),n),
               model.matrix(~fac+0))
    lsm <- .lm.fit(model.matrix(~fac+0),X%*%coef(fit))$coefficients
  }
  rownames(lsm) <- levels(fac)
  lsm
}

# apply.ls.means
# estimates ls.means from models with slopes and factors across random outcomes in a list
# advanced.procD.lm
apply.ls.means <- function(pfit, Yr, g=NULL, data=NULL){
  if(is.null(data)) dat <- pfit$data else dat <- data
  if(!is.null(g)) lsm <- Map(function(y) ls.means(pfit, Y=y, g=g, data=dat), Yr) else
    lsm <- Map(function(y) ls.means(pfit, Y=y, g=NULL, data=dat), Yr)
  lsm
}

# slopes
# estimates slopes from models with slopes and factors
# advanced.procD.lm
slopes = function(pfit, Y=NULL, g = NULL, slope=NULL, data = NULL){ 
  if(is.null(Y)) Y <- pfit$wY
  if(is.null(data)) dat <- pfit$data else dat <- data
  if(!is.null(g)) {
    if(ncol(as.matrix(g)) > 1){
      fac.check <- sapply(g, is.factor)
      facs <- g[,fac.check]
    } else if(is.factor(g)) facs <- g else stop("groups input neither a factor nor factors")
    if(ncol(as.matrix(facs)) > 1) fac <- factor(apply(facs, 1,function(x) paste(x, collapse=":"))) else 
      fac <- as.factor(unlist(facs))
  } else fac <- single.factor(pfit)
  if(!is.null(slope)) covs <- as.matrix(slope) else covs <- as.matrix(cov.extract(pfit))
  if(ncol(covs) > 1) stop("Only one covariate can be used for slope-comparisons")
  if(ncol(covs) == 0) stop("No covariate for which to compare slopes")
  B <- qr.coef(qr(model.matrix(~fac*covs)), Y)
  fac.p <- qr(model.matrix(~fac))$rank
  Bslopes <- as.matrix(B[-(1:fac.p),])
  Bslopes[2:nrow(Bslopes),] <- Bslopes[2:nrow(Bslopes),] + Bslopes[1,]
  if(ncol(Y)==1) Bslopes <- cbind(1, Bslopes)
  rownames(Bslopes) <- levels(fac)
  Bslopes
}

# apply.slopes
# estimates slopes from models with slopes and factors across random outcomes in a list
# advanced.procD.lm
apply.slopes <- function(pfit, Yr, g=NULL, slope=NULL, data=NULL){
  Y <- pfit$wY
  if(is.null(data)) dat <- pfit$data else dat <- data
  if(!is.null(g)) 
    slopes <- Map(function(y) slopes(pfit, Y=y, g=g, slope=slope, data=dat), Yr) else
      slopes <- Map(function(y) slopes(pfit, Y=y, g=NULL, slope=slope, data=dat), Yr)
  if(ncol(Y)==1) slopes <- Map(function(s) cbind(1,s), slopes)
  slopes 
}

# vec.cor.matrix
# the vector correlations among multiple slopes
# advanced.procD.lm
vec.cor.matrix <- function(M) {
  M = as.matrix(M)
  w = 1/sqrt(diag(tcrossprod(M)))
  vc = tcrossprod(M*w)
  options(warn = -1)
  vc
}

# vec.ang.matrix
# converts the vector correlations from vec.cor.matrix to angles
# advanced.procD.lm
vec.ang.matrix <- function(M, type = c("rad", "deg", "r")){
  M= as.matrix(M)
  type= match.arg(type)
  if(type == "r") {
    vc = vec.cor.matrix(M)
  } else {
    vc = vec.cor.matrix(M)
    vc = acos(vc)
    diag(vc) = 0
  }
  if(type == "deg") vc = vc*180/pi
  vc
}

# description still needed
#
#
pls <- function(x,y, RV=FALSE, verbose = FALSE){
  x <- as.matrix(x); y <- as.matrix(y)
  px <- ncol(x); py <- ncol(y); pmin <- min(px,py)
  S <-var(cbind(x,y))
  S12 <- matrix(S[1:px,-(1:px)], px,py)
  pls <- .Internal(La_svd("S", S12, double(pmin), matrix(0,px,pmin), matrix(0,pmin,py)))
  U <- pls$u; V <- t(pls$vt)
    XScores <- x %*% U
    YScores <- y %*% V
  r.pls <- cor(XScores[,1],YScores[,1])
  if(RV==TRUE){
    S11 <- S[1:px,1:px]
    S22 <- S[-(1:px),-(1:px)]
    RV <- sum(colSums(S12^2))/sqrt(sum(S11^2)*sum(S22^2))
  } else
      RV <- NULL
  if(verbose==TRUE){
    XScores <- as.matrix(XScores); Y <- as.matrix(YScores)
    rownames(U)  = colnames(x); rownames(V) = colnames(y)
    out <- list(pls.svd = pls, r.pls = r.pls, RV=RV, left.vectors=U, 
                right.vectors=V, XScores=XScores,YScores=YScores)
  } else out <- r.pls
  out
}

# description still needed
#
#
quick.pls <- function(x,y, px, py, pmin) {# no RV; no verbose output
  # assume parameters already found
  S <-var(cbind(x,y))
  S12 <- matrix(S[1:px,-(1:px)], px,py)
  pls <- .Internal(La_svd("S", S12, double(pmin), matrix(0,px,pmin), matrix(0,pmin,py)))
  U<-pls$u; V <- t(pls$vt)
    XScores <- x %*% U
    YScores <- y %*% V
  cor(XScores[,1],YScores[,1])
}

# description still needed
#
#
apply.pls <- function(x,y, RV=FALSE, iter){
  x <- as.matrix(x); y <- as.matrix(y)
  px <- ncol(x); py <- ncol(y)
  pmin <- min(px,py)
  ind <- perm.index(nrow(x), iter)
  y.rand <-lapply(1:(iter+1), function(j) y[ind[[j]],])
  if(RV == TRUE) RV.rand <- sapply(1:(iter+1), function(j) pls(x,y.rand[[j]], RV=TRUE, verbose = TRUE)$RV) else
    r.rand <- sapply(1:(iter+1), function(j) quick.pls(x,y.rand[[j]], px,py,pmin))
  if(RV == TRUE) RV.rand else r.rand
}

# description still needed
#
#
plsmulti<-function(x,gps){
  gp.names <- levels(gps)
  g<-factor(as.numeric(gps))
  ngps<-nlevels(g)
  S<-var(x)
  gps.combo <- combn(ngps, 2)
  pls.gp <- sapply(1:ncol(gps.combo), function(j){ # no loops
    S12<-S[which(g==gps.combo[1,j]),which(g==gps.combo[2,j])]
    px <- nrow(S12); py <- ncol(S12); pmin <- min(px,py)
    pls<-.Internal(La_svd("S", S12, double(pmin), matrix(0,px,pmin), matrix(0,pmin,py)))
    U<-pls$u; V<-t(pls$vt)
    XScores<-x[,which(g==gps.combo[1,j])]%*%U[,1]; YScores<-x[,which(g==gps.combo[2,j])]%*%V[,1]
    cor(XScores,YScores)
  })
  if(length(pls.gp) > 1) pls.mat <- dist(matrix(0, ngps,)) else 
    pls.mat = 0 
  for(i in 1:length(pls.mat)) pls.mat[[i]] <- pls.gp[i]
  pls.obs <- mean(pls.gp) 
  list(r.pls = pls.obs, r.pls.mat=pls.mat)
}

# description still needed
#
#
quick.plsmulti <- function(x,gps){
  gp.names <- levels(gps)
  g<-factor(as.numeric(gps))
  ngps<-nlevels(g)
  S<-var(x)
  gps.combo <- combn(ngps, 2)
  pls.gp <- sapply(1:ncol(gps.combo), function(j){ # no loops
    S12<-S[which(g==gps.combo[1,j]),which(g==gps.combo[2,j])]
    px <- nrow(S12); py <- ncol(S12); pmin <- min(px,py)
    pls<-.Internal(La_svd("S", S12, double(pmin), matrix(0,px,pmin), matrix(0,pmin,py)))
    U<-pls$u; V<-t(pls$vt)
    XScores<-x[,which(g==gps.combo[1,j])]%*%U[,1]; YScores<-x[,which(g==gps.combo[2,j])]%*%V[,1]
    cor(XScores,YScores)
  })
  mean(pls.gp) 
}

# description still needed
#
#
apply.plsmulti <- function(x,gps, iter){
  ngps<-nlevels(gps)
  S <-var(x)
  r.obs <- plsmulti(x,gps)$r.pls
  ind <- perm.index(nrow(x), iter)
  x.r<-lapply(1:(iter+1), function(j) x[ind[[j]],which(gps==levels(gps)[1])]) #shuffle 1st block
  r.rand<-sapply(1:(iter+1), function(j) quick.plsmulti(cbind(x.r[[j]],
                                                            x[,which(gps!=levels(gps)[1])]),gps)) 
  r.rand
}

# description still needed
#
#
CR<-function(x,gps){
  gp.names <- levels(gps)
  g<-factor(as.numeric(gps))
  ngps<-nlevels(g)
  S<-var(x)
  diag(S)<-0
  gps.combo <- combn(ngps, 2)
  CR.gp <- sapply(1:ncol(gps.combo), function(j){ # no loops
    S11<-S[which(g==gps.combo[1,j]),which(g==gps.combo[1,j])]
    S22<-S[which(g==gps.combo[2,j]),which(g==gps.combo[2,j])]
    S12<-S[which(g==gps.combo[1,j]),which(g==gps.combo[2,j])]
    sqrt(sum(colSums(S12^2))/sqrt(sum(S11^2)*sum(S22^2)))
  })
  if(length(CR.gp) > 1) CR.mat <- dist(matrix(0, length(CR.gp),)) else 
    CR.mat = 0 # may not be necessary
  for(i in 1:length(CR.mat)) CR.mat[[i]] <- CR.gp[i]
  
  CR.obs <- mean(CR.gp) 
  list(CR = CR.obs, CR.mat=CR.mat)
}

# description still needed
#
#
quick.CR <-function(x,gps){ # no CR.mat made
  gp.names <- levels(gps)
  g<-factor(as.numeric(gps))
  ngps<-nlevels(g)
  S<-var(x)
  diag(S)<-0
  gps.combo <- combn(ngps, 2)
  CR.gp <- sapply(1:ncol(gps.combo), function(j){
    S11<-S[which(g==gps.combo[1,j]),which(g==gps.combo[1,j])]
    S22<-S[which(g==gps.combo[2,j]),which(g==gps.combo[2,j])]
    S12<-S[which(g==gps.combo[1,j]),which(g==gps.combo[2,j])]
    sqrt(sum(colSums(S12^2))/sqrt(sum(S11^2)*sum(S22^2)))
  })
  mean(CR.gp)
}

# description still needed
#
#
apply.CR <- function(x,gps, iter){
  CR.obs <- CR(x,gps)$CR
  ind <- perm.index(length(gps), iter)
  x.r<-lapply(1:(iter+1), function(j) x[,ind[[j]]])
  CR.rand<-sapply(1:(iter+1), function(j) quick.CR(x.r[[j]],gps)) 
  CR.rand
}

# description still needed
#
#
boot.CR <- function(x,gps, iter){
  x<-as.matrix(x)
  boot <- boot.index(nrow(x)-1, iter)
  x.r<-lapply(1:(iter+1), function(j) x[boot[[j]],])
  CR.boot<-sapply(1:(iter+1), function(j) quick.CR(x.r[[j]],gps)) 
  CR.boot
}

# description still needed
#
#
CR.phylo<-function(x,invC,D.mat,gps){
  gp.names <- levels(gps)
  g<-factor(as.numeric(gps))
  ngps<-nlevels(g)
  one<-matrix(1,nrow(x),1)  
  a<-colSums(invC %*% x)*sum(invC)^-1  
  R<- t(x-one%*%a)%*%invC%*%(x-one%*%a)*(nrow(x)-1)^-1 
  diag(R)<-0
  gps.combo <- combn(ngps, 2)
  CR.gp <- sapply(1:ncol(gps.combo), function(j){ 
    R11<-R[which(g==gps.combo[1,j]),which(g==gps.combo[1,j])]
    R22<-R[which(g==gps.combo[2,j]),which(g==gps.combo[2,j])]
    R12<-R[which(g==gps.combo[1,j]),which(g==gps.combo[2,j])]
    sqrt(sum(colSums(R12^2))/sqrt(sum(R11^2)*sum(R22^2)))
  })
  if(length(CR.gp) > 1) CR.mat <- dist(matrix(0, length(CR.gp),)) else 
    CR.mat = 0 
  for(i in 1:length(CR.mat)) CR.mat[[i]] <- CR.gp[i]
  
  CR.obs <- mean(CR.gp) 
  list(CR = CR.obs, CR.mat=CR.mat)
}

# description still needed
#
#
quick.CR.phylo <- function(x,invC,gps){
  x <- as.matrix(x); invC <- as.matrix(invC)
  gp.names <- levels(gps)
  g<-factor(as.numeric(gps))
  ngps<-nlevels(g)
  one<-matrix(1,nrow(x),1)  
  a<-colSums(invC %*% x)*sum(invC)^-1  
  R<- t(x-one%*%a)%*%invC%*%(x-one%*%a)*(nrow(x)-1)^-1 
  diag(R)<-0
  gps.combo <- combn(ngps, 2)
  CR.gp <- sapply(1:ncol(gps.combo), function(j){ 
    R11<-R[which(g==gps.combo[1,j]),which(g==gps.combo[1,j])]
    R22<-R[which(g==gps.combo[2,j]),which(g==gps.combo[2,j])]
    R12<-R[which(g==gps.combo[1,j]),which(g==gps.combo[2,j])]
    sqrt(sum(colSums(R12^2))/sqrt(sum(R11^2)*sum(R22^2)))
  })
  
  mean(CR.gp) 
}

# description still needed
#
#
apply.phylo.CR <- function(x,invC,gps, iter){
  ind <- perm.index(length(gps), iter)
  x.r<-lapply(1:(iter+1), function(j) x[,ind[[j]]])
  CR.rand<-sapply(1:(iter+1), function(j) quick.CR.phylo(x.r[[j]],invC=invC,gps=gps)) 
  CR.rand
}

# description still needed
#
#
boot.phylo.CR <- function(x,invC,gps, iter){
  x<-as.matrix(x)
  boot <- boot.index(nrow(x)-1, iter)
  x.r<-lapply(1:(iter+1), function(j) x[boot[[j]],])
  invC.r <- lapply(1:(iter+1), function(j) invC[boot[[j]],boot[[j]]])
  CR.boot<-sapply(1:(iter+1), function(j) quick.CR.phylo(x=x.r[[j]],invC=invC.r[[j]],gps=gps)) 
  CR.boot
}

# description still needed
#
#
phylo.mat<-function(x,phy){
  C<-vcv.phylo(phy,anc.nodes=FALSE) 
  C<-C[rownames(x),rownames(x)] 
  invC <-fast.solve(C)
  eigC <- eigen(C)
  lambda <- zapsmall(eigC$values)
  if(any(lambda == 0)){
    warning("Singular phylogenetic covariance matrix. Proceed with caution")
    lambda = lambda[lambda > 0]
  }
  eigC.vect = eigC$vectors[,1:(length(lambda))]
  D.mat <- qr.solve(eigC.vect*sqrt(lambda),eigC.vect) 
  list(invC = invC, D.mat = D.mat,C = C)
}

#####-----------------------------------------------------------------------------------
### All print functions

print.procD.lm <- function (x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call), "\n\n")
  cat("\nType I (Sequential) Sums of Squares and Cross-products\n")
  if(x$perm.method == "RRPP") cat ("Randomized Residual Permutation Procedure Used\n") else
    cat("Randomization of Raw Values used\n")
  cat(paste(x$permutations, "Permutations"))
  cat("\n\n")
  print(x$aov.table)
  invisible(x)
}

print.advanced.procD.lm <- function (x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call), "\n\n")
  cat("\nRandomized Residual Permutation Procedure Used\n")
  cat(paste(x$permutations, "Permutations"))
  cat("\nANOVA Table")
  cat("\n\n")
  print(x$anova.table); cat("\n\n")
  if(!is.null(x$LS.means)) {cat("LS means\n"); print(x$LS.means); cat("\n")}
  if(!is.null(x$slopes)) {cat("Slopes\n");print(x$slopes); cat("\n\n")}
  if(!is.null(x$LS.means.dist)) {cat("LS means distance matrix\n");print(x$LS.means.dist); cat("\n")}
  if(!is.null(x$Z.means.dist)) {cat("Effect sizes (Z)\n");print(x$Z.means.dist); cat("\n")}
  if(!is.null(x$P.means.dist)) {cat("P-values\n");print(x$P.means.dist); cat("\n\n")}
  if(!is.null(x$slopes.dist)) {cat("Contrasts in slope vector length\n");print(x$slopes.dist); cat("\n")}
  if(!is.null(x$Z.slopes.dist)) {cat("Effect sizes (Z)\n");print(x$Z.slopes.dist); cat("\n")}
  if(!is.null(x$P.slopes.dist)) {cat("P-values\n");print(x$P.slopes.dist); cat("\n\n")}
  if(!is.null(x$slopes.cor)) {cat("Correlations between slope vectors\n");print(x$slopes.cor); cat("\n")}
  if(!is.null(x$Z.slopes.cor)) {cat("Effects sizes (Z)\n");print(x$Z.slopes.cor); cat("\n")}
  if(!is.null(x$P.slopes.cor)) {cat("P-values\n");print(x$P.slopes.cor); cat("\n\n")}
  if(!is.null(x$slopes.angles)) {cat("Angles between slope vectors\n");print(x$slopes.angles); cat("\n")}
  if(!is.null(x$Z.angles)) {cat("Effects sizes (Z)\n");print(x$Z.angles); cat("\n")}
  if(!is.null(x$P.angles)) {cat("P-values\n");print(x$P.angles); cat("\n\n")}
  invisible(x)
}

print.morphol.disparity <- function (x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call), "\n\n")
  cat("\nRandomized Residual Permutation Procedure Used\n")
  cat(paste(x$permutations, "Permutations\n"))
  cat("\nProcrustes variances for defined groups\n")
  print(x$Procrustes.var)
  cat("\n")
  cat("\nPairwise absolute differences between variances\n")
  print(x$PV.dist)
  cat("\n")
  cat("\nP-Values\n")
  print(x$PV.dist.Pval)
  cat("\n\n")
  invisible(x)
}

print.pls <- function (x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call), "\n\n")
  if(x$method=="RV") {
    cat(paste("\nRV:", round(x$RV, nchar(x$permutations)-1)))
    cat(paste("\n\nP-value:", round(x$P.value, nchar(x$permutations)-1)))
    cat(paste("\n\nBased on", x$permutations, "random permutations"))
  }
  if(x$method=="PLS") {
    cat(paste("\nr-PLS:", round(x$r.pls, nchar(x$permutations)-1)))
    cat(paste("\n\nP-value:", round(x$P.value, nchar(x$permutations)-1)))
    cat(paste("\n\nBased on", x$permutations, "random permutations"))
  }
  invisible(x)
}

print.gpagen <- function (x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call), "\n\n")
  cat("\nGeneralized Procrustes Analysis\n")
  cat("with Partial Procrustes Superimposition\n\n")
  cat(paste(x$p-x$nsliders, "fixed landmarks\n"))
  cat(paste(x$nsliders, "semilandmarks (sliders)\n"))
  cat(paste(x$k,"-dimensional landmarks\n",sep=""))
  cat(paste(x$iter, "GPA iterations to converge\n"))
  if(!is.null(x$slide.method)) sm <- match.arg(x$slide.method, c("BE", "ProcD")) else
    sm <- "none"
  if(sm == "ProcD") cat("Minimized squared Procrustes Distance used\n")
  if(sm == "BE") cat("Minimized Bending Energy used\n")
  cat("\n\nConsensus (mean) Configuration\n\n")
  print(x$consensus)
  invisible(x)
}

print.bilat.symmetry <- function (x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call), "\n\n")
  cat(paste("Symmetry (data) type:", x$data.type), "\n")
  cat("\nType I (Sequential) Sums of Squares and Cross-products\n")
  if(x$perm.method == "RRPP") cat ("Randomized Residual Permutation Procedure Used\n") else
    cat("Randomization of Raw Values used\n")
  cat(paste(x$permutations, "Permutations"))
  cat("\n\nShape ANOVA\n")
  print(x$shape.anova)
  if(x$data.type == "Matching") {
    cat("\n\nCentroid Size ANOVA\n")
    print(x$size.anova)
  }
  invisible(x)
}

print.procD.lm <- function (x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call), "\n\n")
  cat("\nType I (Sequential) Sums of Squares and Cross-products\n")
  if(x$perm.method == "RRPP") cat ("Randomized Residual Permutation Procedure Used\n") else
    cat("Randomization of Raw Values used\n")
  cat(paste(x$permutations, "Permutations"))
  cat("\n\n")
  print(x$aov.table)
  invisible(x)
}

print.CR <- function (x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call), "\n\n")
  cat(paste("\nCR:", round(x$CR, nchar(x$permutations))))
  cat(paste("\n\nP-value:", round(x$P.value, nchar(x$permutations))))
  cat(paste("\n\nBased on", x$permutations, "random permutations"))
  cat(paste("\n\nConfidence Intervals", round(x$CInterval,nchar(x$permutations))))
  invisible(x)
}

print.physignal <- function(x){
  cat("\nCall:\n")
  cat(deparse(x$call), "\n\n")
  cat(paste("\nObserved Phylogenetic Signal (K):", round(x$phy.signal, nchar(x$permutations))))
  cat(paste("\n\nP-value:", round(x$pvalue, nchar(x$permutations))))
  cat(paste("\n\nBased on", x$permutations, "random permutations"))
  invisible(x)
}

#####-----------------------------------------------------------------------------------

### All plot functions

plot.gpagen <- function(g, ...){
  plotAllSpecimens(g$coords)
}

plot.het <- function(r,f){
  r <- center(r)
  f <- center(f)
  r <- sqrt(diag(tcrossprod(r)))
  f <- sqrt(diag(tcrossprod(f)))
  lfr <- loess(r~f)
  lfr <- cbind(lfr$x, lfr$y, lfr$fitted)
  lfr <- lfr[order(lfr[,1]),]
  plot(lfr, pch=19, asp=1, 
       xlab = "Procrustes Distance Fitted Values",
       ylab = "Procrustes Distance Residuals", 
       main = "Residuals vs. Fitted")
  points(lfr[,1], lfr[,3], type="l", col="red")
}

plot.QQ <- function(r){
  r <- center(r)
  r <- sqrt(diag(tcrossprod(r)))
  r <- sort(r)
  n <- length(r)
  tq <- (seq(1,n)-0.5)/n
  tq <- qnorm(tq)
  plot(tq, r, pch=19, xlab = "Theoretical Quantiles",
       ylab = "Procrustes Distance Residuals", 
       main = "Q-Q plot")
}

plot.procD.lm <- function(x, outliers=FALSE){
  r <- x$wResiduals
  f <- x$wFitted
  pca.r <- prcomp(r)
  var.r <- round(pca.r$sdev^2/sum(pca.r$sdev^2)*100,2)
  plot(pca.r$x, pch=19, asp =1,
       xlab = paste("PC 1", var.r[1],"%"),
       ylab = paste("PC 2", var.r[2],"%"),
       main = "PCA Residuals")
  pca.f <- prcomp(f)
  var.f <- round(pca.f$sdev^2/sum(pca.f$sdev^2)*100,2)
  dr <- sqrt(diag(tcrossprod(center(r))))
  plot.QQ(r)
  plot(pca.f$x[,1], dr, pch=19, asp =1,
       xlab = paste("PC 1", var.f[1],"%"),
       ylab = "Procrustes Distance Residuals",
       main = "Residuals vs. PC 1 fitted")
  lfr <- loess(dr~pca.f$x[,1])
  lfr <- cbind(lfr$x, lfr$fitted); lfr <- lfr[order(lfr[,1]),]
  points(lfr, type="l", col="red")
  plot.het(r,f)
  p <- ncol(r)
  if(outliers==TRUE){
    if(p/3 == round(p/3)) ra <- arrayspecs(r,p/3,3) else 
      ra <- arrayspecs(r,p/2,2)
    plotOutliers(ra)
  }
}

plot.advanced.procD.lm <- function(x) plot.procD.lm(x)

plot.pls <- function(p, label = NULL, warpgrids=TRUE){
  A1 <- p$A1; A2 <- p$A2
  XScores <- p$XScores; YScores <- p$YScores
  if (length(dim(A1)) == 3) {
    A1.ref <- mshape(A1)
    pls1.min <- A1[, , which.min(XScores[, 1])]
    pls1.max <- A1[, , which.max(XScores[, 1])]
  }
  if (length(dim(A2)) == 3) {
    A2.ref <- mshape(A2)
    pls2.min <- A2[, , which.min(XScores[, 1])]
    pls2.max <- A2[, , which.max(XScores[, 1])]
  }
    if (length(dim(A1)) != 3 && length(dim(A2)) != 3) {
      plot(XScores[, 1], YScores[, 1], pch = 21, bg = "black", 
           main = "PLS Plot", xlab = "PLS1 Block 1", ylab = "PLS1 Block 2")
      if (length(label != 0)) {
        text(XScores[, 1], YScores[, 1], label, adj = c(-0.7, -0.7))
      }
    }
    if (length(dim(A1)) == 3 || length(dim(A2)) == 3) {
      
      par(mar = c(1, 1, 1, 1) + 0.1)
      split.screen(matrix(c(0.22, 1, 0.22, 1, 0.19, 0.39, 0, 
                            0.19, 0.8, 1, 0, 0.19, 0, 0.19, 0.19, 0.39, 0, 0.19, 
                            0.8, 1), byrow = T, ncol = 4))
      screen(1)
      plot(XScores[, 1], YScores[, 1], pch = 21, bg = "black", 
           main = "PLS1 Plot: Block 1 (X) vs. Block 2 (Y) ", 
           xlab = "PLS1 Block 1", ylab = "PLS1 Block 2")
      if (length(label != 0)) {
        text(XScores[, 1], YScores[, 1], label, adj = c(-0.7, 
                                                        -0.7))    
      }
      if (warpgrids == TRUE) {
        if (length(dim(A1)) == 3 && dim(A1)[2] == 2) {
          screen(2)
          tps(A1.ref, pls1.min, 20, sz = 0.7)
          screen(3)
          tps(A1.ref, pls1.max, 20, sz = 0.7)
        }
        if (length(dim(A2)) == 3 && dim(A2)[2] == 2) {
          screen(4)
          tps(A2.ref, pls2.min, 20, sz = 0.7)
          screen(5)
          tps(A2.ref, pls2.max, 20, sz = 0.7)
        }
      }
      close.screen(all.screens = TRUE)
      par(mar = c(5.1, 4.1, 4.1, 2.1))
    }
    if (length(dim(A1)) == 3 && dim(A1)[2] == 3) {
      plot(XScores[, 1], YScores[, 1], pch = 21, bg = "black", 
           main = "PLS Plot", xlab = "PLS1 Block 1", ylab = "PLS1 Block 2")
      if (length(label != 0)) {
        text(XScores[, 1], YScores[, 1], label, adj = c(-0.7, 
                                                        -0.7))
      }
      open3d()
      plot3d(pls1.min, type = "s", col = "gray", main = paste("PLS Block1 negative"), 
             size = 1.25, aspect = FALSE)
      open3d()
      plot3d(pls1.max, type = "s", col = "gray", main = paste("PLS Block1 positive"), 
             size = 1.25, aspect = FALSE)
    }
    if (length(dim(A2)) == 3 && dim(A2)[2] == 3) {
      open3d()
      plot3d(pls2.min, type = "s", col = "gray", main = paste("PLS Block2 negative"), 
             size = 1.25, aspect = FALSE)
      open3d()
      plot3d(pls2.max, type = "s", col = "gray", main = paste("PLS Block2 positive"), 
             size = 1.25, aspect = FALSE)
    } 
}

plot.bilat.symmetry <- function(b, warpgrids = TRUE, mesh= NULL){
  k <- dim(b$symm.component)[[2]]
  if(b$data.type == "Matching"){
    if(k==2){  
      par(mfrow=c(2,2),oma=c(1.5,0,1.5,0))
      plotAllSpecimens(b$symm.component)
      plotAllSpecimens(b$asymm.component)
      plotRefToTarget(b$DA.mns[,,1],b$DA.mns[,,2],method="TPS",main="Directional Asymmetry")
      plotRefToTarget(b$FA.mns[,,1],b$FA.mns[,,2],method="TPS",main="Fluctuating Asymmetry")
      mtext("Symmetric Shape Component (left) and Asymmetric Shape Component (right)",outer = TRUE,side=3)
      mtext("Mean directional (left) and fluctuating (right) asymmetry",side = 1, outer = TRUE)
      par(mfrow=c(1,1))
    }
    if (k==3){
      if (is.null(mesh)){
        open3d()
        plotRefToTarget(b$DA.mns[,,1],b$DA.mns[,,2],method="points",main="Directional Asymmetry")
        open3d()
        plotRefToTarget(b$FA.mns[,,1],b$FA.mns[,,2],method="points",main="Fluctuating Asymmetry")
      } 
      if(!is.null(mesh)){
        plotRefToTarget(b$DA.mns[,,1],b$DA.mns[,,2],mesh,method="surface")
        title3d(main="Directional Asymmetry")
        plotRefToTarget(b$FA.mns[,,1],b$FA.mns[,,2],mesh,method="surface")
        title3d(main="Fluctuating Asymmetry")
      }
    }
    layout(1) 
  }
  if(b$data.typ == "Object"){
    if(warpgrids==TRUE){
      if(k==2){  
        par(mfrow=c(2,2),oma=c(1.5,0,1.5,0))
        plotAllSpecimens(b$symm.component)
        plotAllSpecimens(b$asymm.component)
        plotRefToTarget(b$DA.mns[,,1],b$DA.mns[,,2],method="TPS",main="Directional Asymmetry")
        plotRefToTarget(b$FA.mns[,,1],b$FA.mns[,,2],method="TPS",main="Fluctuating Asymmetry")
        mtext("Symmetric Shape Component (left) and Asymmetric Shape Component (right)",outer = TRUE,side=3)
        mtext("Mean directional (left) and fluctuating (right) asymmetry",side = 1, outer = TRUE)
      }
      if (k==3){
        if(is.null(mesh)) {
          open3d()
          plotRefToTarget(b$DA.mns[,,1],b$DA.mns[,,2],method="points",main="Directional Asymmetry")
          open3d()
          plotRefToTarget(b$FA.mns[,,1],b$FA.mns[,,2],method="points",main="Fluctuating Asymmetry")
        } 
        if(!is.null(mesh)){
          plotRefToTarget(b$DA.mns[,,1],b$DA.mns[,,2],mesh,method="surface")
          title3d(main="Directional Asymmetry")
          plotRefToTarget(b$FA.mns[,,1],b$FA.mns[,,2],mesh,method="surface")
          title3d(main="Fluctuating Asymmetry")
        }  
      }
      layout(1) 
    } 
  }
}

plot.CR <- function(cr){
  CR.val <- cr$random.CR
  CR.obs <- cr$CR
  p <- cr$P.value
  ndec <- nchar(p)-2
  CR.obs <- round(CR.obs, ndec)
  main.txt <- paste("Observed CR =",CR.obs,";", "P-value =", p)
  hist(CR.val,30,freq=TRUE,col="gray",xlab="CR Coefficient",xlim=c(0,max(c(2,CR.val))),
       main=main.txt, cex.main=0.8)
  arrows(CR.obs,50,CR.obs,5,length=0.1,lwd=2)
}

plot.physignal <- function(ps){
  K.val <- ps$random.K
  K.obs <- ps$phy.signal
  p <- ps$pvalue
  ndec <- nchar(p)-2
  K.obs <- round(K.obs, ndec)
  main.txt <- paste("Observed K =",K.obs,";", "P-value =", p)
  hist(K.val,30,freq=TRUE,col="gray",xlab="Phylogenetic Signal, K",
       main=main.txt, cex.main=0.8)
  arrows(K.obs,50,K.obs,5,length=0.1,lwd=2)
}

plot.evolrate <- function(ER){
  Rate.val <- ER$random.sigma
  Rate.obs <- ER$random.sigma[1]
  p <- ER$P.value
  ndec <- nchar(ER$permutations)
  Rate.obs <- round(Rate.obs, ndec)
  p <- round(p, ndec)
  main.txt <- paste("Observed Rate Ratio =",Rate.obs,";", "P-value =", p)
  hist(Rate.val,30,freq=TRUE,col="gray",xlab="Rate Ratios",xlim=c(0,max(c(2,Rate.val))),
       main=main.txt, cex.main=0.8)
  arrows(Rate.obs,50,Rate.obs,5,length=0.1,lwd=2)
}

#####-----------------------------------------------------------------------------------

### All summary functions (mostly the same as print functions)

summary.gpagen <- function(x) print.gpagen(x)
summary.procD.lm <- function(x) print.procD.lm(x)
summary.advanced.procD.lm <- function(x) print.advanced.procD.lm(x)
summary.bilat.symmetry <- function(x) print.bilat.symmetry(x)
summary.morphol.disparity <- function(x) print.morphol.disparity(x)
summary.pls <- function(x) print.pls(x)
summary.CR <- function(x) print.CR(x)

#####-----------------------------------------------------------------------------------

### geomorph-specific logicals

is.gpagen <- function(x) class(x) == "gpagen"
is.phylo <- function(x) class(x) == "phylo"
is.geomorph.data.frame <- function(x) class(x) == "geomorph.data.frame"

#####-----------------------------------------------------------------------------------

### retained from old geomorph support code
### need to update and merge, or replace with new functions


scan.to.ref<-function(scandata,specland,refland){  	#DCA
  ref.scan<-tps2d3d(scandata,specland,refland)
  ref.scan}

refscan.to.spec<-function(refscan,refland,specland){ 	#DCA
  unwarp.scan<-tps2d3d(refscan,refland,specland)
  unwarp.scan}


#' Estimate mean shape for a set of aligned specimens
#'
#' Estimate the mean shape for a set of aligned specimens
#'
#' The function estimates the average landmark coordinates for a set of aligned specimens. It is assumed 
#' that the landmarks have previously been aligned using Generalized Procrustes Analysis (GPA) 
#'  [e.g., with \code{\link{gpagen}}]. This function is described in Claude (2008).
#'
#' @param A An array (p x k x n) containing GPA-aligned coordinates for a set of specimens
#' @keywords utilities
#' @export
#' @author Julien Claude 
#' @references Claude, J. 2008. Morphometrics with R. Springer, New York.
#' @examples
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment   
#'
#' mshape(Y.gpa$coords)   #mean (consensus) configuration
mshape<-function(A){apply(A,c(1,2),mean)}	


# Trajectory Size: Pathlength Distance
pathdist<-function(M) {as.matrix(dist(M))} 
trajsize<-function(M,n,p){
  traj.pathdist<-array(0,dim=c(n,1))   		
  for (i in 1:n){
    temp<-pathdist(M[,,i])
    for (j in 1:(p-1)){
      traj.pathdist[i]<-traj.pathdist[i]+temp[j,j+1]
    }
  }
  traj.size.dist<-as.matrix(dist(traj.pathdist))		
}

# Trajectory Orientation
trajorient<-function(M,n,k){
  traj.orient<-array(NA,dim=c(n,k))   
  check.1<-array(NA,dim=c(n))
  for (i in 1:n){
    temp<-svd(var(M[,,i]))$v[1:k,1]
    traj.orient[i,]<-temp
    check.1[i]<-M[1,,i]%*%traj.orient[i,]  
    check.1[i]<-check.1[i]/abs(check.1[i])
    if(check.1[i]==-1) traj.orient[i,]<--1*traj.orient[i,]
  }
  options(warn=-1)				
  traj.ang.diff<-(180/pi)*acos(traj.orient%*%t(traj.orient))
}

# Trajectory Shape
trajshape<-function(M){
  x<-pgpa(M)
  traj.shape.dist<-as.matrix(x$intereucl.dist) 
}

# general plotting function for phenotypic trajectories
trajplot<-function(Data,M, groups, group.cols = NULL, ...){
  n<-dim(M)[3]; p<-dim(M)[1]
  pmax <- max(Data[,1]); pmin <- min(Data[,1])
  plot(Data[,1:2],type="n",
       xlim = c(2*pmin, pmax),
       xlab="PC I", ylab="PC II",
       main="Two Dimensional View  of Phenotypic Trajectories",asp=1)
  if(!is.null(group.cols)) col.index <- group.cols else col.index <-1:length(levels(groups))
  if(length(groups) == n) {
    gp.index <- levels(groups)
    col.temp <-array(,length(groups))
    for(i in 1:length(col.temp)) col.temp[i] = col.index[which(match(gp.index,gp.index[gp.index==groups[i]])==1)]
    col.index=col.temp
  }
  
  points(Data[,1:2],pch=21,bg="gray",cex=.75)
  for (i in 1:n){  	 	
    for (j in 1:(p-1)){		
      points(M[(j:(j+1)),1,i],M[(j:(j+1)),2,i],type="l",pch=21,col=col.index[i], lwd=2)  #was black    
    }
  }
  for (i in 1:n){		 	
    for (j in 2:(p-1)){		
      points(M[j,1,i],M[j,2,i],pch=21,bg="gray",col="black",cex=1.5)
    }
  }
  for (i in 1:n){
    points(M[1,1,i],M[1,2,i],pch=21,bg="white",col="black",cex=1.5)
  }
  for (i in 1:n){
    points(M[p,1,i],M[p,2,i],pch=21,bg="black",col="black",cex=1.5)
  }
  legend("topleft", levels(groups), lwd=2, col=levels(as.factor(col.index)))
}

# Write .nts file for output of digitize2d(), buildtemplate() digit.fixed() and digitsurface()
# A is an nx2 or nx3 matrix of the output coordinates. To be used internally only.

writeland.nts <- function(A, spec.name, comment=NULL){
  ntsfile=paste(spec.name,".nts",sep="")
  file.create(file=ntsfile)
  if(is.null(comment)){
    cat(paste('"',spec.name,sep=""),file= ntsfile,sep="\n",append=TRUE)
  }
  else if(!is.null(comment)){
    cat(paste('"',spec.name,sep=""),file= ntsfile,sep="\n")
    cat(paste('"',comment,sep=""),file= ntsfile,sep="\n",append=TRUE)
  }
  dims <- dim(A)
  if (dims[2] == 2){
    cat(paste(1,dims[1],2,0,"dim=2"),file= ntsfile,sep="\n",append=TRUE)
  }
  else if (dims[2] == 3){
    cat(paste(1,dims[1],3,0, "dim=3"),file= ntsfile,sep="\n",append=TRUE)
  }
  write.table(A ,file= ntsfile,col.names = FALSE, row.names = FALSE,sep="  ",append=TRUE)
}

# picscale is called by digitize2d
picscale<- function(scale){
  digscale<-NULL
  digscale<-locator(2,type="o",lwd=2,col="red",lty="11")
  cat(paste("Keep scale (y/n)?"), "\n")
  ans <- readLines(n = 1)
  if (ans == "n") {
    cat(paste("Set scale again"), "\n")
  }
  while (ans == "n") {
    digscale<-NULL
    digscale<-locator(2,type="o",lwd=2,col="red",lty="11")
    cat(paste("Keep scale (y/n)?"), "\n")
    ans <- readLines(n = 1)
    if (ans == "y") { 
    }
    if (ans == "n") {
      cat(paste("Set scale again"), "\n")
    }
  }
  scale/sqrt(sum(diff(digscale$x)^2+diff(digscale$y)^2))      
}

# Function written by person who wrote identify() - called by define.modules
identifyPch <- function(x, y = NULL, n = length(x), pch = 19, col="red", ...)
{
  xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
  sel <- rep(FALSE, length(x)); res <- integer(0)
  while(sum(sel) < n) {
    ans <- identify(x[!sel], y[!sel], n = 1, plot = FALSE, ...)
    if(!length(ans)) break
    ans <- which(!sel)[ans]
    points(x[ans], y[ans], pch = pch, col=col)
    sel[ans] <- TRUE
    res <- c(res, ans)
  }
  res
} 

# allometry data frame coerces covariate = Size
allometry.data.frame <- function(f1){
  dat <- procD.data.frame(f1)
  if(!is.vector(dat[[2]])) stop("First formula must contain only a single size covariate")
  names(dat)[[2]] <- "Size"
  dat
}

# function for generating random SS for submodels, using RRPP for trajectories only
random.trajectories <- function(pf, Yalt = c("resample", "RRPP"), iter, pca=TRUE){ # like anova.parts, but faster for resampling
  k <- length(pf$Terms)
  Y <- as.matrix(pf$Y)
  if(pca==TRUE) Y<-prcomp(Y)$x else Y <- Y
  n <- nrow(Y)
  p<-ncol(Y) 
  if(!is.null(pf$weights)) w <- as.vector(pf$weights) else w <-rep(1,n)
  if(any(w < 0)) stop("Weights cannot be negative")
  Xs <- pf$Xs
  dat <- pf$mf
  Terms <- terms(dat)
  facs <- dat[which(attr(Terms,"dataClasses")=="factor")]
  if(ncol(facs) != 2) stop("Model must contain two factors")
  if(k < 3) stop("Model does not appear to be factorial.  Check model formula (see help file).") 
  int.term<-grep(":", pf$Terms[k])
  if(int.term!=1) stop("Last col of X-matrix does not contain interaction between main effects (see help file).")        
  n1<-length(levels(facs[,1]))
  k1<-length(levels(facs[,2]))
  fac12<-single.factor(pf$call)
  covs <- cov.extract(pf$call)
  Plm <- PSize <- POrient <-PShape <- as.list(array(,iter+1))
  E <- Yh <- Yhw <- as.list(array(,k+1))
  for(i in 1:(k+1)){
    wfit <- lmfit(as.matrix(Xs[[i]]),Y)
    Yhw[[i]] <- as.matrix(wfit$fitted)
    Yh[[i]] <- as.matrix(wfit$fitted/w)
    E[[i]] <- as.matrix(wfit$residuals)
  }
  if(ncol(covs) == 0) covs <- NULL else covs <- model.matrix(~covs)
  lsmeans.obs <- ls.means(fac12, cov.mf=covs, Y)
  traj.specs.obs<- aperm(array(t(lsmeans.obs), c(p,k1,n1)), c(2,1,3)) 
  trajsize.obs<-trajsize(traj.specs.obs,n1,k1) 
  trajdir.obs<-trajorient(traj.specs.obs,n1,p); diag(trajdir.obs)<-0 
  trajshape.obs<-trajshape(traj.specs.obs) 
  SSEs.obs <- SSE(E)
  ind <- perm.index(n,iter)
  if(iter > 0) {for(i in 1:(iter+1)){
    Er <- Map(function(x) x[ind[[i]],], E)
    Yr <- as.list(array(,k+1))
    if(Yalt == "RRPP") for(ii in 1:(k+1)) Yr[[ii]] <- Reduce("+",list(Er[[ii]], Yh[[ii]])) else
      Yr <- lapply(Yr,function(x) as.matrix(Y[ind[[i]],]))
    Yr <- lapply(Yr, function(x) as.matrix(x)*w)
    SSEs.null <- SSE(mod.resids(Xs,Yr))
    SSEs.r <- SSE(mod.resids(Xs[-1],Yr[1:k]))
    lsmeans.r <- ls.means(fac12, cov.mf=covs, as.matrix(Yr[[k]]))
    traj.specs.r<- aperm(array(t(lsmeans.r), c(p,k1,n1)), c(2,1,3)) 
    trajsize.r<-trajsize(traj.specs.r,n1,k1) 
    trajdir.r<-trajorient(traj.specs.r,n1,p); diag(trajdir.r)<-0 
    trajshape.r<-trajshape(traj.specs.r) 
    Plm[[i]] <- SSEs.null[1:k]-SSEs.r
    PSize[[i]] <- trajsize.r
    POrient[[i]] <- trajdir.r
    PShape[[i]] <- trajshape.r
  }}
  if(iter > 0) Y=Yr[[k]] else Y=Y
  list(Plm=Plm, PSize=PSize, POrient=POrient,PShape=PShape,Y=Y, 
       traj.pts = k1, trajectories = traj.specs.obs, groups = facs[,1])
}


Gower.center <- function(D, calc.dist=FALSE){
  D <- as.matrix(D)
  if(calc.dist == TRUE) D = as.matrix(dist(D))
  n <- nrow(D)
  A <- -0.5*D^2
  Id <- diag(1,n)
  one <- matrix(1,n)
  G <- (Id-(1/n)*one%*%t(one))%*%A%*%(Id-(1/n)*one%*%t(one))
  G
}

Hat.SSE <- function(G, X){
  I <- diag(1,nrow(G))
  H <- X%*%solve(t(X)%*%X)%*%t(X)
  SS <- sum(diag((I-H)%*%G%*%(I-H)))
  SS
}

Hat.SS.model <- function(G, X){
  I <- diag(1,nrow(G))
  H <- X%*%solve(t(X)%*%X)%*%t(X)
  SS <- sum(diag(H%*%G%*%H))
  SS
}

Hat.anova.tab <- function(D, f1, keep.order=TRUE){ # assumes dependent is distance matrix
  form.in <- formula(f1)
  Terms <- terms(form.in, keep.order = keep.order)
  G <- Gower.center(D)
  n <- nrow(D)
  I <- diag(1,n)
  Xn <- matrix(1,n)
  X <- model.matrix(form.in)
  SSE <-Hat.SSE(G,X)
  SSM <-Hat.SS.model(G,X)
  SST <-Hat.SSE(G,Xn)
  dfM <- qr(X)$rank - 1
  dfE <- n - qr(X)$rank
  dfT <- n -1 
  MSM <- SSM/dfM
  MSE <- SSE/dfE
  Fs <- MSM/MSE
  R2 <- SSM/SST
  df <- c(dfM,dfE,dfT)
  SS <- c(SSM,SSE,SST)
  MS <- c(MSM,MSE,NA)
  R2 <- c(R2,NA,NA)
  Fs <- c(Fs,NA,NA)
  a.tab <- data.frame(df,SS,MS,Rsq=R2,F=Fs)
  rownames(a.tab) <- c(attr(Terms, "term.labels"), "Residuals", "Total")
  a.tab
}

