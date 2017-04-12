#' @name geomorph-package
#' @docType package
#' @aliases geomorph
#' @title Geometric morphometric analyses for 2D/3D data
#' @author Dean C. Adams, Michael Collyer, Antigoni Kaliontzopoulou & Emma Sherratt
#' 
#' @description Functions in this package allow one to read, manipulate, and digitize landmark data; generate shape
#'  variables via Procrustes analysis for points, curves and surface data, perform statistical analyses
#'  of shape variation and covariation, and provide graphical depictions of shapes and patterns of
#'  shape variation.
#'  
#' @import ape
#' @import rgl
#' @import stats
#' @import utils
#' @import graphics
#' @import grDevices
#' @importFrom geiger sim.char
#' @importFrom jpeg readJPEG
#' @importFrom Matrix nearPD
#' 
#' @section geomorph TOC:
#' geomorph-package
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

#' Landmarks on mosquito wings
#'
#' @name mosquito
#' @docType data
#' @author Dean Adams
#' @keywords datasets
NULL

#' Landmarks on pupfish
#'
#' @name pupfish
#' @docType data
#' @author Michael Collyer
#' @keywords datasets
#' @description Landmark data from Cyprindon pecosensis body shapes, with indication of Sex and 
#' Population from which fish were sampled (Marsh or Sinkhole).  
#' @details These data were previously aligned 
#' with GPA.  Centroid size (CS) is also provided.
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic 
#' change for phenotypes described by high-dimensional data. Heredity. 113: doi:10.1038/hdy.2014.75.
NULL

#' Tail-shapes of larval salamanders
#'
#' @name larvalTails
#' @docType data
#' @author Michael Collyer
#' @keywords datasets
#' @description Landmark data from tails of larval Salamanders exposed to different treatments of herbicides. 
#' @details  
#' Data set includes tail landmarks (landmarks), index for identifying semilandmarks (sliders), and vectors 
#' for variables including herbicide treatment (Treatment) and Family (Family).  The latter variable indicates
#' the egg masses (clutches) sampled in the wild from which individual eggs were randomly assigned to treatment.  
#' See Levis et al. (2016) for more experimental details.
#' @references Levis, N.A, M.L. Schooler, J.R. Johnson, and M.L. Collyer. 2016. The effects of terrestrial and aquatic herbicides on 
#' larval salamander morphology and swim speed. Biological Journal of the Linnean Society.  Accepted.
NULL

#' Estimate mean shape for a set of aligned specimens
#'
#' Estimate the mean shape for a set of aligned specimens
#'
#' The function estimates the average landmark coordinates for a set of aligned specimens. It is assumed 
#' that the landmarks have previously been aligned using Generalized Procrustes Analysis (GPA) 
#'  [e.g., with \code{\link{gpagen}}]. This function is described in Claude (2008).
#'
#' @param A Either a list (length n, ach p x k), A 3D array (p x k x n), or a matrix (pk X n) containing GPA-aligned coordinates for a set of specimens
#' @keywords utilities
#' @export
#' @author Julien Claude 
#' @references Claude, J. 2008. Morphometrics with R. Springer, New York.
#' @examples
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment   
#'
#' mshape(Y.gpa$coords)   #mean (consensus) configuration
mshape<-function(A){
  if(is.array(A)) res <- apply(A,c(1,2),mean)
  if(is.list(A)) res <- Reduce("+", A)/length(A)
  if(is.matrix(A)) res <- colMeans(A)
  if(!is.array(A) && !is.list(A) && !is.matrix(A)) stop("There are not multiple configurations from which to obtain a mean.")
  return(res)
}	

#####----------------------------------------------------------------------------------------------------

# SUPPORT FUNCTIONS

# center
# centers a matrix faster than scale()
# used in other functions for gpagen; digitsurface
center <- function(x){
  if(is.vector(x)) x- mean(x) else {
    x <- as.matrix(x)
    x - rep(colMeans(x), rep.int(nrow(x), ncol(x)))
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

# orp
# projection in GPA
# used in gpagen functions
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

# rotate.mat
# simple rotation matrix via svd
# digitsurface
rotate.mat <- function(M,Y){
  k <- ncol(M)
  M <- cs.scale(M); Y <- cs.scale(Y)
  MY <- crossprod(M,Y)
  sv <- La.svd(MY,k,k)
  u <- sv$u; u[,k] <- u[,k]*determinant(MY)$sign
  v <- t(sv$vt)
  tcrossprod(v,u)
}

# apply.pPsup
# applies a partial Procrustes superimposition to matrices in a list
# used in gpagen functions
apply.pPsup<-function(M, Ya) {	# M = mean (reference); Ya all Y targets
  k <- ncol(Ya[[1]]); p <- nrow(Ya[[1]]); n <- length(Ya)
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
  Xsvd <- La.svd(X, k, k)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  rtu <-((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
  v <-t(Xsvd$vt)[, Positive, drop = FALSE]
  v%*%rtu
}

# fast.solve
# chooses between fast.ginv or qr.solve, when det might or might not be 0
# used in any function requiring a matrix inverse where the certainty of
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
    Up1 <- Up2 <- array(0,dim=c(k*p,p))  
    PC <- getSurfPCs(y, surf)
    z11 <- z12 <- cbind(surf,surf); z21 <- z22 <- cbind(p+surf, surf)
    if(k==3) z31 <- z32 <- cbind(2*p+surf, surf)
    pc11 <- PC$p1x; pc12 <- PC$p1y; pc21 <- PC$p2x; pc22 <- PC$p2y
    diag(Up1[1:p,1:p]) <- pc11; diag(Up2[1:p,1:p]) <- pc21
    diag(Up1[(1+p):(2*p),1:p]) <- pc12; diag(Up2[(1+p):(2*p),1:p]) <- pc22
    if(k==3) {pc13 <- PC$p1z; pc23 <- PC$p2z
    diag(Up1[(1+2*p):(3*p),1:p]) <- pc13 
    diag(Up2[(1+2*p):(3*p),1:p]) <- pc23}
  }                   
  U <- cbind(U,Up1,Up2)    

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
  Q <- cbind(1,Mr)
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
  iter <- 0
  pb <- txtProgressBar(min = 0, max = max.iter, initial = 0, style=3) 
  setTxtProgressBar(pb,iter)
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
  setTxtProgressBar(pb,iter)
  while(Q > 0.0001){
    iter <- iter+1
    Ya <- apply.pPsup(M, Ya)
    M <- Reduce("+",Ya)/n
    ss2 <- n*(1-sum(M^2))
    Q <- abs(ss-ss2)
    ss <- ss2
    M <- cs.scale(M)
    setTxtProgressBar(pb,iter)
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
  if(iter < max.iter) setTxtProgressBar(pb,max.iter)
  close(pb)
  list(coords= Ya, CS=CS, iter=iter, consensus=M, Q=Q, nsliders=NULL)
}

# .pGPA
# same as pGPA, but without progress bar option
# used in gpagen
.pGpa <- function(Y, PrinAxes = FALSE, Proj = FALSE, max.iter = 5){
  iter <- 0
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
  V <- La.svd(center(y), nu=0)$vt
  p <- nrow(y); k <- ncol(y)
  pc.match <- 1:p; pc.match[-surf] = NA
  nearpts <- lapply(1:p, function(j) {
    nn <- pc.match[j]
    if(is.na(nn)) 0 else 
      c(nearest(y,nn, k=k+1),nn)})
  tmp.pts <- lapply(1:p, function(j) {
    k <- nearpts[[j]]
    if(sum(k) > 0) x <- center(y[k,]) else x <- NA
    x})
  pc.dir <- lapply(1:p, function(j) {
    x <- tmp.pts[[j]]
    if(is.matrix(x)) {
      pc <- La.svd(x, nu=0)$vt
      s=sign(diag(crossprod(V,pc)))
      pc*s
      } else 0
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

# Based on the equation y -U%*%solve(crossprod(U))%*%crossprod(U,(y-ref))
#                          LP             MP                RP
# left part (LP), middle part(MP) and right part (RP) are accomplished faster in the code in this function.

semilandmarks.slide.tangents.procD <- function(y,tans, ref){
  yc <- y - ref
  p <- nrow(yc); k <-ncol(yc)
  if(k==3) {ycx <- yc[,1]; ycy <- yc[,2]; ycz <- yc[,3 ]} else {ycx <- yc[,1]; ycy <- yc[,2]}
  if(k==3) {tx <- tans[,1]; ty <- tans[,2]; tz <- tans[,3 ]} else {tx <- tans[,1]; ty <- tans[,2]}
  if(k==3){
    RP = tx*ycx+ty*ycy+tz*ycz
    MP = rowSums(tans^2)
    RMP = RP/MP; RMP[!is.finite(RMP)] = 0
    y - tans*RMP
  } else {
    RP = tx*ycx+ty*ycy 
    MP = rowSums(tans^2)
    RMP = RP/MP; RMP[!is.finite(RMP)] = 0
    y - tans*RMP
  }
}

# semilandmarks.slide.surf.procD
# slides landmarks within PC planes tangent to surfaces using minimized ProcD
# used in pGpa.wSliders

# Based on the equation y -U%*%solve(crossprod(U))%*%crossprod(U,(y-ref))
#                          LP             MP                RP
# left part (LP), middle part(MP) and right part (RP) are accomplished faster in the code in this function.

semilandmarks.slide.surf.procD <- function(y,surf, ref){
  yc <- y - ref
  p <- nrow(yc); k <-ncol(yc)
  PC <- getSurfPCs(y, surf)
  p1x <- PC$p1x; p1y <- PC$p1y; p1z <- PC$p1z; p2x <- PC$p2x; p2y <- PC$p2y; p2z <- PC$p2z
  if(k==3) {ycx <- yc[,1]; ycy <- yc[,2]; ycz <- yc[,3 ]} else {ycx <- yc[,1]; ycy <- yc[,2]}
  if(k==3){
    RP = c(p1x*ycx+p1y*ycy+p1z*ycz,p2x*ycx+p2y*ycy+p2z*ycz)
    MP = c(p1x^2+p1y^2+p1z^2,p2x^2+p2y^2+p2z^2)
    RMP = RP/MP; RMP[!is.finite(RMP)] = 0
    y - cbind(p1x*RMP[1:p]+p2x*RMP[-(1:p)], p1y*RMP[1:p]+p2y*RMP[-(1:p)], 
              p1z*RMP[1:p]+p2z*RMP[-(1:p)])
  } else {
    RP = c(p1x*ycx+p1y*ycy, p2x*ycx+p2y*ycy)
    MP = c(p1x^2+p1y^2, p2x^2+p2y^2)
    RMP = RP/MP; RMP[!is.finite(RMP)] = 0
    y - cbind(p1x*RMP[1:p]+p2x*RMP[-(1:p)], p1y*RMP[1:p]+p2y*RMP[-(1:p)])
  }
}

# semilandmarks.slide.tangents.surf.procD
# slides landmarks along tangents of curves and within tangent planes on surfaces using minimized ProcD
# used in pGpa.wSliders

# Based on the equation y -U%*%solve(crossprod(U))%*%crossprod(U,(y-ref))
#                          LP             MP                RP
# left part (LP), middle part(MP) and right part (RP) are accomplished faster in the code in this function.

semilandmarks.slide.tangents.surf.procD <- function(y,tans,surf, ref){
  yc <- y - ref
  p <- nrow(yc); k <-ncol(yc)
  if(k==3) {tx <- tans[,1]; ty <- tans[,2]; tz <- tans[,3 ]} else {tx <- tans[,1]; ty <- tans[,2]}
  PC <- getSurfPCs(y, surf)
  p1x <- PC$p1x; p1y <- PC$p1y; p1z <- PC$p1z; p2x <- PC$p2x; p2y <- PC$p2y; p2z <- PC$p2z
  if(k==3) {ycx <- yc[,1]; ycy <- yc[,2]; ycz <- yc[,3 ]} else {ycx <- yc[,1]; ycy <- yc[,2]}
  if(k==3){
    RPt = tx*ycx+ty*ycy+tz*ycz
    MPt = rowSums(tans^2)
    RMPt = RPt/MPt; RMPt[!is.finite(RMPt)] = 0
    RPs = c(p1x*ycx+p1y*ycy+p1z*ycz,p2x*ycx+p2y*ycy+p2z*ycz)
    MPs = c(p1x^2+p1y^2+p1z^2,p2x^2+p2y^2+p2z^2)
    RMPs = RPs/MPs; RMPs[!is.finite(RMPs)] = 0
    y - (tans*RMPt+cbind(p1x*RMPs[1:p]+p2x*RMPs[-(1:p)], p1y*RMPs[1:p]+p2y*RMPs[-(1:p)], 
                         p1z*RMPs[1:p]+p2z*RMPs[-(1:p)]))
  } else {
    RPt = tx*ycx+ty*ycy
    MPt = rowSums(tans^2)
    RMPt = RPt/MPt; RMPt[!is.finite(RMPt)] = 0
    RPs = c(p1x*ycx+p1y*ycy,p2x*ycx+p2y*ycy)
    MPs = c(p1x^2+p1y^2,p2x^2+p2y^2)
    RMPs = RPs/MPs; RMPs[!is.finite(RMPs)] = 0
    y - (tans*RMPt+cbind(p1x*RMPs[1:p]+p2x*RMPs[-(1:p)], p1y*RMPs[1:p]+p2y*RMPs[-(1:p)]))
  }
}

# BE.slide
# performs sliding iterations using bending energy
# used in pGpa.wSliders
BE.slide <- function(curves, surf, Ya, ref, max.iter=5){# see pGpa.wCurves for variable meaning
  n <- length(Ya); p <- nrow(Ya[[1]]); k <- ncol(Ya[[1]])
  iter <- 1 # from initial rotation of Ya
  pb <- txtProgressBar(min = 0, max = max.iter, initial = 0, style=3) 
  slid0 <- Ya
  Q <- ss0 <- sum(Reduce("+",Ya)^2)/n
  setTxtProgressBar(pb,iter)
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
    ss <- sum(Reduce("+",slid)^2)/n
    slid0 <- apply.pPsup(ref,slid)
    ref = cs.scale(Reduce("+", slid0)/n)
    Q <- abs(ss0-ss)
    ss0 <- ss
    setTxtProgressBar(pb,iter)
    if(iter >= max.iter) break
  }
  if(iter < max.iter) setTxtProgressBar(pb,max.iter)
  close(pb)
  list(coords=slid0, consensus=ref, iter=iter+1, Q=Q)
}

# .BE.slide
# same as BE.slide, but without progress bar option
# used in pGpa.wSliders
.BE.slide <- function(curves, surf, Ya, ref, max.iter=5){# see pGpa.wCurves for variable meaning
  n <- length(Ya); p <- nrow(Ya[[1]]); k <- ncol(Ya[[1]])
  iter <- 1 # from initial rotation of Ya
  slid0 <- Ya
  Q <- ss0 <- sum(Reduce("+",Ya)^2)/n
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
    ss <- sum(Reduce("+",slid)^2)/n
    slid0 <- apply.pPsup(ref,slid)
    ref = cs.scale(Reduce("+", slid0)/n)
    Q <- abs(ss0-ss)
    ss0 <- ss
    if(iter >= max.iter) break
  }
  list(coords=slid0, consensus=ref, iter=iter+1, Q=Q)
}

# procD.slide
# performs sliding iterations using minimized ProcD
# used in pGpa.wSliders
procD.slide <- function(curves, surf, Ya, ref, max.iter=5){# see pGpa.wCurves for variable meaning
  n <- length(Ya); p <- nrow(Ya[[1]]); k <- ncol(Ya[[1]])
  iter <- 1 # from initial rotation of Ya
  pb <- txtProgressBar(min = 0, max = max.iter, initial = 0, style=3) 
  slid0 <- Ya
  Q <- ss0 <- sum(Reduce("+",Ya)^2)/n
  setTxtProgressBar(pb,iter)
  while(Q > 0.0001){
    iter <- iter+1
    if(!is.null(curves)) tans <- Map(function(y) tangents(curves, y, scaled=TRUE), slid0)
    if(is.null(surf) & !is.null(curves))
      slid <- Map(function(tn,y) semilandmarks.slide.tangents.procD(y, tn, ref), tans, slid0)
    if(!is.null(surf) & is.null(curves))
      slid <- Map(function(y) semilandmarks.slide.surf.procD(y, surf, ref), slid0)
    if(!is.null(surf) & !is.null(curves))
      slid <- Map(function(tn,y) semilandmarks.slide.tangents.surf.procD(y, tn, surf, ref), tans, slid0)
    ss <- sum(Reduce("+",slid)^2)/n
    slid0 <- apply.pPsup(ref,slid)
    ref = cs.scale(Reduce("+", slid0)/n)
    Q <- abs(ss0-ss)
    ss0 <- ss
    setTxtProgressBar(pb,iter)
    if(iter >=max.iter) break
  }
  if(iter < max.iter) setTxtProgressBar(pb,max.iter)
  close(pb)
  list(coords=slid0, consensus=ref, iter=iter+1, Q=Q)
}

# .procD.slide
# same as procD.slide, but without progress bar option
# used in pGpa.wSliders
.procD.slide <- function(curves, surf, Ya, ref, max.iter=5){# see pGpa.wCurves for variable meaning
  n <- length(Ya); p <- nrow(Ya[[1]]); k <- ncol(Ya[[1]])
  iter <- 1 # from initial rotation of Ya
  slid0 <- Ya
  Q <- ss0 <- sum(Reduce("+",Ya)^2)/n
  while(Q > 0.0001){
    iter <- iter+1
    if(!is.null(curves)) tans <- Map(function(y) tangents(curves, y, scaled=TRUE), slid0)
    if(is.null(surf) & !is.null(curves))
      slid <- Map(function(tn,y) semilandmarks.slide.tangents.procD(y, tn, ref), tans, slid0)
    if(!is.null(surf) & is.null(curves))
      slid <- Map(function(y) semilandmarks.slide.surf.procD(y, surf, ref), slid0)
    if(!is.null(surf) & !is.null(curves))
      slid <- Map(function(tn,y) semilandmarks.slide.tangents.surf.procD(y, tn, surf, ref), tans, slid0)
    ss <- sum(Reduce("+",slid)^2)/n
    slid0 <- apply.pPsup(ref,slid)
    ref = cs.scale(Reduce("+", slid0)/n)
    Q <- abs(ss0-ss)
    ss0 <- ss
    if(iter >=max.iter) break
  }
  list(coords=slid0, consensus=ref, iter=iter+1, Q=Q)
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
  list(coords= Ya, CS=CS, iter=iter, consensus=M, Q=Q, nsliders=NULL)
}

# .pGPA.wSliders
# same as pGPA.wSliders, without option for progress bar
# used in gpagen
.pGpa.wSliders <- function(Y, curves, surf, ProcD = TRUE, PrinAxes = FALSE, Proj = FALSE, max.iter = 5){
  n <- length(Y); p <- nrow(Y[[1]]); k <- ncol(Y[[1]])
  Yc <- Map(function(y) center.scale(y), Y)
  CS <- sapply(Yc,"[[","CS")
  Ya <- lapply(Yc,"[[","coords")
  Ya <- apply.pPsup(Ya[[1]], Ya)
  M <- Reduce("+", Ya)/n
  if(ProcD == FALSE) gpa.slide <- .BE.slide(curves, surf, Ya, ref=M, max.iter=max.iter) else 
    gpa.slide <- .procD.slide(curves, surf, Ya, ref=M, max.iter=max.iter)
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
  list(coords= Ya, CS=CS, iter=iter, consensus=M, Q=Q, nsliders=NULL)
}

# tps
#
#
tps<-function(matr, matt, n,sz=1.5, pt.bg="black",
              grid.col="black", grid.lwd=1, grid.lty=1, refpts=FALSE){		#DCA: altered from J. Claude: 2D only	
  xm<-min(matr[,1])
  ym<-min(matr[,2])
  xM<-max(matr[,1])
  yM<-max(matr[,2])
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

# tps2d
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
coefx<-fast.solve(L)%*%m2[,1]
coefy<-fast.solve(L)%*%m2[,2]
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

# tps2d3d
#
#
tps2d3d<-function(M, matr, matt, PB=TRUE){		#DCA: altered from J. Claude 2008  
  p<-dim(matr)[1]; k<-dim(matr)[2];q<-dim(M)[1]
  Pdist<-as.matrix(dist(matr))
  ifelse(k==2,P<-Pdist^2*log(Pdist^2),P<- Pdist) 
  P[which(is.na(P))]<-0
  Q<-cbind(1, matr)
  L<-rbind(cbind(P,Q), cbind(t(Q),matrix(0,k+1,k+1)))
  m2<-rbind(matt, matrix(0, k+1, k))   
  coefx<-fast.solve(L)%*%m2[,1]
  coefy<-fast.solve(L)%*%m2[,2]
  if(k==3){coefz<-fast.solve(L)%*%m2[,3]}
  fx<-function(matr, M, coef, step){
    Xn<-numeric(q)
    for (i in 1:q){
      Z<-apply((matr-matrix(M[i,],p,k,byrow=TRUE))^2,1,sum)  
      ifelse(k==2,Z1<-Z*log(Z),Z1<-sqrt(Z)); Z1[which(is.na(Z1))]<-0
      ifelse(k==2,Xn[i]<-coef[p+1]+coef[p+2]*M[i,1]+coef[p+3]*M[i,2]+sum(coef[1:p]*Z1),
             Xn[i]<-coef[p+1]+coef[p+2]*M[i,1]+coef[p+3]*M[i,2]+coef[p+4]*M[i,3]+sum(coef[1:p]*Z1))
      if(PB==TRUE){setTxtProgressBar(pb, step + i)}
    }    
    Xn}
  matg<-matrix(NA, q, k)
  if(PB==TRUE){pb <- txtProgressBar(min = 0, max = q*k, style = 3) }
  matg[,1]<-fx(matr, M, coefx, step = 1)
  matg[,2]<-fx(matr, M, coefy, step=q)
  if(k==3){matg[,3]<-fx(matr, M, coefz, step=q*2)
  } 
  if(PB==TRUE) close(pb)
  matg
}

# pcoa
# acquires principal coordinates from distance matrices
# used in all functions with 'procD.lm" via procD.fit
pcoa <- function(D){
  options(warn=-1)
  if(class(D) != "dist") stop("function only works with distance matrices")
  cmd <- cmdscale(D, k=attr(D, "Size") -1, eig=TRUE)
  options(warn=0)
  p <- length(cmd$eig[zapsmall(cmd$eig) > 0])
  Yp <- cmd$points[,1:p]
  Yp
}

# In development for various functions

model.matrix.g <- function(f1, data = NULL) {
  f1 <- as.formula(f1)
  Terms <- terms(f1)
  labs <- attr(Terms, "term.labels")
  if(!is.null(data)) {
    matches <- na.omit(match(labs, names(data)))
    dat <- as.data.frame(data[matches]) 
  } else dat <- NULL
  model.matrix(f1, data=dat)
}

# procD.fit
# lm-like fit modified for Procrustes residuals 
# general workhorse for all 'procD.lm' functions
# used in all 'procD.lm' functions
procD.fit <- function(f1, keep.order=FALSE, pca=TRUE, data=NULL, 
                      SS.type = c("I", "II", "III"), ...){
  form.in <- formula(f1)
  if(any(class(f1)=="lm")) {
    weights <- f1$weights
    contrasts <- f1$contrasts
    offset <-f1$offset
    data <- dat <- model.frame(f1)
    Y <- as.matrix(data[1])
  } else {
    if(!is.null(data)) {
      dat <- data
      Y <- eval(form.in[[2]], envir = data)
    } else {
      Y <- eval(form.in[[2]], parent.frame())
      dat <- model.frame(form.in[-2])
    }
    dat <- droplevels(dat)
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
  if(!is.null(SS.type)) SS.type <- match.arg(SS.type) else SS.type <- "I"
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
  if(form.terms[[2]] == 1) {
    X <- model.matrix(~1, data = data.frame(x=rep(1,n)))
    Terms <- terms(form.terms, keep.order=keep.order, data = data.frame(x=rep(1,n)))
  } else {
    df.check <- c("matrix", "vector", "numeric", "integer", "character", "logical", "factor")
    df.check <- sapply(1:length(dat), function(j) any(class(dat[[j]]) == df.check))
    mf <- as.data.frame(dat[df.check])
    Terms <- terms(form.terms, keep.order=keep.order, data = mf)
    if(nrow(Y) != nrow(mf)) stop("Different numbers of specimens in dependent and independent variables")
    X <- model.matrix(Terms, data=mf)
  }

  if(is.null(weights) & !is.null(dots$weights)) w <- dots$weights else 
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
  if(SS.type == "III"){
    Xs <- lapply(1:length(uk), function(j)  Xj <- X[, X.k %in% uk[-j]])
    Xs <- c(Xs[-1], list(X))
  } else if(SS.type == "II") {
    factors <- attr(Terms, "factors")
    fac.guide <- colSums(factors)
    Xs <- lapply(1:length(fac.guide), function(j){
      x <- fac.guide[j]
      if(x == 1){
        fc <- which(fac.guide == 1)
        model.matrix(Terms[fc[fc != j]])
      } else {
        model.matrix(Terms[fac.guide < x])
      }
    })
    Xs <- c(Xs, list(X))
  } else
      Xs <- lapply(1:length(uk), function(j)  Xj <- X[, X.k %in% uk[1:j]])
  QRs <- lapply(Xs, function(x) qr(x))
  fitted <- lapply(QRs, function(x) qr.fitted(x,Y))
  residuals <- lapply(QRs, function(x) qr.resid(x,Y))
  coefficients <- lapply(QRs, function(x) qr.coef(x,Y))
  if(sum(w)==n){
    wXs <- Xs; wQRs <- QRs; wFitted <- fitted; wResiduals <- residuals
    wCoefficients <- coefficients
  } else {
    wXs <- lapply(Xs, function(x) x*sqrt(w))
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
  mf.out <- droplevels(data.frame(Y=Y, mf.out))
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

perm.index <-function(n, iter, seed=NULL){
  if(is.null(seed)) seed = iter else
    if(seed == "random") seed = sample(1:iter,1) else
      if(!is.numeric(seed)) seed = iter
      set.seed(seed)
      ind <- c(list(1:n),(Map(function(x) sample.int(n,n), 1:iter)))
      rm(.Random.seed, envir=globalenv())
      ind
}

# perm.CR.index
# creates a permutation index for resampling, shuffling landmarks
# used in all functions utilizing CR (modularity)

perm.CR.index <- function(g, k, iter, seed=NULL){ # g is numeric partititon.gp
  if(is.null(seed)) seed = iter else
    if(seed == "random") seed = sample(1:iter,1) else
      if(!is.numeric(seed)) seed = iter
      set.seed(seed)
      p <- length(g)
      ind <- c(list(1:p),(Map(function(x) sample.int(p,p), 1:iter)))
      ind <- Map(function(x) g[x], ind)
      ind <- Map(function(x) as.factor(rep(x,k,each = k, length=p*k)), ind)
      rm(.Random.seed, envir=globalenv())
      ind
}

# boot.index
# creates a bootstrap index for resampling
# used in modularity test functions
boot.index <-function(n, iter, seed=NULL){
  if(is.null(seed)) seed = iter else
    if(seed == "random") seed = sample(1:iter,1) else
      if(!is.numeric(seed)) seed = iter
      set.seed(seed)
      ind <- c(list(1:n),(Map(function(x) sample.int(n, n, replace = TRUE), 1:iter)))
      rm(.Random.seed, envir=globalenv())
      ind
}

# fastFit
# calculates fitted values for a linear model, after decomoposition of X to get U
# used in SS.iter and Fpgls.iter; future need: advanced.procD.lm 
fastFit <- function(U,y,n,p){
  if(!is.matrix(y)) y <- as.matrix(y)
  if(p > n) tcrossprod(U)%*%y else 
    U%*%crossprod(U,y) 
}

# fastLM
# calculates fitted values and residuals, after fastFit
# placeholder in case needed later
fastLM<- function(U,y){
  p <- dim(y)[2]; n <- dim(y)[1]
  yh <- fastFit(U,y,n,p)
  list(fitted = yh, residuals = y-yh)
}

# SS.iter
# calculates SS in random iterations of a resampling procedure
# used in nearly all 'procD.lm' functions, unless pgls in used
SS.iter = function(pfit,iter, seed = NULL, Yalt="RRPP", SS.type= NULL){
  Y <- as.matrix(pfit$Y)
  Xf <- as.matrix(pfit$X)
  k <- length(pfit$Xs)
  n <- dim(Y)[1]; p <- dim(Y)[2]
  Yh <- pfit$fitted
  E <- pfit$residuals
  w <- pfit$weights
  if(is.null(SS.type)) SS.type <- "I"
  if(is.na(match(SS.type, c("I","III")))) SS.type <- "I"
  if(SS.type == "III") {
    Ur <- lapply(pfit$wQRs[1:(k-1)], function(x) qr.Q(x))
    Uf <- qr.Q(qr(Xf))
    Uf <- lapply(1:(k-1), function(.) Uf)
  } else if(SS.type == "II") {
    Ur <- lapply(pfit$wQRs[1:(k-1)], function(x) qr.Q(x))
    Terms <- pfit$Terms
    factors <- attr(Terms, "factors")
    fac.guide <- colSums(factors)
    Xfs <- lapply(1:length(fac.guide), function(j){
      x <- fac.guide[j]
      if(x == 1){
        model.matrix(Terms[fac.guide==1])
      } else {
        fc <- c(which(fac.guide < x), j)
        model.matrix(Terms[fc])
      }
    })
    if(sum(w) == n) wXs <- Xs else wXs <- lapply(Xfs, function(x) x*sqrt(w))
    Uf <- lapply(wXs[2:k], function(x) qr.Q(qr(x)))
    } else {
      Ur <- lapply(pfit$wQRs[1:(k-1)], function(x) qr.Q(x))
      Uf <- lapply(pfit$wQRs[2:k], function(x) qr.Q(x))
      }
  q <- qr(Xf)$rank
  if(q/n > 0.4) Pr <- Map(function(uf,ur) tcrossprod(uf) - tcrossprod(ur), Uf, Ur) else
    Pr <- NULL
  ind = perm.index(n,iter, seed=seed)
  SS <- NULL
  pb <- txtProgressBar(min = 0, max = ceiling(iter/100), initial = 0, style=3) 
  jj <- iter+1
  step <- 1
  if(jj > 100) j <- 1:100 else j <- 1:jj
  while(jj > 0){
    ind.j <- ind[j]
    if(Yalt=="RRPP") {
      if(sum(w)==n) {
        Yr = Map(function(x) (Map(function(y,e) e[x,]+y, Yh[1:(k-1)], E[1:(k-1)])),ind.j)
      } else {
        Yr = Map(function(x) (Map(function(y,e) (e[x,]+y)*sqrt(w), Yh[1:(k-1)], E[1:(k-1)])),ind.j) 
      }} else {
        if(sum(w)==n) {
          Yr = Map(function(x) Map(function(y) y[x,], lapply(1:(k-1),function(.) Y)),ind.j)
        } else {
          Yr = Map(function(x) Map(function(y) (y[x,])*sqrt(w), lapply(1:(k-1),function(.) Y)),ind.j)
        }
      }
    if(!is.null(Pr)) {
      SS.temp <- lapply(1:length(j), function(j){ 
        mapply(function(p,y) sum((p%*%y)^2), 
               Pr,Yr[[j]])})
    } else
    SS.temp <- lapply(1:length(j), function(j){ 
      mapply(function(ur,uf,y) sum((fastFit(uf,y,n,p) - fastFit(ur,y,n,p))^2), 
             Ur, Uf,Yr[[j]])})
    SS <- c(SS, SS.temp)
    jj <- jj-length(j)
    if(jj > 100) kk <- 1:100 else kk <- 1:jj
    j <- j[length(j)] +kk
    setTxtProgressBar(pb,step)
    step <- step+1
  }
  close(pb)
  simplify2array(SS)
}

# .SS.iter
# same as SS.iter, but without progress bar option
# used in nearly all 'procD.lm' functions, unless pgls in used
.SS.iter = function(pfit,iter, seed = NULL, Yalt="RRPP", SS.type=NULL){
  Y <- as.matrix(pfit$Y)
  Xf <- as.matrix(pfit$X)
  k <- length(pfit$Xs)
  n <- dim(Y)[1]; p <- dim(Y)[2]
  Yh <- pfit$fitted
  E <- pfit$residuals
  w<- pfit$weights
  if(is.null(SS.type)) SS.type <- "I"
  if(is.na(match(SS.type, c("I","III")))) SS.type <- "I"
  if(SS.type == "III") {
    Ur <- lapply(pfit$wQRs[1:(k-1)], function(x) qr.Q(x))
    Uf <- qr.Q(qr(Xf))
    Uf <- lapply(1:(k-1), function(.) Uf)
  } else if(SS.type == "II") {
    Ur <- lapply(pfit$wQRs[1:(k-1)], function(x) qr.Q(x))
    Terms <- pfit$Terms
    factors <- attr(Terms, "factors")
    fac.guide <- colSums(factors)
    Xfs <- lapply(1:length(fac.guide), function(j){
      x <- fac.guide[j]
      if(x == 1){
        model.matrix(Terms[fac.guide==1])
      } else {
        fc <- c(which(fac.guide < x), j)
        model.matrix(Terms[fc])
      }
    })
    if(sum(w) == n) wXs <- Xs else wXs <- lapply(Xfs, function(x) x*sqrt(w))
    Uf <- lapply(wXs[2:k], function(x) qr.Q(x))
  } else {
    Ur <- lapply(pfit$wQRs[1:(k-1)], function(x) qr.Q(x))
    Uf <- lapply(pfit$wQRs[2:k], function(x) qr.Q(x))
  }
  q <- qr(Xf)$rank
  if(q/n > 0.4) Pr <- Map(function(uf,ur) tcrossprod(uf) - tcrossprod(ur), Uf, Ur) else
    Pr <- NULL
  ind = perm.index(n,iter, seed=seed)
  SS <- NULL
  jj <- iter+1
  if(jj > 100) j <- 1:100 else j <- 1:jj
  while(jj > 0){
    ind.j <- ind[j]
    if(Yalt=="RRPP") {
      if(sum(w)==n) {
        Yr = Map(function(x) (Map(function(y,e) e[x,]+y, Yh[1:(k-1)], E[1:(k-1)])),ind.j)
      } else {
        Yr = Map(function(x) (Map(function(y,e) (e[x,]+y)*sqrt(w), Yh[1:(k-1)], E[1:(k-1)])),ind.j) 
      }} else {
        if(sum(w)==n) {
          Yr = Map(function(x) Map(function(y) y[x,], lapply(1:(k-1),function(.) Y)),ind.j)
        } else {
          Yr = Map(function(x) Map(function(y) (y[x,])*sqrt(w), lapply(1:(k-1),function(.) Y)),ind.j)
        }
      }
    
    if(!is.null(Pr)) {
      SS.temp <- lapply(1:length(j), function(j){ 
        mapply(function(p,y) sum((p%*%y)^2), 
               Pr,Yr[[j]])})
    } else
      SS.temp <- lapply(1:length(j), function(j){ 
        mapply(function(ur,uf,y) sum((fastFit(uf,y,n,p) - fastFit(ur,y,n,p))^2), 
               Ur, Uf,Yr[[j]])})
    SS <- c(SS, SS.temp)
    jj <- jj-length(j)
    if(jj > 100) kk <- 1:100 else kk <- 1:jj
    j <- j[length(j)] +kk
  }
  simplify2array(SS)
}

# SS.iter.bilat.symmetry
# calculates SS in random iterations of a resampling procedure
# used in nearly all 'procD.lm' functions, unless pgls in used
SS.iter.bilat.symmetry = function(pfit,iter, seed = NULL, Yalt="RRPP", SS.type= NULL){
  Y <- as.matrix(pfit$Y)
  Xf <- as.matrix(pfit$X)
  k <- length(pfit$Xs)
  n <- dim(Y)[1]; p <- dim(Y)[2]
  Yh <- pfit$fitted
  E <- pfit$residuals
  w<- pfit$weights
  if(is.null(SS.type)) SS.type <- "I"
  if(is.na(match(SS.type, c("I","III")))) SS.type <- "I"
  if(SS.type == "III") {
    Ur <- lapply(pfit$wQRs[1:(k-1)], function(x) qr.Q(x))
    Uf <- qr.Q(qr(Xf))
    Uf <- lapply(1:(k-1), function(.) Uf)
  } else {
    Ur <- lapply(pfit$wQRs[1:(k-1)], function(x) qr.Q(x))
    Uf <- lapply(pfit$wQRs[2:k], function(x) qr.Q(x))
  }
  Pr <- Map(function(uf,ur) tcrossprod(uf) - tcrossprod(ur), Uf, Ur)
  ind = perm.index(n,iter, seed=seed)
  SS <- NULL
  pb <- txtProgressBar(min = 0, max = ceiling(iter/100), initial = 0, style=3) 
  jj <- iter+1
  step <- 1
  if(jj > 100) j <- 1:100 else j <- 1:jj
  while(jj > 0){
    ind.j <- ind[j]
    if(Yalt=="RRPP") {
      if(sum(w)==n) {
        Yr = Map(function(x) (Map(function(y,e) e[x,]+y, Yh[1:(k-1)], E[1:(k-1)])),ind.j)
      } else {
        Yr = Map(function(x) (Map(function(y,e) (e[x,]+y)*sqrt(w), Yh[1:(k-1)], E[1:(k-1)])),ind.j) 
      }} else {
        if(sum(w)==n) {
          Yr = Map(function(x) Map(function(y) y[x,], lapply(1:(k-1),function(.) Y)),ind.j)
        } else {
          Yr = Map(function(x) Map(function(y) (y[x,])*sqrt(w), lapply(1:(k-1),function(.) Y)),ind.j)
        }
      }
    SS.temp <- lapply(1:length(j), function(j){ 
      mapply(function(p,y) sum((p%*%y)^2), 
             Pr,Yr[[j]])})
    SS <- c(SS, SS.temp)
    jj <- jj-length(j)
    if(jj > 100) kk <- 1:100 else kk <- 1:jj
    j <- j[length(j)] +kk
    setTxtProgressBar(pb,step)
    step <- step+1
  }
  close(pb)
  simplify2array(SS)
}

# .SS.iter.bilat.symmetry
# same as SS.iter, but without progress bar option
# used in nearly all 'procD.lm' functions, unless pgls in used
.SS.iter.bilat.symmetry = function(pfit,iter, seed = NULL, Yalt="RRPP", SS.type=NULL){
  Y <- as.matrix(pfit$Y)
  Xf <- as.matrix(pfit$X)
  k <- length(pfit$Xs)
  n <- dim(Y)[1]; p <- dim(Y)[2]
  Yh <- pfit$fitted
  E <- pfit$residuals
  w<- pfit$weights
  if(is.null(SS.type)) SS.type <- "I"
  if(is.na(match(SS.type, c("I","III")))) SS.type <- "I"
  if(SS.type == "III") {
    Ur <- lapply(pfit$wQRs[1:(k-1)], function(x) qr.Q(x))
    Uf <- qr.Q(qr(Xf))
    Uf <- lapply(1:(k-1), function(.) Uf)
  } else {
    Ur <- lapply(pfit$wQRs[1:(k-1)], function(x) qr.Q(x))
    Uf <- lapply(pfit$wQRs[2:k], function(x) qr.Q(x))
  }
  Pr <- Map(function(uf,ur) tcrossprod(uf) - tcrossprod(ur), Uf, Ur)
  ind = perm.index(n,iter, seed=seed)
  SS <- NULL
  jj <- iter+1
  if(jj > 100) j <- 1:100 else j <- 1:jj
  while(jj > 0){
    ind.j <- ind[j]
    if(Yalt=="RRPP") {
      if(sum(w)==n) {
        Yr = Map(function(x) (Map(function(y,e) e[x,]+y, Yh[1:(k-1)], E[1:(k-1)])),ind.j)
      } else {
        Yr = Map(function(x) (Map(function(y,e) (e[x,]+y)*sqrt(w), Yh[1:(k-1)], E[1:(k-1)])),ind.j) 
      }} else {
        if(sum(w)==n) {
          Yr = Map(function(x) Map(function(y) y[x,], lapply(1:(k-1),function(.) Y)),ind.j)
        } else {
          Yr = Map(function(x) Map(function(y) (y[x,])*sqrt(w), lapply(1:(k-1),function(.) Y)),ind.j)
        }
      }
    
    SS.temp <- lapply(1:length(j), function(j){ 
      mapply(function(p,y) sum((p%*%y)^2), 
             Pr,Yr[[j]])})
    SS <- c(SS, SS.temp)
    jj <- jj-length(j)
    if(jj > 100) kk <- 1:100 else kk <- 1:jj
    j <- j[length(j)] +kk
  }
  simplify2array(SS)
}

# SS.pgls.iter
# calculates F values in random iterations of a resampling procedure, with pgls involved
# used in the 'procD.lm' functions where pgls is used
# Formerly Fpgls.iter
SS.pgls.iter = function(pfit,Pcor,iter, seed=NULL, Yalt="RRPP"){
  Y <- as.matrix(pfit$Y)
  k <- length(pfit$QRs)
  n <- dim(Y)[1]; p <- dim(Y)[2]
  Yh <- pfit$fitted
  E <- pfit$residuals
  w<- pfit$weights
  wQRs <- pfit$wQRs
  q <- sapply(1:k, function(j) wQRs[[j]]$rank)
  df <- q[-1] - q[1:(k-1)]
  q <- q[length(q)]
  if(any(rownames(Pcor) != rownames(Y))) Pcor <- Pcor[rownames(Y),rownames(Y)]
  PwXs <- lapply(pfit$wXs, function(x) crossprod(Pcor,as.matrix(x)))
  Xr <- lapply(PwXs[1:(k-1)], function(x) as.matrix(x))
  Xf <- lapply(PwXs[2:k], function(x) as.matrix(x))
  Ur <- lapply(Xr, function(x) qr.Q(qr(x)))
  Uf <- lapply(Xf, function(x) qr.Q(qr(x)))
  if(p > n) Ptransr <- lapply(Ur, function(x) tcrossprod(x)%*%Pcor) else
    Ptransr <- lapply(Ur, function(x) x%*%crossprod(x,Pcor)) 
  if(p > n) Ptransf <- lapply(Uf, function(x) tcrossprod(x)%*%Pcor) else
    Ptransf <- lapply(Uf, function(x) x%*%crossprod(x,Pcor)) 
  Ptrans <- Map(function(r,f) f-r, Ptransr, Ptransf)
  ind = perm.index(n,iter, seed=seed)
  SS <- SSEs <- NULL
  pb <- txtProgressBar(min = 0, max = ceiling(iter/100), initial = 0, style=3) 
  jj <- iter+1
  step <- 1
  if(jj > 100) j <- 1:100 else j <- 1:jj
  while(jj > 0){
    ind.j <- ind[j]
    if(Yalt=="RRPP") {
      if(sum(w)==n) {
        Yr = Map(function(x) (Map(function(y,e) e[x,]+y, Yh[1:(k-1)], E[1:(k-1)])),ind.j)
      } else {
        Yr = Map(function(x) (Map(function(y,e) (e[x,]+y)*sqrt(w), Yh[1:(k-1)], E[1:(k-1)])),ind.j) 
      }} else {
        if(sum(w)==n) {
          Yr = Map(function(x) Map(function(y) y[x,], lapply(1:(k-1),function(.) Y)),ind.j)
        } else {
          Yr = Map(function(x) Map(function(y) (y[x,])*sqrt(w), lapply(1:(k-1),function(.) Y)),ind.j)
        }
      }
    SS.temp <- lapply(1:length(j), function(j){ 
      mapply(function(p,y) sum((p%*%y)^2), Ptrans,Yr[[j]])})
    SSEs.temp <- Map(function(y) sum((Pcor%*%y[[k-1]]-Ptransf[[k-1]]%*%y[[k-1]])^2), Yr)
    SS <- c(SS,SS.temp)
    SSEs <- c(SSEs,SSEs.temp)
    jj <- jj-length(j)
    if(jj > 100) kk <- 1:100 else kk <- 1:jj
    j <- j[length(j)] +kk
    setTxtProgressBar(pb,step)
    step <- step+1
  }
  close(pb)
  SS <- as.matrix(rbind(simplify2array(SS), unlist(SSEs)))
  SS
}

# .SS.pgls.iter
# same as Fpgls.iter, but without progress bar option
# used in the 'procD.lm' functions where pgls is used
# Formerly .Fpgls.iter
.SS.pgls.iter = function(pfit,Pcor,iter, seed=NULL, Yalt="RRPP"){
  Y <- as.matrix(pfit$Y)
  k <- length(pfit$QRs)
  n <- dim(Y)[1]; p <- dim(Y)[2]
  Yh <- pfit$fitted
  E <- pfit$residuals
  w<- pfit$weights
  wQRs <- pfit$wQRs
  q <- sapply(1:k, function(j) wQRs[[j]]$rank)
  df <- q[-1] - q[1:(k-1)]
  q <- q[length(q)]
  if(any(rownames(Pcor) != rownames(Y))) Pcor <- Pcor[rownames(Y),rownames(Y)]  
  PwXs <- lapply(pfit$wXs, function(x) crossprod(Pcor,as.matrix(x)))
  Xr <- lapply(PwXs[1:(k-1)], function(x) as.matrix(x))
  Xf <- lapply(PwXs[2:k], function(x) as.matrix(x))
  Ur <- lapply(Xr, function(x) qr.Q(qr(x)))
  Uf <- lapply(Xf, function(x) qr.Q(qr(x)))
  if(p > n) Ptransr <- lapply(Ur, function(x) tcrossprod(x)%*%Pcor) else
    Ptransr <- lapply(Ur, function(x) x%*%crossprod(x,Pcor)) 
  if(p > n) Ptransf <- lapply(Uf, function(x) tcrossprod(x)%*%Pcor) else
    Ptransf <- lapply(Uf, function(x) x%*%crossprod(x,Pcor)) 
  Ptrans <- Map(function(r,f) f-r, Ptransr, Ptransf)
  ind = perm.index(n,iter, seed=seed)
  SS <- SSEs <- NULL
  jj <- iter+1
  if(jj > 100) j <- 1:100 else j <- 1:jj
  while(jj > 0){
    ind.j <- ind[j]
    if(Yalt=="RRPP") {
      if(sum(w)==n) {
        Yr = Map(function(x) (Map(function(y,e) e[x,]+y, Yh[1:(k-1)], E[1:(k-1)])),ind.j)
      } else {
        Yr = Map(function(x) (Map(function(y,e) (e[x,]+y)*sqrt(w), Yh[1:(k-1)], E[1:(k-1)])),ind.j) 
      }} else {
        if(sum(w)==n) {
          Yr = Map(function(x) Map(function(y) y[x,], lapply(1:(k-1),function(.) Y)),ind.j)
        } else {
          Yr = Map(function(x) Map(function(y) (y[x,])*sqrt(w), lapply(1:(k-1),function(.) Y)),ind.j)
        }
      }
    SS.temp <- lapply(1:length(j), function(j){ 
      mapply(function(p,y) sum((p%*%y)^2), Ptrans,Yr[[j]])})
    SSEs.temp <- Map(function(y) sum(fastLM(Uf[[k-1]],Pcor%*%y[[k-1]])$residuals^2), Yr)
    SS <- c(SS,SS.temp)
    SSEs <- c(SSEs,SSEs.temp)
    jj <- jj-length(j)
    if(jj > 100) kk <- 1:100 else kk <- 1:jj
    j <- j[length(j)] +kk
  }
  SS <- as.matrix(rbind(simplify2array(SS), unlist(SSEs)))
  SS
}

# anova.parts
# makes an ANOVA table based on SS from SS.iter
# used in nearly all 'procD.lm' functions
anova.parts <- function(pfit, SS, SS.type = NULL){ # SS from SS.iter
  Y <- pfit$wY
  k <- length(pfit$term.labels)
  if(is.matrix(SS) && NROW(SS) > k) SS <- SS[1:k,]
  if(is.null(SS.type)) SS.type <- "I"
  if(is.na(match(SS.type, c("I","III")))) SS.type <- "I"
  dfe <-sapply(pfit$wQRs, function(x) x$rank)
  if(SS.type == "III") df <- dfe[length(dfe)] - dfe[1:(length(dfe)-1)] else
    df <- dfe[2:(k+1)] - dfe[1:k]
  if(k==1) SS <- SS[1] else SS <- SS[,1]
  anova.terms <- pfit$term.labels
  SSY <- sum(qr.resid(qr(pfit$wX[,1]),pfit$wY)^2)
  MS <- SS/df
  R2 <- SS/SSY
  SSE.model <- sum(qr.resid(qr(pfit$wXs[[k+1]]),pfit$wY)^2)
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

# anova.parts.pgls
# makes an ANOVA table based on SS from Fpgls.iter
# used in nearly only in procD.pgls
anova.parts.pgls <- function(pfit, SS, SS.type = NULL){
  Y <- pfit$wY
  k <- length(pfit$term.labels)
  SSE <- SS[nrow(SS),][1]
  SS <- SS[1:k,]
  if(is.matrix(SS)) SS <- SS[,1] else SS <- SS[1]
  if(is.null(SS.type)) SS.type <- "I"
  if(is.na(match(SS.type, c("I","III")))) SS.type <- "I"
  dfe <-sapply(pfit$wQRs, function(x) x$rank)
  if(SS.type == "III") df <- dfe[length(dfe)] - dfe[1:(length(dfe)-1)] else
    df <- dfe[2:(k+1)] - dfe[1:k]
  anova.terms <- pfit$term.labels
  SSY <- sum(SS)+SSE
  MS <- SS/df
  R2 <- SS/SSY
  dfE <- nrow(Y)-(sum(df)+1)
  MSE <- SSE/dfE
  Fs <- MS/MSE
  df <- c(df,dfE,nrow(Y)-1)
  SS <- c(SS,SSE, SSY)
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
  anova.terms <- pfit$term.labels
  MS <- SS/df
  SSY <- sum(qr.resid(qr(pfit$wX[,1]),pfit$wY)^2)
  if(object.sym==TRUE){
    SS<-SS/2;SSY<-SSY/2
  }
  MS <- SS/df
  R2 <- SS/SSY
  F1 <- MS[1,]/MS[3,]
  F2 <- MS[2,]/MS[3,]
  if(k==4) F3 <- MS[3,]/MS[4,] else F3 <- NULL
  Fs <- rbind(F1,F2,F3, NA)
  rownames(Fs) <- rownames(MS) <- rownames(SS) <- rownames(R2) <- anova.terms
  Ps <- apply(Fs,1,pval)
  if(k==4) Ps[4] <- NA else Ps[3] <-NA
  Zs <- apply(Fs,1,effect.size)
  if(k==4) Zs[4] <- NA else Zs[3] <-NA
  anova.tab <- data.frame(df = df,SS = SS[,1],MS = MS[,1],Rsq=R2[,1],
                          F=Fs[,1], Z=Zs, Pr = Ps)
  rownames(anova.tab) <- anova.terms
  out <- list(anova.table = anova.tab, anova.terms = anova.terms,
              SS = SS[,1], df = df, R2 = R2[1:k,1], F = Fs[1:(k-1),1],
              random.Fs = Fs)
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

# effect.size
# Effect sizes (standard deviates) form random outcomes
# any analytical function
effect.size <- function(x, center = TRUE) {
  z = scale(x, center=center)
  n <- length(z)
  z[1]*sqrt((n-1)/(n))
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

# single.factor
# converts factorial designs to single-factor variables
# advanced.procD.lm

# helpers for single.factor
leveler <- function(x){ # x = data.frame of 2 columns
  a <-levels(x[,1])
  b <- levels(x[,2])
  na <- length(a); nb <- length(b)
  res <- NULL
  for(i in 1:na){
    ab <- c(a[i],b)
    res <- c(res, combn(ab, 2, simplify=FALSE)[1:nb])
  }
  res <- lapply(res, paste, collapse=":")
  simplify2array(res)
}

multileveler <- function(x){ # x = data.frame of 2 or more columns
  if(NCOL(x) == 0) y <- levels(x)
  if(NCOL(x) == 1) y <- levels(x[,1])
  if(NCOL(x) == 2) y <- leveler(x)
  if(NCOL(x) > 2){
    a <- leveler(x[,1:2])
    for(i in 3:NCOL(x)){
      b <- levels(x[,i])
      na <- length(a); nb <- length(b)
      res <- NULL
      for(i in 1:na){
        ab <- c(a[i],b)
        res <- c(res, combn(ab, 2, simplify=FALSE)[1:nb])
      }
      res <- lapply(res, paste, collapse=":")
      a <- simplify2array(res)
    }
    y <- a
  }
  return(y)
}

single.factor <- function(pfit) {# pfit = Procrustes fit
  Terms <- pfit$Terms
  dat <- pfit$data
  datClasses <- sapply(dat, function(x) data.class(x))
  facs <- dat[which(datClasses == "factor")]
  facs <- as.data.frame(facs)
  if(ncol(facs) > 1) fac <- factor(apply(facs, 1,function(x) paste(x, collapse=":"))) else 
    fac <- as.factor(unlist(facs))
  faclevels <- multileveler(facs)
  factor(fac, levels=faclevels)
}

# cov.extract
# Extracts covariates from design matrices
# advanced.procD.lm
cov.extract <- function(pfit) {
  Terms <- pfit$Terms
  mf <- model.frame(Terms, data = pfit$data)
  if(is.null(.getXlevels(Terms, mf))) covs <- NULL else
  {
    datClasses <- sapply(mf, function(x) data.class(x))
    datClasses <- which(datClasses == "numeric")
    covs <- mf[datClasses]
  }
  if(length(covs) == 0) covs <- NULL else covs <- as.matrix(covs)
  covs
}

# ls.means
# estimates ls.means from models with slopes and factors
# advanced.procD.lm
# THIS FUNCTION IS NO LONGER USED.  It is retained, however, for future debugging
# purposes, in the event updates have errors
ls.means = function(pfit, Y=NULL, g=NULL, data=NULL, Pcor = NULL) { 
  if(!is.null(Pcor) && is.null(Y)) stop("If Pcor is provided, Y cannot be null")
  if(is.null(Y)) Y <- pfit$wY 
  n <- nrow(Y)
  if(is.null(data)) dat <- pfit$data else dat <- data
  if(!is.null(g)) {
    if(is.data.frame(g)){
      fac.check <- sapply(g, is.factor)
      facs <- g[,fac.check]
    } else if(is.factor(g)) facs <- g else stop("groups input neither a factor nor factors")
    if(ncol(as.matrix(facs)) > 1) fac <- factor(apply(facs, 1,function(x) paste(x, collapse=":"))) else 
      fac <- as.factor(unlist(facs))
  } else fac <- single.factor(pfit)
  covs <- cov.extract(pfit) 
  if(is.null(covs)) X0 <- model.matrix(~fac) else {
    if(ncol(as.matrix(covs)) == 1) covs <- as.vector(covs)
    X0 <- model.matrix(~covs+fac)
  }
  X <- X0*sqrt(pfit$weights)
  if(!is.null(Pcor)) fit <- .lm.fit(Pcor%*%X,Y)$coefficients else fit <- .lm.fit(X,Y)$coefficients
  if(!is.null(covs)){
    Xcov.mean <- as.matrix(X[,1:ncol(model.matrix(~covs))])
    Xcov.mean <- sapply(1:ncol(Xcov.mean), 
                        function(j) rep(mean(Xcov.mean[,j]),nrow(Xcov.mean)))
    Xnew <- cbind(Xcov.mean, X[,-(1:ncol(Xcov.mean))])
  } else Xnew <- X
  lsm <- as.matrix(lm.fit(model.matrix(~fac+0),Xnew%*%fit)$coefficients)
  rownames(lsm) <- levels(fac)
  lsm
}

# apply.ls.means
# estimates ls.means from models with slopes and factors across random outcomes in a list
# via helper functions, quick.ls.means.setup, and quick.ls.means
# advanced.procD.lm
quick.ls.means.set.up <- function(pfit, g=NULL, data=NULL) {
  if(is.null(data)) dat <- pfit$data else dat <- data
  if(!is.null(g)) {
    if(is.data.frame(g)){
      fac.check <- sapply(g, is.factor)
      facs <- g[,fac.check]
    } else if(is.factor(g)) facs <- g else stop("groups input neither a factor nor factors")
    if(ncol(as.matrix(facs)) > 1) fac <- factor(apply(facs, 1,function(x) paste(x, collapse=":"))) else 
      fac <- as.factor(unlist(facs))
  } else fac <- single.factor(pfit)
  covs <- cov.extract(pfit) 
  if(is.null(covs)) X0 <- model.matrix(~fac) else {
    if(ncol(as.matrix(covs)) == 1) covs <- as.vector(covs)
    X0 <- model.matrix(~covs+fac)
  }
  X <- X0*sqrt(pfit$weights)
  if(!is.null(covs)){
    Xcov.mean <- as.matrix(X[,1:ncol(model.matrix(~covs))])
    Xcov.mean <- sapply(1:ncol(Xcov.mean), 
                        function(j) rep(mean(Xcov.mean[,j]),nrow(Xcov.mean)))
    Xnew <- cbind(Xcov.mean, X[,-(1:ncol(Xcov.mean))])
  } else Xnew <- X
  list(X0=X, X=Xnew, fac=fac)
}

quick.ls.means <- function(X0, X, Y, fac, Pcor=NULL){
  if(!is.null(Pcor)) fit <- .lm.fit(Pcor%*%X0,Y)$coefficients else fit <- .lm.fit(X0,Y)$coefficients
  lsm <- as.matrix(lm.fit(model.matrix(~fac+0),X%*%fit)$coefficients)
  rownames(lsm) <- levels(fac)
  lsm
}

apply.ls.means <- function(pfit, Yr, g=NULL, data=NULL, Pcor=NULL){
  setup<-quick.ls.means.set.up(pfit, g=g, data=data)
  X0 <- setup$X0; X <- setup$X; fac <- setup$fac
  Map(function(y) quick.ls.means(X0,X, Y=y, fac=fac, Pcor=Pcor), Yr)
}

# slopes
# estimates slopes from models with slopes and factors
# advanced.procD.lm
# THIS FUNCTION IS NO LONGER USED.  It is retained, however, for future debugging
# purposes, in the event updates have errors
slopes = function(pfit, Y=NULL, g = NULL, slope=NULL, data = NULL, Pcor = NULL){ 
  if(!is.null(Pcor) && is.null(Y)) stop("If Pcor is provided, Y cannot be null")
  if(is.null(Y)) Y <- pfit$wY
  if(is.null(data)) dat <- pfit$data else dat <- data
  if(!is.null(g)) {
    if(is.data.frame(g)){
      fac.check <- sapply(g, is.factor)
      facs <- g[,fac.check]
    } else if(is.factor(g)) facs <- g else stop("groups input neither a factor nor factors")
    if(ncol(as.matrix(facs)) > 1) fac <- factor(apply(facs, 1,function(x) paste(x, collapse=":"))) else 
      fac <- as.factor(unlist(facs))
  } else fac <- single.factor(pfit)
  if(!is.null(slope)) covs <- as.matrix(slope) else covs <- as.matrix(cov.extract(pfit))
  if(ncol(covs) > 1) stop("Only one covariate can be used for slope-comparisons")
  if(ncol(covs) == 0) stop("No covariate for which to compare slopes")
  if(!is.null(Pcor)) {
    B <- qr.coef(qr(Pcor%*%model.matrix(~fac*covs)), Y)
    } else B <- qr.coef(qr(model.matrix(~fac*covs)), Y)
  fac.p <- qr(model.matrix(~fac))$rank
  if(is.matrix(B)) {
    Bslopes <- as.matrix(B[-(1:fac.p),])
    Badjust <- rbind(0,t(sapply(2:fac.p, function(x) matrix(Bslopes[1,]))))
    Bslopes <- Bslopes + Badjust} else {
      Bslopes <- B[-(1:fac.p)]
      Badjust <- c(0,rep(Bslopes[1],(fac.p-1)))
      Bslopes <- matrix(Bslopes + Badjust)
    }
  if(ncol(Bslopes)==1) Bslopes <- cbind(1, Bslopes)
  rownames(Bslopes) <- levels(fac)
  Bslopes
}

# apply.slopes
# estimates slopes from models with slopes and factors across random outcomes in a list
# advanced.procD.lm
quick.slopes.set.up <- function(pfit, g=NULL, slope=NULL, data=NULL) {
  if(is.null(data)) dat <- pfit$data else dat <- data
  if(!is.null(g)) {
    if(is.data.frame(g)){
      fac.check <- sapply(g, is.factor)
      facs <- g[,fac.check]
    } else if(is.factor(g)) facs <- g else stop("groups input neither a factor nor factors")
    if(ncol(as.matrix(facs)) > 1) fac <- factor(apply(facs, 1,function(x) paste(x, collapse=":"))) else 
      fac <- as.factor(unlist(facs))
  } else fac <- single.factor(pfit)
  if(!is.null(slope)) covs <- as.matrix(slope) else covs <- as.matrix(cov.extract(pfit))
  if(ncol(covs) > 1) stop("Only one covariate can be used for slope-comparisons")
  if(ncol(covs) == 0) stop("No covariate for which to compare slopes")
  list(covs=covs, fac=fac)
}

quick.slopes <- function(covs, fac, Y, Pcor=NULL){
  if(!is.null(Pcor)) {
    B <- qr.coef(qr(Pcor%*%model.matrix(~fac*covs)), Y)
  } else B <- qr.coef(qr(model.matrix(~fac*covs)), Y)
  fac.p <- qr(model.matrix(~fac))$rank
  ones <- matrix(1,fac.p-1)
  if(is.matrix(B)) {
    Bslopes <- as.matrix(B[-(1:fac.p),])
    Int <- Bslopes[1,]
    Badjust <- rbind(0,ones%*%Int)
    Bslopes <- Bslopes + Badjust} else {
      Bslopes <- B[-(1:fac.p)]
      Badjust <- c(0,rep(Bslopes[1],(fac.p-1)))
      Bslopes <- matrix(Bslopes + Badjust)
    }
  if(ncol(Bslopes)==1) Bslopes <- cbind(1, Bslopes)
  rownames(Bslopes) <- levels(fac)
  Bslopes
}
  
apply.slopes <- function(pfit, Yr, g=NULL, slope=NULL, data=NULL, Pcor=NULL){
  setup <- quick.slopes.set.up(pfit, g=g, slope=slope, data=data)
  covs <- setup$covs; fac <- setup$fac
  slopes <- Map(function(y) quick.slopes(covs, fac,y, Pcor=Pcor), Yr)
  if(ncol(slopes[[1]])==1) slopes <- Map(function(s) cbind(1,s), slopes)
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

# pls
# performs PLS analysis
# Used in two.b.pls, integration.test, phylo.integration, apply.pls
pls <- function(x,y, RV=FALSE, verbose = FALSE){
  x <- as.matrix(x); y <- as.matrix(y)
  px <- dim(x)[2]; py <- dim(y)[2]; pmin <- min(px,py)
  S12 <- crossprod(center(x),center(y))/(dim(x)[1] - 1)
  if(length(S12) == 1) {
    r.pls <- cor(x,y) 
    pls <- NULL
    U <- NULL
    V <- NULL
    XScores <- x 
    YScores <- y 
    }
  else {
    pls <- La.svd(S12, pmin, pmin)
    U <- pls$u; V <- t(pls$vt)
    XScores <- x %*% U
    YScores <- y %*% V
    r.pls <- cor(XScores[,1],YScores[,1])
  }
  if(RV == TRUE){
    S11 <- var(x)
    S22 <- var(y)
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

# quick.pls
# a streamlines pls code
# used in: apply.pls
quick.pls <- function(x,y) {# no RV; no verbose output
  # assume parameters already found and assume x and y are centered
  S12 <- crossprod(x,y)/(dim(x)[1] - 1)
  if(length(S12) == 1) res <- cor(x,y) else {
    pls <- La.svd(S12, 1, 1)
    U<-pls$u; V <- as.vector(pls$vt)
    res <- cor(x%*%U,y%*%V)
  }
  res
}

# apply.pls 
# run permutations of pls analysis
# used in: two.b.pls, integration.test
apply.pls <- function(x,y, RV=FALSE, iter, seed = NULL){
  x <- as.matrix(x); y <- as.matrix(y)
  px <- ncol(x); py <- ncol(y)
  pmin <- min(px,py)
  ind <- perm.index(nrow(x), iter, seed=seed)
  RV.rand <- r.rand <- NULL
  pb <- txtProgressBar(min = 0, max = ceiling(iter/100), initial = 0, style=3) 
  jj <- iter+1
  step <- 1
  if(jj > 100) j <- 1:100 else j <- 1:jj
  while(jj > 0){
    ind.j <- ind[j]
    y.rand <-lapply(1:length(j), function(i) y[ind.j[[i]],])
    if(RV == TRUE) RV.rand <- c(RV.rand,sapply(1:length(j), function(i) pls(x,y.rand[[i]], RV=TRUE, verbose = TRUE)$RV)) else
      r.rand <- c(r.rand, sapply(1:length(j), function(i) quick.pls(x,y.rand[[i]])))
    jj <- jj-length(j)
    if(jj > 100) kk <- 1:100 else kk <- 1:jj
    j <- j[length(j)] +kk
    setTxtProgressBar(pb,step)
    step <- step+1
  }
  close(pb)
  if(RV == TRUE) RV.rand else r.rand
}

# .apply.pls 
# same as apply.pls, but without progress bar option
# used in: two.b.pls, integration.test
.apply.pls <- function(x,y, RV=FALSE, iter, seed = NULL){
  x <- as.matrix(x); y <- as.matrix(y)
  px <- ncol(x); py <- ncol(y)
  pmin <- min(px,py)
  ind <- perm.index(nrow(x), iter, seed=seed)
  RV.rand <- r.rand <- NULL
  jj <- iter+1
  if(jj > 100) j <- 1:100 else j <- 1:jj
  while(jj > 0){
    ind.j <- ind[j]
    y.rand <-lapply(1:length(j), function(i) y[ind.j[[i]],])
    if(RV == TRUE) RV.rand <- c(RV.rand,sapply(1:length(j), function(i) pls(x,y.rand[[i]], RV=TRUE, verbose = TRUE)$RV)) else
      r.rand <- c(r.rand, sapply(1:length(j), function(i) quick.pls(x,y.rand[[i]])))
    jj <- jj-length(j)
    if(jj > 100) kk <- 1:100 else kk <- 1:jj
    j <- j[length(j)] +kk
  }
  if(RV == TRUE) RV.rand else r.rand
}

# pls.multi
# obtain average of pairwise PLS analyses for 3+modules
# used in: apply.plsmulti, integration.test
plsmulti<-function(x,gps){
  g<-factor(as.numeric(gps))
  ngps<-nlevels(g)
  S<-var(x)
  gps.combo <- combn(ngps, 2)
  pls.gp <- sapply(1:ncol(gps.combo), function(j){ # no loops
    S12<-S[which(g==gps.combo[1,j]),which(g==gps.combo[2,j])]
    px <- nrow(S12); py <- ncol(S12); pmin <- min(px,py)
    pls<-La.svd(S12, pmin, pmin)
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

# quick.plsmulti
# a streamlined plsmulti
# used in apply.plsmulti
quick.plsmulti <- function(x,g,gps.combo){
  # assumed x is already centered
  pls.gp <- sapply(1:ncol(gps.combo), function(j){ # no loops
    xx <- x[,g==gps.combo[1,j]]
    yy <- x[,g==gps.combo[2,j]]
    S12<-crossprod(xx,yy)/(dim(xx)[1] - 1)
    pls<-La.svd(S12, 1, 1)
    U<-pls$u; V<-as.vector(pls$vt)
    cor(xx%*%U, yy%*%V)
  })
  mean(pls.gp) 
}

# apply.plsmulti
# permutation for multipls
# used in: integration.test
apply.plsmulti <- function(x,gps, iter, seed = NULL){
  g<-factor(as.numeric(gps))
  ngps<-nlevels(g)
  gps.combo <- combn(ngps, 2)
  ind <- perm.index(nrow(x), iter, seed=seed)
  pb <- txtProgressBar(min = 0, max = ceiling(iter/100), initial = 0, style=3) 
  jj <- iter+1
  step <- 1
  if(jj > 100) j <- 1:100 else j <- 1:jj
  r.rand <- NULL
  while(jj > 0){
    ind.j <- ind[j]
    x.r<-lapply(1:length(j), function(i) x[ind.j[[i]], g==1]) 
    r.rand<-c(r.rand, sapply(1:length(j), function(i) quick.plsmulti(cbind(x.r[[i]],
                                                  x[,g!=1]), g, gps.combo))) 
    jj <- jj-length(j)
    if(jj > 100) kk <- 1:100 else kk <- 1:jj
    j <- j[length(j)] +kk
    setTxtProgressBar(pb,step)
    step <- step+1
  }
  close(pb)
  r.rand
}

# .apply.plsmulti
# same as apply.plsmulti, but without progress bar option
# used in: integration.test
.apply.plsmulti <- function(x,gps, iter, seed = NULL){
  g<-factor(as.numeric(gps))
  ngps<-nlevels(g)
  gps.combo <- combn(ngps, 2)
  ind <- perm.index(nrow(x), iter, seed=seed)
  jj <- iter+1
  if(jj > 100) j <- 1:100 else j <- 1:jj
  r.rand <- NULL
  while(jj > 0){
    ind.j <- ind[j]
    x.r<-lapply(1:length(j), function(i) x[ind.j[[i]], g==1]) 
    r.rand<-c(r.rand, sapply(1:length(j), function(i) quick.plsmulti(cbind(x.r[[i]],
                                 x[,g!=1]), g, gps.combo)))
    jj <- jj-length(j)
    if(jj > 100) kk <- 1:100 else kk <- 1:jj
    j <- j[length(j)] +kk
  }
  r.rand
}

# CR
# Function to estimate CR coefficient
# used in: modularity.test, apply.CR
CR<-function(x,gps){
  g <- gps
  ngps <- length(unique(g))
  S<-var(x)
  diag(S)<-0
  gps.combo <- combn(ngps, 2)
  CR.gp <- sapply(1:ncol(gps.combo), function(j){ # no loops
    S11<-S[which(g==gps.combo[1,j]),which(g==gps.combo[1,j])]
    S22<-S[which(g==gps.combo[2,j]),which(g==gps.combo[2,j])]
    S12<-S[which(g==gps.combo[1,j]),which(g==gps.combo[2,j])]
    sqrt(sum(colSums(S12^2))/sqrt(sum(S11^2)*sum(S22^2)))
  })
  if(length(CR.gp) > 1) CR.mat <- dist(matrix(0, ngps,)) else 
    CR.mat = 0 # may not be necessary
  for(i in 1:length(CR.mat)) CR.mat[[i]] <- CR.gp[i]

  CR.obs <- mean(CR.gp) 
  list(CR = CR.obs, CR.mat=CR.mat)
}

# quick.CR
# streamlined CR
# used in: apply.CR, boot.CR
quick.CR <-function(x,gps){ # no CR.mat made
  g <- gps
  ngps <- length(unique(g))
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

# apply.CR
# permutation for CR
# used in: modularity.test
apply.CR <- function(x,g,k, iter, seed = NULL){# g = partition.gp
  ind <- perm.CR.index(g,k, iter, seed=seed)
  pb <- txtProgressBar(min = 0, max = ceiling(iter/100), initial = 0, style=3) 
  jj <- iter+1
  step <- 1
  if(jj > 100) j <- 1:100 else j <- 1:jj
  CR.rand <- NULL
  while(jj > 0){
    ind.j <- ind[j]
    CR.rand<-c(CR.rand, sapply(1:length(j), function(i) quick.CR(x, gps=ind.j[[i]])))
    jj <- jj-length(j)
    if(jj > 100) kk <- 1:100 else kk <- 1:jj
    j <- j[length(j)] +kk
    setTxtProgressBar(pb,step)
    step <- step+1
  }
  close(pb)
  CR.rand
}

# .apply.CR
# same as apply.CR, but without progress bar option
# used in: modularity.test
.apply.CR <- function(x,g,k, iter, seed = NULL){# g = partition.gp
  ind <- perm.CR.index(g,k, iter, seed=seed)
  jj <- iter+1
  if(jj > 100) j <- 1:100 else j <- 1:jj
  CR.rand <- NULL
  while(jj > 0){
    ind.j <- ind[j]
    CR.rand<-c(CR.rand, sapply(1:length(j), function(i) quick.CR(x, gps=ind.j[[i]])))
    jj <- jj-length(j)
    if(jj > 100) kk <- 1:100 else kk <- 1:jj
    j <- j[length(j)] +kk
  }
  CR.rand
}

# boot.CR
# Bootstrap for CR
# used in: modularity.test
boot.CR <- function(x,gps, k,iter, seed = NULL){
  x<-as.matrix(x)
  boot <- boot.index(nrow(x), iter, seed=seed)
  if(k==1){
    jj <- iter+1
    if(jj > 100) j <- 1:100 else j <- 1:jj
    CR.boot <- NULL
    while(jj > 0){
      boot.j <- boot[j]
      x.r<-lapply(1:length(j), function(i) x[boot.j[[i]],])
      CR.boot<-c(CR.boot, sapply(1:length(j), function(i) quick.CR(x.r[[i]],gps)))
      jj <- jj-length(j)
      if(jj > 100) kk <- 1:100 else kk <- 1:jj
      j <- j[length(j)] +kk
    }
  }
  
  if(k>1){  #for GM data
    angle <- seq(-44,44,2)
    if(k==2){
      rot.mat<-lapply(1:(length(angle)), function(i) matrix(c(cos(angle[i]*pi/180),
                                                              sin(angle[i]*pi/180),-sin(angle[i]*pi/180),cos(angle[i]*pi/180)),ncol=2))      
    }
    if(k==3){
      rot.mat<-lapply(1:(length(angle)), function(i) matrix(c(cos(angle[i]*pi/180),
                                                              sin(angle[i]*pi/180),0,-sin(angle[i]*pi/180),cos(angle[i]*pi/180), 0,0,0,1),ncol=3))      
    }
    jj <- iter+1
    if(jj > 100) j <- 1:100 else j <- 1:jj
    CR.boot <- NULL
    while(jj > 0){
      boot.j <- boot[j]
      x.r<-lapply(1:length(j), function(i) x[boot.j[[i]],])
      Alist.r <-lapply(1:length(x.r), function(i) { y <- x.r[[i]]; arrayspecs(y, ncol(y)/k,k)})
      CR.boot <- c(CR.boot, lapply(1:length(Alist.r), function(i){
        A <- Alist.r[[i]]
        A <- lapply(1:dim(A)[[3]], function(ii) A[,,ii])
        B <- Map(function(r) t(mapply(function(a) matrix(t(a%*%r)), A)), rot.mat)
        CRs <- Map(function(x) quick.CR(x,gps), B)
        Reduce("+", CRs)/length(CRs)
      }))
      jj <- jj-length(j)
      if(jj > 100) kk <- 1:100 else kk <- 1:jj
      j <- j[length(j)] +kk
    }
  }
  unlist(CR.boot)
}

# CR.phylo
# phylogenetic CR analysis
# used in: phylo.modularity
CR.phylo<-function(x,invC,gps){
  g <- gps
  ngps <- length(unique(g))
  gps.combo <- combn(ngps, 2)
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
  if(length(CR.gp) > 1) CR.mat <- dist(matrix(0, ngps,)) else 
    CR.mat = 0 
  for(i in 1:length(CR.mat)) CR.mat[[i]] <- CR.gp[i]
  
  CR.obs <- mean(CR.gp) 
  list(CR = CR.obs, CR.mat=CR.mat)
}

# quick.CR.phylo
# streamlined phylo CR
# used in: apply.phylo.CR
quick.CR.phylo <- function(x,invC,gps){
  x <- as.matrix(x); invC <- as.matrix(invC)
  g <- gps
  ngps <- length(unique(g))
  gps.combo <- combn(ngps, 2)
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

# apply.phylo.CR
# permutation for phylo.CR
# used in: phylo.modularity
apply.phylo.CR <- function(x,invC,gps, k, iter, seed=NULL){
  ind <- perm.CR.index(g=gps,k, iter, seed=seed)
  pb <- txtProgressBar(min = 0, max = ceiling(iter/100), initial = 0, style=3) 
  jj <- iter+1
  step <- 1
  if(jj > 100) j <- 1:100 else j <- 1:jj
  CR.rand <- NULL
  while(jj > 0){
    ind.j <- ind[j]
    CR.rand<-c(CR.rand, sapply(1:length(j), function(i) quick.CR.phylo(x,invC=invC,gps=ind.j[[i]]))) 
    jj <- jj-length(j)
    if(jj > 100) kk <- 1:100 else kk <- 1:jj
    j <- j[length(j)] +kk
    setTxtProgressBar(pb,step)
    step <- step+1
  }
  close(pb)
  CR.rand
}

# .apply.phylo.CR
# same as apply.phylo.CR, but without progress bar option
# used in: phylo.modularity
.apply.phylo.CR <- function(x,invC,gps, k, iter, seed=NULL){
  ind <- perm.CR.index(g=gps,k, iter, seed=seed)
  jj <- iter+1
  if(jj > 100) j <- 1:100 else j <- 1:jj
  CR.rand <- NULL
  while(jj > 0){
    ind.j <- ind[j]
    CR.rand<-c(CR.rand, sapply(1:length(j), function(i) quick.CR.phylo(x,invC=invC,gps=ind.j[[i]]))) 
    jj <- jj-length(j)
    if(jj > 100) kk <- 1:100 else kk <- 1:jj
    j <- j[length(j)] +kk
  }
  CR.rand
}

# boot.phylo.CR
# bootstrap for phylo.CR
# phylo.modularity 
boot.phylo.CR <- function(x, invC, gps, k,iter, seed = NULL){
  x<-as.matrix(x)
  boot <- boot.index(nrow(x), iter, seed=seed)
  if(k==1){
    jj <- iter+1
    if(jj > 100) j <- 1:100 else j <- 1:jj
    CR.boot <- NULL
    while(jj > 0){
      boot.j <- boot[j]
      x.r<-lapply(1:length(j), function(i) x[boot.j[[i]],])
      invC.r <- lapply(1:length(j), function(i) invC[boot.j[[i]],boot.j[[i]]])
      CR.boot<-c(CR.boot, sapply(1:length(j), function(i) quick.CR.phylo(x=x.r[[i]],invC=invC.r[[i]],gps=gps)))
      jj <- jj-length(j)
      if(jj > 100) kk <- 1:100 else kk <- 1:jj
      j <- j[length(j)] +kk
    }
  }
  
  if(k>1){  #for GM data
    angle <- seq(-44,44,2)
    if(k==2){
      rot.mat<-lapply(1:(length(angle)), function(i) matrix(c(cos(angle[i]*pi/180),
                                                              sin(angle[i]*pi/180),-sin(angle[i]*pi/180),cos(angle[i]*pi/180)),ncol=2))      
    }
    if(k==3){
      rot.mat<-lapply(1:(length(angle)), function(i) matrix(c(cos(angle[i]*pi/180),
                                                              sin(angle[i]*pi/180),0,-sin(angle[i]*pi/180),cos(angle[i]*pi/180), 0,0,0,1),ncol=3))      
    }
    jj <- iter+1
    if(jj > 100) j <- 1:100 else j <- 1:jj
    CR.boot <- NULL
    while(jj > 0){
      boot.j <- boot[j]
      x.r<-lapply(1:length(j), function(i) x[boot.j[[i]],])
      Alist.r <-lapply(1:length(x.r), function(i) { y <- x.r[[i]]; arrayspecs(y, ncol(y)/k,k)})
      invC.r <- lapply(1:length(j), function(i) invC[boot.j[[i]],boot.j[[i]]])
      CR.boot <- c(CR.boot, lapply(1:length(Alist.r), function(i){
        A <- Alist.r[[i]]
        A <- lapply(1:dim(A)[[3]], function(ii) A[,,ii])
        B <- Map(function(r) t(mapply(function(a) matrix(t(a%*%r)), A)), rot.mat)
        CRs <- Map(function(x) quick.CR.phylo(x=x.r[[i]],invC=invC.r[[i]],gps=gps), B)
        Reduce("+", CRs)/length(CRs)
      }))
      jj <- jj-length(j)
      if(jj > 100) kk <- 1:100 else kk <- 1:jj
      j <- j[length(j)] +kk
    }
  }
  unlist(CR.boot)
}

# phylo.mat 
# estimate BM phylo.cov.matrix and transform from phylogeny
# used in: compare.evol.rates, compare.multi.evol.rates, phylo.integration, phylo.modularity, physignal
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
  D.mat <- fast.solve(eigC.vect%*% diag(sqrt(lambda)) %*% t(eigC.vect)) 
  rownames(D.mat) <- colnames(D.mat) <- colnames(C)
  list(invC = invC, D.mat = D.mat,C = C)
}

# pls.phylo
# phylogenetic pls
# used in: phylo.integration, apply.pls.phylo
pls.phylo <- function(x,y, Ptrans, verbose = FALSE){
  x <- as.matrix(x); y <- as.matrix(y)
  px <- ncol(x); py <- ncol(y); pmin <- min(px,py)
  x <- Ptrans%*%x
  y <- Ptrans%*%y
  R12<-  crossprod(x,y)/(nrow(x)-1)
  pls <- La.svd(R12, pmin, pmin)
  U <- pls$u; V <- t(pls$vt)
  XScores <- x %*% U 
  YScores <- y %*% V
  r.pls <- cor(XScores[, 1], YScores[, 1]) 
  if(verbose==TRUE){
    XScores <- as.matrix(XScores); Y <- as.matrix(YScores)
    rownames(U)  = colnames(x); rownames(V) = colnames(y)
    out <- list(pls.svd = pls, r.pls = r.pls, left.vectors=U, 
                right.vectors=V, XScores=XScores,YScores=YScores)
  } else out <- r.pls
  out
}

# apply.pls.phylo
# permutation for phylo.pls
# used in: phylo.integration
# CURRENTLY NOT IN USE - USING apply.pls BECAUSE OF TYPE I ERROR ISSUES
apply.pls.phylo <- function(x,y,Ptrans, iter, seed = NULL){
  n.x <- dim(x)[2]
  ind <- perm.index(nrow(x), iter, seed=seed)
  x <- Ptrans%*%x
  pb <- txtProgressBar(min = 0, max = ceiling(iter/100), initial = 0, style=3) 
  jj <- iter+1
  step <- 1
  if(jj > 100) j <- 1:100 else j <- 1:jj
  r.rand <- NULL
  while(jj > 0){
    ind.j <- ind[j]
    y.rand <-lapply(1:length(j), function(i) y[ind.j[[i]],])
    r.rand <- c(r.rand, sapply(1:length(j), function(i) {
      yy <- Ptrans%*%y.rand[[i]] 
      quick.pls(x,yy)
    }))
    jj <- jj-length(j)
    if(jj > 100) kk <- 1:100 else kk <- 1:jj
    j <- j[length(j)] +kk
    setTxtProgressBar(pb,step)
    step <- step+1
  }
  close(pb)
  r.rand
}

# .apply.pls.phylo
# same as apply.phylo.pls, but without progress bar option
# used in: phylo.integration
# CURRENTLY NOT IN USE - USING .apply.pls BECAUSE OF TYPE I ERROR ISSUES
.apply.pls.phylo <- function(x,y,Ptrans,iter, seed = NULL){
  n.x <- dim(x)[2]
  ind <- perm.index(nrow(x), iter, seed=seed)
  x <- Ptrans%*%x
  jj <- iter+1
  if(jj > 100) j <- 1:100 else j <- 1:jj
  r.rand <- NULL
  while(jj > 0){
    ind.j <- ind[j]
    y.rand <-lapply(1:length(j), function(i) y[ind.j[[i]],])
    r.rand <- c(r.rand, sapply(1:length(j), function(i) {
      yy <- Ptrans%*%y.rand[[i]] 
      quick.pls(x,yy)
    }))
    jj <- jj-length(j)
    if(jj > 100) kk <- 1:100 else kk <- 1:jj
    j <- j[length(j)] +kk
  }
  r.rand
}

# plsmulti.phylo
# average pairwise phylo.pls
# used in: phylo.integration, apply.plsmulti.phylo
plsmulti.phylo<-function(x,gps, Ptrans){
  g<-factor(as.numeric(gps))
  ngps<-nlevels(g)
  x <- Ptrans%*%x
  gps.combo <- combn(ngps, 2)
  R<-  crossprod(x) * (nrow(x)-1)^-1 
  pls.gp <- sapply(1:ncol(gps.combo), function(j){ 
    R12<-R[which(g==gps.combo[1,j]),which(g==gps.combo[2,j])]
    px <- nrow(R12); py <- ncol(R12); pmin <- min(px,py)
    pls<-La.svd(R12, pmin, pmin)
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

# apply.plsmulti.phylo
# permutations for plsmulti.phylo
# used in: phylo.integration
# CURRENTLY NOT IN USE - USING apply.plsmulti BECAUSE OF TYPE I ERROR ISSUES
apply.plsmulti.phylo <- function(x, gps,Ptrans, iter=iter, seed=seed){
  g<-factor(as.numeric(gps))
  ngps<-nlevels(g)
  gps.combo <- combn(ngps, 2)
  xx <- Ptrans%*%x[,g==1]; yy <- Ptrans%*%x[,g!=1]
  ind <- perm.index(nrow(x), iter, seed=seed)
  pb <- txtProgressBar(min = 0, max = ceiling(iter/100), initial = 0, style=3) 
  jj <- iter+1
  step <- 1
  if(jj > 100) j <- 1:100 else j <- 1:jj
  r.rand <- NULL
  while(jj > 0){
    ind.j <- ind[j]
    x.r <-lapply(1:length(j), function(i) xx[ind.j[[i]],])
    r.rand <- c(r.rand, sapply(1:length(j), function(i) {
      quick.pls(x.r[[i]],yy)
      }))
    jj <- jj-length(j)
    if(jj > 100) kk <- 1:100 else kk <- 1:jj
    j <- j[length(j)] +kk
    setTxtProgressBar(pb,step)
    step <- step+1
  }
  close(pb)
  r.rand
}

# .apply.plsmulti.phylo
# same as apply.plsmulti.phylo, but without progress bar option
# used in: phylo.integration
# CURRENTLY NOT IN USE - USING .apply.plsmulti BECAUSE OF TYPE I ERROR ISSUES
.apply.plsmulti.phylo <- function(x, gps,Ptrans, iter=iter, seed=seed){
  g<-factor(as.numeric(gps))
  ngps<-nlevels(g)
  gps.combo <- combn(ngps, 2)
  xx <- Ptrans%*%x[,g==1]; yy <- Ptrans%*%x[,g!=1]
  ind <- perm.index(nrow(x), iter, seed=seed)
  jj <- iter+1
  if(jj > 100) j <- 1:100 else j <- 1:jj
  r.rand <- NULL
  while(jj > 0){
    ind.j <- ind[j]
    x.r <-lapply(1:length(j), function(i) xx[ind.j[[i]],])
    r.rand <- c(r.rand, sapply(1:length(j), function(i) {
      quick.pls(x.r[[i]],yy)
    }))
    jj <- jj-length(j)
    if(jj > 100) kk <- 1:100 else kk <- 1:jj
    j <- j[length(j)] +kk
  }
  r.rand
}

# sigma.d
# multivariate evolutionary rate
# used in: compare.evol.rates
sigma.d<-function(x,invC,D.mat,gp){
  N<-dim(x)[1];p<-dim(x)[2]
  g<-factor(as.numeric(gp))
  ngps<-nlevels(g)
  ones<-matrix(1,N,N) 
  x.c<-x-crossprod(ones,invC)%*%x/sum(invC) 
  R<-crossprod(x.c, crossprod(invC,x.c))/N
  vec.d2<-diag(tcrossprod(D.mat%*%(x.c)))
  sigma.d.all<-sum(vec.d2)/N/p
  sigma.d.gp<-sapply(split(vec.d2, gp), mean)/p  
  sigma.d.ratio<-sigma.d.rat<-sigma.d.rat.mat<-rate.mat<-NULL
  if(ngps>1){
    gps.combo <- combn(ngps, 2)
    sigma.d.rat <- sapply(1:ncol(gps.combo), function(j){ 
      rates<-c(sigma.d.gp[levels(g)==gps.combo[1,j]],sigma.d.gp[levels(g)==gps.combo[2,j]])
      rate.rat<-max(rates)/min(rates)
    })
    if(length(sigma.d.rat) > 1) rate.mat <- dist(matrix(0, length(sigma.d.gp),)) else 
      rate.mat = 0 
    for(i in 1:length(rate.mat)) rate.mat[[i]] <- sigma.d.rat[i]
    sigma.d.ratio<-max(rate.mat)
  }
  list(sigma.d.ratio = sigma.d.ratio, sigma.d.all = sigma.d.all, 
       sigma.d.gp = sigma.d.gp, sigma.d.gp.ratio = rate.mat,R = R)  
}

# fast.sigma.d
# same as sigma.d but only calculates sigma.d.ratio - fast in loops
# used in: compare.evol.rates
fast.sigma.d<-function(x,Ptrans,g, ngps, gps.combo, N,p){
  xp <- Ptrans%*%x
  vec.d2<-diag(tcrossprod(xp))
  sigma.d.all<-sum(vec.d2)/N/p
  sigma.d.gp<-sapply(split(vec.d2, g), mean)/p  
  sigma.d.ratio<-sigma.d.rat<-sigma.d.rat.mat<-rate.mat<-NULL
  sigma.d.rat <- sapply(1:ncol(gps.combo), function(j){ 
    rates<-c(sigma.d.gp[levels(g)==gps.combo[1,j]],sigma.d.gp[levels(g)==gps.combo[2,j]])
    max(rates)/min(rates)
  })
  sigma.d.rat
}

# sigma.d.multi
# multiple trait multivariate evolutionary rates
# used in: compare.multi.evol.rates
sigma.d.multi<-function(x,invC,D.mat,gps,Subset){
  sig.calc<-function(x.i,invC.i,D.mat.i,Subset){
    x.i<-as.matrix(x.i)
    N<-dim(x.i)[1];p<-dim(x.i)[2]
    ones<-matrix(1,N,N) 
    x.c<- x.i - crossprod(ones,invC.i)%*%x.i/sum(invC.i) 
    R<-crossprod(x.c, crossprod(invC.i,x.c))/N
    if(Subset==FALSE) sigma<-sigma<-sum((D.mat.i%*%x.c)^2)/N  else 
      sigma<-sum((D.mat.i%*%x.c)^2)/N/p
    return(list(sigma=sigma,R=R))
  }
  g<-factor(as.numeric(gps))
  ngps<-nlevels(g)  
  gps.combo <- combn(ngps, 2)
  global<-sig.calc(x,invC,D.mat,Subset)
  rate.global<-global$sigma; R<-global$R
  ngps<-nlevels(gps)
  rate.gps<-sapply(1:ngps, function(j){ sig.calc(x[,g==j],
                                                 invC,D.mat,Subset)$sigma  })
  sigma.d.ratio<-max(rate.gps)/min(rate.gps)
  sigma.d.rat <- sapply(1:ncol(gps.combo), function(j){ 
    rates<-c(rate.gps[levels(g)==gps.combo[1,j]],rate.gps[levels(g)==gps.combo[2,j]])
    max(rates)/min(rates)
  })
  if(length(sigma.d.rat) > 1) rate.mat <- dist(matrix(0, length(rate.gps),)) else 
    rate.mat = 0 
  for(i in 1:length(rate.mat)) rate.mat[[i]] <- sigma.d.rat[i]
  list(sigma.d.ratio = sigma.d.ratio, rate.global = rate.global, 
       rate.gps = rate.gps, sigma.d.gp.ratio = rate.mat,R = R)  
}


##Fast version of compare.multi.rates for permutations

sig.calc<-function(x.i,Ptrans.i,Subset, N, p){
  xp.i <- Ptrans.i%*%x.i
  if(Subset==FALSE) sigma<-sum((xp.i)^2)/N else 
    sigma<-sum((xp.i)^2)/N/p
  return(sigma)
}

fast.sigma.d.multi<-function(x,Ptrans,Subset, gindx, ngps, gps.combo, N, p){ 
  rate.gps<-lapply(1:ngps, function(j){ sig.calc(x[,gindx[[j]]], Ptrans,Subset, N, p)})
  sapply(1:ncol(gps.combo), function(j){ 
    a <- gps.combo[1,j]
    b <- gps.combo[2,j]
    rates<-c(rate.gps[[a]],rate.gps[[b]])
    max(rates)/min(rates)
  })
}

# trajset.int
# set-up trajectories from a model with an interaction
# used in: trajectory.analysis
trajset.int <- function(y, tp,tn) { # assume data in order of variables within points
  pt.list <- sort(rep(1:tn,tp))
  lapply(1:tn, function(j){
    y[which(pt.list == j),]
  })
}

# trajset.gps
# set-up trajectories from a model without an interaction: assumes data are trajectories
# used in: trajectory.analysis
trajset.gps <- function(y, traj.pts) { # assume data in order of variables within points
  lapply(1:nrow(y), function(j) {
    matrix(y[j,], traj.pts, , byrow=T)
  })
}

# trajsize
# find path distance of trajectory
# used in: trajectory.analysis
trajsize <- function(y) {
  k <- nrow(y[[1]])
  tpairs <- cbind(1:(k-1),2:k)
  sapply(1:length(y), function(j) {
    d <- as.matrix(dist(y[[j]]))
    sum(d[tpairs])
  })
}

# trajorient
# sfind trajectory correlations from first PCs
# used in: trajectory.analysis
trajorient <- function(y, tn, p) {
  m <- t(sapply(1:tn, function(j){
    x <- y[[j]]
    La.svd(center.scale(x)$coords, 0, 1)$vt
  }))
  vec.cor.matrix(m)
}

# trajshape
# find shape differences among trajectories
# used in: trajectory.analysis
trajshape <- function(y){
  y <- Map(function(x) center.scale(x)$coords, y)
  M <- Reduce("+",y)/length(y)
  z <- apply.pPsup(M, y)
  z <- t(sapply(z, function(x) as.vector(t(x))))
  as.matrix(dist(z))
}

# traj.w.int
# full PTA stats for trajectories from a model with an interaction
# used in: trajectory.analysis
traj.w.int <- function(ff, fr, data=NULL, iter, seed= NULL){
  pfitf <- procD.fit(ff, data = data, pca=FALSE)
  pfitr <- procD.fit(fr, data = data, pca=FALSE)
  ex.terms <- length(pfitf$term.labels) - 3
  Y <- pfitf$wY
  Yh <- pfitr$fitted[[length(pfitr$fitted)]]
  E <-  pfitr$residuals[[length(pfitr$residuals)]]
  n <- nrow(Y)
  tp <- nlevels(data[[match(attr(terms(ff),"term.labels")[[ex.terms+2]], names(data))]])
  group.levels <- levels(data[[match(attr(terms(ff),"term.labels")[[ex.terms+1]], names(data))]])
  tn <- length(group.levels)
  p <- ncol(Y)
  ind <- perm.index(n, iter, seed=seed)
  Yr <- Map(function(x) Yh + E[x,], ind)
  if(sum(pfitf$weights) != n) Yr <- Map(function(y) y*sqrt(pfitf$weights), Yr)
  means <- apply.ls.means(pfitf, Yr)
  trajs <- lapply(means, function(x) trajset.int(x, tp, tn))
  PD <- lapply(trajs, trajsize) # path distances
  MD <- lapply(PD, function(x) as.matrix(dist(x)))
  Tcor <- lapply(trajs, function(x) trajorient(x,tn)) # trajectory correlations
  Tang <- lapply(Tcor, acos) # trajectory angles
  SD <- lapply(trajs, trajshape) # trajectory shape differences
  
  list(means = means, trajectories = trajs, 
       PD = PD,
       MD=MD, Tcor=Tcor, Tang=Tang, SD=SD,
       P.MD = Pval.matrix(simplify2array(MD)),
       P.angle = Pval.matrix(simplify2array(Tang)),
       P.SD = Pval.matrix(simplify2array(SD)),
       Z.MD = Effect.size.matrix(simplify2array(MD)),
       Z.angle = Effect.size.matrix(simplify2array(Tang)),
       Z.SD = Effect.size.matrix(simplify2array(SD)),
       ngroups = tn, npoints = tp)
}

# traj.by.groups
# full PTA stats for trajectories from a model without an interaction; assume data are trajectories
# used in: trajectory.analysis
traj.by.groups <- function(ff, fr, traj.pts, data=NULL, iter, seed= NULL){
  pfitf <- procD.fit(ff, data = data, pca=FALSE)
  pfitr <- procD.fit(fr, data = data, pca=FALSE)
  ex.terms <- length(pfitf$term.labels) - 1
  Y <- pfitf$wY
  Yh <- pfitr$fitted[[length(pfitr$fitted)]]
  E <-  pfitr$residuals[[length(pfitr$residuals)]]
  n <- nrow(Y)
  tp <- traj.pts
  p <- ncol(Y)/traj.pts
  if(p != floor(p)) stop("The number of variables divided by the number of trajectory points is not an integer")
  gps <- data[[match((pfitf$term.labels)[ex.terms+1], names(data))]]
  group.levels <- levels(gps)
  tn <- length(group.levels)
  ind <- perm.index(n, iter, seed=seed)
  Yr <- Map(function(x) Yh + E[x,], ind)
  if(sum(pfitf$weights) != n) Yr <- Map(function(y) y*sqrt(pfitf$weights), Yr)
  means <- apply.ls.means(pfitf, Yr)
  trajs <- lapply(means, function(x) trajset.gps(x,traj.pts))
  PD <- lapply(trajs, trajsize) # path distances
  MD <- lapply(PD, function(x) as.matrix(dist(x)))
  Tcor <- lapply(trajs, function(x) trajorient(x,tn)) # trajectory correlations
  Tang <- lapply(Tcor, acos) # trajectory angles
  SD <- lapply(trajs, trajshape) # trajectory shape differences
  
  list(means = means, trajectories = trajs, 
       PD = PD,
       MD=MD, Tcor=Tcor, Tang=Tang, SD=SD,
       P.MD = Pval.matrix(simplify2array(MD)),
       P.angle = Pval.matrix(simplify2array(Tang)),
       P.SD = Pval.matrix(simplify2array(SD)),
       Z.MD = Effect.size.matrix(simplify2array(MD)),
       Z.angle = Effect.size.matrix(simplify2array(Tang)),
       Z.SD = Effect.size.matrix(simplify2array(SD)),
       ngroups = tn, npoints = tp)
}

#####-----------------------------------------------------------------------------------

### geomorph-specific logicals

is.gpagen <- function(x) class(x) == "gpagen"
is.phylo <- function(x) class(x) == "phylo"
is.geomorph.data.frame <- function(x) class(x) == "geomorph.data.frame"

#####-----------------------------------------------------------------------------------

### geomorph-specific S3 for internal use (copies of base functions)

droplevels.geomorph.data.frame <- function (x, except = NULL, ...) {
  ix <- vapply(x, is.factor, NA)
  if (!is.null(except)) 
    ix[except] <- FALSE
  x[ix] <- lapply(x[ix], factor)
  x
}

#####-----------------------------------------------------------------------------------

### retained from old geomorph support code
### need to update and merge, or replace with new functions


scan.to.ref<-function(scandata,specland,refland){  	#DCA
  ref.scan<-tps2d3d(scandata,specland,refland)
  ref.scan}

refscan.to.spec<-function(refscan,refland,specland){ 	#DCA
  unwarp.scan<-tps2d3d(refscan,refland,specland)
  unwarp.scan}

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
##Function to read tps file for digitize2d (streamlined for specific use)
readland.tps2 <- function (file, specID = c("None", "ID", "imageID")) 
{
  ignore.case = TRUE
  specID <- match.arg(specID)
  tpsfile <- scan(file = file, what = "char", sep = "\n", quiet = TRUE)
  lmdata <- grep("LM=", tpsfile, ignore.case)
  if (length(lmdata !=0)) {
    nland <- as.numeric(sub("LM=", "", tpsfile[lmdata], ignore.case))
    k <- 2
  }
  if (length(lmdata) == 0) {
    lmdata <- grep("LM3=", tpsfile, ignore.case)
    nland <- as.numeric(sub("LM3=", "", tpsfile[lmdata], ignore.case))
    k <- 3
  }
  n <- nspecs <- length(lmdata)
  if (max(nland) - min(nland) != 0) {
    stop("Number of landmarks not the same for all specimens.")
  }
  p <- nland[1]
  imscale <- as.numeric(sub("SCALE=", "", tpsfile[grep("SCALE", 
                                                       tpsfile, ignore.case)], ignore.case))
  if (is.null(imscale)) {
    imscale = array(0, nspecs)
  }
  if (length(imscale)==0) {
    imscale = array(0, nspecs)
  }
  if (length(imscale) != nspecs) {
    imscale = array(1, nspecs)
  }
  tmp <- tpsfile[-(grep("=", tpsfile))]
  options(warn = -1)
  tmp <- matrix(as.numeric(unlist(strsplit(tmp,"\\s+"))),ncol = k, byrow = T)

  coords <- aperm(array(t(tmp), c(k, p, n)), c(2, 1, 3))
  #  imscale <- aperm(array(rep(imscale, p * k), c(n, k, p)), c(3, 2, 1))
  #  coords <- coords * imscale
  coords<-coords[1:nland,,] 
  if(n==1) coords <- array(coords, c(nland,k,n))
  if (specID == "imageID") {
    imageID <- (sub("IMAGE=", "", tpsfile[grep("IMAGE", tpsfile, ignore.case)], 
                    ignore.case))
    if (length(imageID) != 0) {
      imageID <- sub(".jpg", "", imageID, ignore.case)
      imageID <- sub(".tif", "", imageID, ignore.case)
      imageID <- sub(".bmp", "", imageID, ignore.case)
      imageID <- sub(".tiff", "", imageID, ignore.case)
      imageID <- sub(".jpeg", "", imageID, ignore.case)
      imageID <- sub(".jpe", "", imageID, ignore.case)
      dimnames(coords)[[3]] <- as.list(imageID)
    }
  }
  if (specID == "ID") {
    ID <- sub("ID=", "", tpsfile[grep("ID", tpsfile, ignore.case)], ignore.case)
    if (length(ID) != 0) {
      dimnames(coords)[[3]] <- as.list(ID)
    }
  }
  return(list(coords = coords,scale=imscale)  )                  
}
