#' @name geomorph-package
#' @docType package
#' @aliases geomorph
#' @title Geometric morphometric analyses for 2D/3D data
#' @author Dean C. Adams, Michael Collyer & Antigoni Kaliontzopoulou
#' @description Functions in this package allow one to read, manipulate, and digitize landmark data; generate shape
#'  variables via Procrustes analysis for points, curves and surface data, perform statistical analyses
#'  of shape variation and covariation, and provide graphical depictions of shapes and patterns of
#'  shape variation.
#'
#' @import RRPP
#' @import rgl
#' @import stats
#' @import utils
#' @import graphics
#' @import grDevices
#' @importFrom jpeg readJPEG
#' @importFrom ape multi2di.phylo
#' @importFrom ape root.phylo

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
#' change for phenotypes described by high-dimensional data. Heredity. 115: 357-365.
NULL


#' Head and tail shapes of larval salamanders
#'
#' @name larvalMorph
#' @docType data
#' @author Michael Collyer
#' @keywords datasets
#' @description Landmark data from heads and tails of larval Salamanders exposed to different treatments of herbicides.
#' @details
#' Data set includes tail landmarks (coords), index for identifying semilandmarks (sliders), and vectors
#' for variables including herbicide treatment (Treatment) and Family (Family).  The latter variable indicates
#' the egg masses (clutches) sampled in the wild from which individual eggs were randomly assigned to treatment.
#' See Levis et al. (2016) for more experimental details.
#' @references Levis, N.A, M.L. Schooler, J.R. Johnson, and M.L. Collyer. 2016. The effects of terrestrial and aquatic herbicides on
#' larval salamander morphology and swim speed. Biological Journal of the Linnean Society. 118:569-581.
NULL

#' Estimate mean shape for a set of aligned specimens
#'
#' Estimate the mean shape for a set of aligned specimens
#'
#' The function estimates the average landmark coordinates for a set of aligned specimens. It is assumed
#' that the landmarks have previously been aligned using Generalized Procrustes Analysis (GPA)
#'  [e.g., with \code{\link{gpagen}}]. This function is described in Claude (2008).
#'  
#'  One can then use the generic function \code{\link{plot}} to produce a numbered plot of landmark 
#'  positions and potentially add links, in order to review landmark positions
#'
#' @param A Either a list (length n, p x k), A 3D array (p x k x n), or a matrix (n x pk) containing GPA-aligned coordinates for a set of specimens
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
  if(is.array(A)) {
    dims <- dim(A)
    if(length(dims) == 3) res <- apply(A,c(1,2),mean) else
      if(length(dims) == 2){
        if(dims[[2]] == 2 || dims[[2]] == 3) res <- A else
        {
          cat("\nWarning: It appears that data are in a matrix with specimens as rows.")
          cat("\nMeans are found for each column of the matrix.\n\n")
          res <- colMeans(A)
        }
      }
  }
  if(is.list(A)) res <- Reduce("+", A)/length(A)
  if(!is.array(A) && !is.list(A) && !is.matrix(A)) stop("There are not multiple configurations from which to obtain a mean.")
  class(res) <- c("mshape", "matrix")
  return(res)
}

#####----------------------------------------------------------------------------------------------------

# SUPPORT FUNCTIONS

# scanTPS
# Scans data and other info from TPS files
# used by readland.tps
scanTPS <- function(file) {
  ignore.case = TRUE
  tpsfile <- scan(file = file, what = "char", sep = "\n", quiet = TRUE)
  commline <- grep("COMMENT=", tpsfile, ignore.case)
  if(length(commline) != 0) tpsfile <- tpsfile[-commline] # removes COMMENT= lines
  lmline <- grep("LM=", tpsfile, ignore.case)
  if(length(lmline) == 0) lmline <- grep("LM3=", tpsfile, ignore.case)
  if(length(lmline) == 0) stop("No landmark command provided; e.g., 'LM=' or 'LM3=")
  endline <- lmline[-1] - 1
  endline <- c(endline, length(tpsfile))
  n <- length(lmline)
  spec.list <- lapply(1:n, function(j){
    start <- lmline[j]
    end <- endline[j]
    temp <- tpsfile[start:end]
    lml <- grep("LM", temp)
    crvl <- grep("CURVES", temp)
    cptl <- grep("POINTS", temp)
    scl <- grep("SCALE", temp)
    iml <- grep("IMAGE", temp)
    idl <- grep("ID", temp)
    notlm <- grep("=", temp)
    templm <- strsplit(temp[-notlm], "\\s+")
    lm <- lapply(templm, as.numeric)
    p <- length(lm)
    k <- sapply(lm, length)
    if(length(unique(k)) == 1) k <- unique(k)
    scale <- as.numeric(unlist(strsplit(temp[scl], "SCALE="))[2])
    id <- unlist(strsplit(temp[idl], "ID="))[2]
    image <- unlist(strsplit(temp[iml], "IMAGE="))[2]
    image <- sub(".jpg", "", image, ignore.case)
    image <- sub(".tif", "", image, ignore.case)
    image <- sub(".bmp", "", image, ignore.case)
    image <- sub(".tiff", "", image, ignore.case)
    image <- sub(".jpeg", "", image, ignore.case)
    image <- sub(".jpe", "", image, ignore.case)
    plm <- as.numeric(unlist(strsplit(temp[lml], "="))[2])
    pcv <- p - plm
    if(p > plm) {
      curve.lm <- lm[-(1:plm)] 
      lm <- lm[1:plm]
      curve.pts <- as.vector(na.omit(as.numeric(unlist(strsplit(temp[cptl], "POINTS=")))))
    } else curve.lm <- curve.pts <- NULL
    
    out <- list(lm = lm, curve.lm = curve.lm, p = p, plm = plm, pcv = pcv,
                k = k, curve.pts = curve.pts, scale = scale, id = id, image = image)
    out  
  })
}

# orp
# projection in GPA
# used in gpagen functions
orp<-function(A){
  if(is.array(A)) {
    dims <- dim(A)
    n <- dims[3]; k <- dims[2]; p <- dims[1]
    Y <- lapply(1:n, function(j) A[,,j])
  } else
    if(is.list(A)){
      Y <- A
      n <- length(A); dims <- dim(A[[1]]); k <- dims[2]; p <- dims[1]
    } else stop("Input must be either a list or array")

  Y1 <- as.vector(center.scale((Reduce("+", Y)/n))$coords)
  oo <- matrix(1,n)%*%Y1
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

# tangents
# finds tangents in a matrix based on sliders
# used in all functions associated with pPga.wCurves
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
nearest <- function(X, m, k = 4) {
  a <- X[m,]
  b <- sapply(1:nrow(X), function (j) sum((a-X[j,])^2))
  match(sort(b)[2:(k + 1)], b)
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
# used in BE.slide
Ltemplate <-function(Mr, Mt=NULL){
  p <-nrow(Mr); k <- ncol(Mr)
  if(!is.null(Mt)) P <- as.matrix(dist(Mr-Mt)) else P <- as.matrix(dist(Mr))
  if(k==2) {P <-P^2*log(P); P[is.na(P)] <- 0}
  Q <- cbind(1,Mr)
  L<-rbind(cbind(P,Q), cbind(t(Q),matrix(0,k+1,k+1)))
  Linv <- -fast.solve(L)[1:p,1:p]
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
  Linv <- -fast.solve(L)[1:p,1:p]
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
  n <- length(Y); dims <- dim(Y[[1]]); p <- dims[1]; k <- dims[2]
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
  k <- ncol(y)
  kk <- round(0.05 * length(surf))
  kk <- max(c(k, kk))
  p <- nrow(y)
  pc.match <- 1:p; pc.match[-surf] = NA
  nearpts <- lapply(1:p, function(j) {
    nn <- pc.match[j]
    if(is.na(nn)) 0 else
      c(nearest(y, nn, k = kk+1), nn)})
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
    int.part <- fast.solve(t(t(tx*L)*tx)+t(t(ty*L)*ty)+
                            t(t(tz*L)*tz))%*%cbind(tx*L,ty*L,tz*L)
    Ht <- rbind(tx*int.part, ty*int.part, tz*int.part)
  } else {
    int.part <- fast.solve(t(t(tx*L)*tx)+t(t(ty*L)*ty))%*%cbind(tx*L,ty*L)
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
    int.part <- fast.solve(t(t(p1x*L)*p1x)+t(t(p1y*L)*p1y)+
                            t(t(p1z*L)*p1z))%*%cbind(p1x*L,p1y*L,p1z*L)
    Hp1 <- rbind(p1x*int.part, p1y*int.part, p1z*int.part)
  } else {
    int.part <- fast.solve(t(t(p1x*L)*p1x)+t(t(p1y*L)*p1y))%*%cbind(p1x*L,p1y*L)
    Hp1 <- rbind(p1x*int.part, p1y*int.part)
  }
  if(k==3) {
    int.part <- fast.solve(t(t(p2x*L)*p2x)+t(t(p2y*L)*p2y)+
                            t(t(p2z*L)*p2z))%*%cbind(p2x*L,p2y*L,p2z*L)
    Hp2 <- rbind(p2x*int.part, p2y*int.part, p2z*int.part)
  } else {
    int.part <- fast.solve(t(t(p2x*L)*p2x)+t(t(p2y*L)*p2y))%*%cbind(p2x*L,p2y*L)
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
    int.part <- fast.solve(t(t(tx*L)*tx)+t(t(ty*L)*ty)+
                            t(t(tz*L)*tz))%*%cbind(tx*L,ty*L,tz*L)
    Ht <- rbind(tx*int.part, ty*int.part, tz*int.part)
  } else {
    int.part <- fast.solve(t(t(tx*L)*tx)+t(t(ty*L)*ty))%*%cbind(tx*L,ty*L)
    Ht <- rbind(tx*int.part, ty*int.part)
  }
  PC <- getSurfPCs(y, surf)
  p1x <- PC$p1x; p1y <- PC$p1y; p1z <- PC$p1z; p2x <- PC$p2x; p2y <- PC$p2y; p2z <- PC$p2z
  if(k==3) {
    int.part <- fast.solve(t(t(p1x*L)*p1x)+t(t(p1y*L)*p1y)+
                            t(t(p1z*L)*p1z))%*%cbind(p1x*L,p1y*L,p1z*L)
    Hp1 <- rbind(p1x*int.part, p1y*int.part, p1z*int.part)
  } else {
    int.part <- fast.solve(t(t(p1x*L)*p1x)+t(t(p1y*L)*p1y))%*%cbind(p1x*L,p1y*L)
    Hp1 <- rbind(p1x*int.part, p1y*int.part)
  }
  if(k==3) {
    int.part <- fast.solve(t(t(p2x*L)*p2x)+t(t(p2y*L)*p2y)+
                            t(t(p2z*L)*p2z))%*%cbind(p2x*L,p2y*L,p2z*L)
    Hp2 <- rbind(p2x*int.part, p2y*int.part, p2z*int.part)
  } else {
    int.part <- fast.solve(t(t(p2x*L)*p2x)+t(t(p2y*L)*p2y))%*%cbind(p2x*L,p2y*L)
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
tps <- function(matr, matt, n, sz=1.5, pt.bg="black",
              grid.col="black", grid.lwd=1, grid.lty=1, refpts=FALSE, k3 = FALSE){		#DCA: altered from J. Claude: 2D only
  xm <- min(matr[,1])
  ym <- min(matr[,2])
  xM <- max(matr[,1])
  yM <- max(matr[,2])
  rX <- xM-xm; rY <- yM-ym
  a <- seq(xm - 1/5*rX, xM + 1/5*rX, length=n)
  b <- seq(ym - 1/5*rX, yM + 1/5*rX, by=(xM-xm)*7/(5*(n-1)))
  m <- round(0.5+(n-1)*(2/5*rX + yM-ym)/(2/5*rX + xM-xm))
  M <- as.matrix(expand.grid(a,b))
  ngrid <- tps2d(M, matr, matt)
  if (k3 == FALSE){
    plot.new()
    plot.window(1.05*range(ngrid[,1]), 1.05*range(ngrid[,2]), xaxt="n", yaxt="n", 
                xlab="", ylab="", bty="n", asp = 1)
    for (i in 1:m){
      plot.xy(xy.coords(ngrid[(1:n)+(i-1)*n,]), type = "l",
              col=grid.col, lwd=grid.lwd, lty=grid.lty)
    }
    for (i in 1:n){
      plot.xy(xy.coords(ngrid[(1:m)*n-i+1,]), type = "l",
              col=grid.col, lwd=grid.lwd, lty=grid.lty)
    }
    if(refpts==FALSE) {
      plot.xy(xy.coords(matt), type="p", pch=21, bg=pt.bg, cex=sz) 
    } else {
      plot.xy(xy.coords(matr), type="p", pch=21, bg=pt.bg, cex=sz)
    }
  }
  # added in for shape.predictor; plots in rgl window
  if(k3 == TRUE){
    ngrid <- cbind(ngrid,0)
    plot3d(ngrid, cex=0.2, aspect=FALSE, axes=FALSE, xlab="", ylab="", zlab="")
    for (i in 1:m){
      lines3d(ngrid[(1:n)+(i-1)*n,], col=grid.col, lwd=grid.lwd, lty=grid.lty)
      }
    for (i in 1:n){
      lines3d(ngrid[(1:m)*n-i+1,], col=grid.col, lwd=grid.lwd, lty=grid.lty)
      }
    if(refpts==FALSE) {
      points3d(cbind(matt,0), col=pt.bg, size=sz*10)
      } else {
        points3d(cbind(matr,0),col=pt.bg,size=sz*10)
      }
  }
}

# tps2d
#
#
tps2d <- function(M, matr, matt){
  p <- dim(matr)[1]; q <- dim(M)[1]; n1 <- p+3
  P <- matrix(NA, p, p)
  for (i in 1:p){
    for (j in 1:p){
      r2 <- sum((matr[i,] - matr[j,])^2)
      P[i,j] <- r2*log(r2)
    }
  }
  P[which(is.na(P))] <- 0
  Q <- cbind(1, matr)
  L <- rbind(cbind(P, Q), cbind(t(Q), matrix(0,3,3)))
  m2 <- rbind(matt, matrix(0, 3, 2))
  coefx <- fast.solve(L)%*%m2[,1]
  coefy <- fast.solve(L)%*%m2[,2]
  fx <- function(matr, M, coef){
    Xn <- numeric(q)
    for (i in 1:q){
      Z <- apply((matr - matrix(M[i,], p, 2, byrow=TRUE))^2, 1, sum)
      Xn[i] <- coef[p+1] + coef[p+2]*M[i,1] + coef[p+3]*M[i,2] + sum(coef[1:p]*(Z*log(Z)))
    }
    return(Xn)
  } 
  matg <- matrix(NA, q, 2)
  matg[,1] <- fx(matr, M, coefx)
  matg[,2] <- fx(matr, M, coefy)
  return(matg)
}

# tps2d3d
#
#
tps2d3d <- function(M, matr, matt, PB=TRUE){		#DCA: altered from J. Claude 2008
  p <- dim(matr)[1]; k <- dim(matr)[2]; q <- dim(M)[1]
  Pdist <- as.matrix(dist(matr))
  ifelse(k==2, P <- Pdist^2*log(Pdist^2), P <- Pdist)
  P[which(is.na(P))] <- 0
  Q <- cbind(1, matr)
  L <-rbind(cbind(P, Q), cbind(t(Q), matrix(0,k+1,k+1)))
  m2 <- rbind(matt, matrix(0, k+1, k))
  coefx <- fast.solve(L)%*%m2[,1]
  coefy <- fast.solve(L)%*%m2[,2]
  if(k==3){coefz <- fast.solve(L)%*%m2[,3]}
  fx <- function(matr, M, coef, step){
    Xn <- numeric(q)
    for (i in 1:q){
      Z <- apply((matr-matrix(M[i,], p, k, byrow=TRUE))^2, 1, sum)
      ifelse(k==2, Z1<-Z*log(Z), Z1<-sqrt(Z)); Z1[which(is.na(Z1))] <- 0
      ifelse(k==2, Xn[i] <- coef[p+1] + coef[p+2]*M[i,1] + coef[p+3]*M[i,2] + sum(coef[1:p]*Z1),
             Xn[i] <- coef[p+1] + coef[p+2]*M[i,1] + coef[p+3]*M[i,2] + coef[p+4]*M[i,3] + sum(coef[1:p]*Z1))
      if(PB==TRUE){setTxtProgressBar(pb, step + i)}
    }
    return(Xn)
    }
  matg <- matrix(NA, q, k)
  if(PB==TRUE){pb <- txtProgressBar(min = 0, max = q*k, style = 3) }
  matg[,1] <- fx(matr, M, coefx, step = 1)
  matg[,2] <- fx(matr, M, coefy, step=q)
  if(k==3){matg[,3] <- fx(matr, M, coefz, step=q*2)
  }
  if(PB==TRUE) close(pb)
  return(matg)
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

# gdf.to.df
# attempts to coerce a geomorph data frame to a data frame
# but only for relevant parts
# used in advanced.procD.lm
gdf.to.df <- function(L){
  if(!is.list(L)) stop("Missing list to convert to data frame")
  check1 <- sapply(L, class)
  match1 <- match(check1, c("numeric", "matrix", "vector", "factor"))
  Lnew <- L[!is.na(match1)]
  check2 <- sapply(Lnew, NROW)
  if(length(unique(check2)) == 1) Lnew <- Lnew else {
    check3 <- as.vector(by(check2, check2, length))
    check4 <- check2[which.max(check3)]
    Lnew <- Lnew[check2 == check4]
  }
  as.data.frame(Lnew)
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
    y.rand <-lapply(1:length(j), function(i) as.matrix(y[ind.j[[i]],]))
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
    y.rand <-lapply(1:length(j), function(i) as.matrix(y[ind.j[[i]],]))
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
CR <-function(x,gps){
  x <- center(x)
  g <- gps
  gl <- unique(g)
  ngps <- length(gl)
  gps.combo <- combn(ngps, 2)
  Xs<- lapply(1:ngps, function(j) {
    x[, g %in% gl[j]]
  })
  CR.gp <- sapply(1:ncol(gps.combo), function(j){ # no loops
    ind <- gps.combo[,j]; a <- ind[1]; b <- ind[2]
    S11 <- crossprod(Xs[[a]]); diag(S11) <- 0
    S22 <- crossprod(Xs[[b]]); diag(S22) <- 0
    S12 <- crossprod(Xs[[a]], Xs[[b]])
    sqrt(sum(S12^2)/sqrt(sum(S11^2)*sum(S22^2)))
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
  x <- center(x)
  g <- gps
  gl <- unique(g)
  ngps <- length(gl)
  gps.combo <- combn(ngps, 2)
  Xs<- lapply(1:ngps, function(j) {
    x[, g %in% gl[j]]
  })
    CR.gp <- sapply(1:ncol(gps.combo), function(j){ # no loops
      ind <- gps.combo[,j]; a <- ind[1]; b <- ind[2]
      S11 <- crossprod(Xs[[a]]); diag(S11) <- 0
      S22 <- crossprod(Xs[[b]]); diag(S22) <- 0
      S12 <- crossprod(Xs[[a]], Xs[[b]])
      sqrt(sum(S12^2)/sqrt(sum(S11^2)*sum(S22^2)))
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
  C<-fast.phy.vcv(phy)
  C<-C[rownames(x),rownames(x)]
  invC <-fast.solve(C)
  eigC <- eigen(C)
  lambda <- zapsmall(abs(Re(eigC$values)))
  if(any(lambda == 0)){
    warning("Singular phylogenetic covariance matrix. Proceed with caution")
  }
  D.mat <- fast.solve(eigC$vectors%*% diag(sqrt(abs(eigC$values))) %*% t(eigC$vectors))
  rownames(D.mat) <- colnames(D.mat) <- colnames(C)
  rownames(invC) <- colnames(invC) <- colnames(C)
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
#  xp <- Ptrans%*%x
  xp <- x
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

# GMfromShapes0
# function to read landmarks without sliders from a StereoMorph shapes object
# used in: readland.shapes
GMfromShapes0 <- function(Shapes, scaled = TRUE){ # No curves
  scaling <- Shapes$scaling
  sp.names <- dimnames(Shapes$landmarks.pixel)[[3]]
  lm.names <- dimnames(Shapes$landmarks.pixel)[[1]]
  if(is.null(scaling)) {
    landmarks <- Shapes$landmarks.pixel
    cat("\nWarning: No specimens have scaling")
    cat("\nUnscaled landmarks imported, as a result\n")
    scaled = FALSE
    } else {
      if(any(is.na(scaling))){
        sp.na <- which(is.na(scaling))
        cat("\nWarning: Some specimens have no scaling\n")
        cat(sp.names[sp.na])
        cat("\nUnscaled landmarks imported, as a result\n")
        landmarks <- Shapes$landmarks.pixel
        scaled = FALSE
      } else {
        if(scaled) landmarks <- Shapes$landmarks.scaled else 
          landmarks <- Shapes$landmarks.pixel
      }
    }
  dims <- dim(landmarks)
  n <- dims[[3]]; p <- dims[[1]]; k <- dims[[2]]
  landmarks <- lapply(1:n, function(j){
    landmarks[,,j]
  })
  names(landmarks) <- sp.names
  
  out <- list(landmarks = landmarks, fixed = 1:dims[[1]],
              sliders = NULL, curves = NULL, n = n, 
              p = p, k = k, scaled = scaled)
  class(out) <- "geomorphShapes"
  invisible(out)
}

# evenPts
# basic function for spacing out curve points via linear interpolation
# simple form of pointsAtEvenSpacing from StereoMorph 
# used in: readland.shapes and difit.curves
evenPts <- function(x, n){
  x <- as.matrix(na.omit(x))
  N <- NROW(x); p <- NCOL(x)
  if(N == 1) stop("x must be a matrix")
  if(n < 3) {
    n <- 2
    nn <- 3 # so lapply function works
  } else nn <- n
  
  if(N == 2) {
    x <- rbind(x, x[2,])
    N <- 3 # third row cut off later
  }
    xx <- x[2:N, ] - x[1:(N - 1), ]
    ds <- sqrt(rowSums(xx^2))
    cds <- c(0, cumsum(ds))
    cuts <- cumsum(rep(cds[N]/(n-1), n-1))
    targets <- lapply(1:(nn-2), function(j){
      dtar <- cuts[j]
      ll <- which.max(cds[cds < dtar])
      ul <- ll + 1
      adj <- (dtar- cds[ll])/(cds[[ul]] - cds[ll])
      x[ll,] + adj * (x[ul,] - x[ll,])
    })
    
  out <- matrix(c(x[1,], unlist(targets), x[N,]), n, p, byrow = TRUE)
  out
}

# GMfromShapes1
# function to read landmarks with sliders from a StereoMorph shapes object
# used in: readland.shapes
GMfromShapes1 <- function(Shapes, nCurvePts, continuous.curve = NULL, scaled = TRUE){ # with Curves
  out <- GMfromShapes0(Shapes, scaled = scaled)
  scaled = out$scaled
  if(scaled) curves <- Shapes$curves.scaled else 
    curves <- Shapes$curves.pixel
  sp.names <- names(out$landmarks)
  n <- out$n; p <- out$p; k <- out$k
  
  # define fixed landmarks
  fixedLM  <- out$landmarks
  
  # check for matching fixedLM and curves
  fLM.names <- names(fixedLM)
  cv.names <- names (curves)
  
  if(length(fLM.names) != length(cv.names)) {
    mismatch <- setdiff(fLM.names, cv.names)
    prob <- paste("\nThese specimens are missing either landmarks or curves:", mismatch, "\n")
    stop(prob)
  }
  
  # define curves
  curves.check <- sapply(1:n, function(j) length(curves[[j]]))
  if(length(unique(curves.check)) > 1) {
    names(curves.check) <- names(curves)
    cat("\nSpecimens have different numbers of curves\n")
    Mode <- as.numeric(names(which.max(table(curves.check))))
    cat("\nTypical number of curves:", Mode, "\n")
    cat("\nCheck these specimens (number of curves indicated)\n")
    print(curves.check[curves.check != Mode])
    stop("\nCannot proceed until this issue is resolved")
  }
    
    
  curve.n <- curves.check[1]
  curve.nms <- names(curves[[1]])
  nCurvePts <- unlist(nCurvePts)
  if(length(nCurvePts) != curve.n) stop("The number of curves and length of nCurvePts do not match")
  curves <- lapply(1:n, function(j){
    x <- curves[[j]]
    res <- lapply(1:curve.n, function(jj) evenPts(x[[jj]], nCurvePts[[jj]]))
    names(res) <- curve.nms
    res
  })
  
  # index fixed landmarks in curves (anchors)
  anchors <- lapply(1:curve.n, function(j){
    cv <- curves[[1]][[j]]
    lm <- fixedLM[[1]]
    if(!is.matrix(lm)) lm <- matrix(lm, nrow = 1)
    ends <- rbind(cv[1,], cv[nrow(cv),])
    a <- which(apply(lm, 1,function(x) identical(x, ends[1,])))
    b <- which(apply(lm, 1,function(x) identical(x, ends[2,])))
    c(a,b)
  })
  
  # determine which curve points become landmarks
  curve.refs <- list()
  
  for(i in 1:curve.n) {
    cp <- nCurvePts[i]
    ce <- c(2, cp+1) - 1
    cvr <- (1:cp)[-ce]
    curve.refs[[i]] <- cvr
  }
  
  # create indicator values for semilandmarks
  lm.curve.refs<- list()
  pp <- p
  for(i in 1:curve.n){
    cr <- curve.refs[[i]]
    if(length(cr) == 1 && cr == 0) cp <- NA else
      cp <- cr + pp - 1
    lm.curve.refs[[i]] <- cp
    if(is.null(cp)) pp <- pp else pp <- max(cp)
  }
  
  # assemble the landmarks, both fixed and semi
  landmarks <- lapply(1:n, function(j){
    x <- fixedLM[[j]]
    cv <- curves[[j]]
    for(i in 1:curve.n) x <- rbind(x, cv[[i]][curve.refs[[i]],])
    rownames(x)[(p+1):(nrow(x))] <- paste("curveLM", (p+1):(nrow(x)), sep=".")
    x
  })
  names(landmarks) <- sp.names
  
  # make an index in case any curves are continuous
  cc <- rep(0, curve.n)
  if(!is.null(continuous.curve)){
    continuous.curve <- unlist(continuous.curve)
    cc[continuous.curve] = 1
  } 
  
  # create a curves matrix
  curves.mat.parts <- lapply(1:curve.n, function(j){
    lcr <- lm.curve.refs[[j]]
    if(all(is.na(lcr))) {
      mat.part <- NULL
      } else {
        anch <- anchors[[j]]
        strp <- c(anch[[1]], lcr, anch[[2]])
        if(cc[j] == 1) strp <- c(strp, strp[2])
        n.mp <- length(strp) - 2
        mat.part <- matrix(NA, n.mp, 3)
        for(i in 1:n.mp) mat.part[i,] <- strp[i:(i+2)]
      }
    mat.part
  })
  
  curves.mat <- curves.mat.parts[[1]]
  if(curve.n >= 2){
    for(i in 2:curve.n) if(!is.null(curves.mat.parts[[i]]))
      curves.mat <- rbind(curves.mat, curves.mat.parts[[i]])
  }
  single.anchors <- sapply(lapply(anchors, unique), length) == 1
  fixed <- out$fixed
  sliders <- (1:nrow(landmarks[[1]]))[-fixed]
  makeSemi <- as.numeric(single.anchors) * cc
  for(i in 1:length(makeSemi)){
    if(makeSemi[i] > 0) {
      targ <- unique(anchors[[i]])
      fixed <- fixed[-targ]
      sliders <- c(targ, sliders)
    }
  }
  
  curve.mat.nms <- rownames(landmarks[[1]])[curves.mat[,2]]
  rownames(curves.mat) <- curve.mat.nms
  
  # Assemble the output
  out$landmarks <- landmarks
  out$fixed <- fixed
  out$sliders <- sliders
  out$curves <- curves.mat
  out$p <- length(fixed) + length(sliders)
  
  invisible(out)
}




#####-----------------------------------------------------------------------------------

### geomorph-specific logicals

is.gpagen <- function(x) inherits(x, "gpagen")
is.phylo <- function(x) inherits(x, "phylo")
is.geomorph.data.frame <- function(x) inherits(x, "geomorph.data.frame")

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


### All functions below perform in a constrained way the same as functions in ape and geiger:
### -----------------------------------------------------------------------------------------


# sim.char.BM
# An alternative to sim.char for BM
# compared to sim.char, uses better memoization


# sim.char has a similar function to this but it is called in every simulation
# and defers to C for help.  This is done once only here. (Produces a projection matrix)
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

# A function to extract eigen values, as per mvrnorm
# But does not incorporate other traps and options, as in mvnorm.
Sig.eigen <- function(Sig, tol = 1e-06){
  E <- eigen(Sig)
  ev <- E$values
  if(any(Im(ev)==0)) ev<-Re(ev)
  if (!all(ev >= -tol * abs(ev[1L]))) {
    k <- which(ev >= -tol * abs(ev[1L]))
    E$values <- ev[k]
    E$vectors <- E$vectors[,k]
  } else k <- NCOL(Sig)
  E$scaled.vectors <- as.matrix(E$vectors)[, 1:k] %*% diag(sqrt(E$values), k)
  rownames(E$vectors) <- rownames(E$scaled.vectors) <- 
    names(E$values) <- rownames(Sig)
  E$p <- length(ev)
  E
}

# For a single permutation, simulate traits
# Requires estabished scaled eigenvectors (E)
# and a projection matrix (M),  sim.char re-estimates
# E and M in every permutation
fast.sim.BM <- function(E, M, R, root = 1) {
  N <- ncol(M)
  p <- E$p
  v <- E$scaled.vectors
  X <- v %*% t(matrix(R, N, p))
  res <- tcrossprod(M, X) + root
  rownames(res) <- rownames(M)
  res
}

sim.set.up <- function(E, M, nsim, seed = NULL){
  if(is.null(seed)) seed = nsim else
    if(seed == "random") seed = sample(1:nsim, 1) else
      if(!is.numeric(seed)) seed = nsim 
      N <- ncol(M)
      p <- E$p
      n <- N * p
      NN <- n * nsim
      set.seed(seed)
      R <- rnorm(NN)
      rm(.Random.seed, envir=globalenv())
      if(nsim > 1) {
        r.start <- seq(1,(NN-n+1), n)
        r.stop <- seq(n, NN, n)
      } else {
        r.start <- 1
        r.stop <- NN
      }
      
      lapply(1:nsim, function(j) R[r.start[j] : r.stop[j]])
}

# Put everything together to make a list of results
sim.char.BM <- function(phy, par, nsim = 1, root = 1, seed = NULL){
  M <- phy.sim.mat(phy)
  E <- Sig.eigen(par)
  sim.set <- sim.set.up(E, M, nsim, seed = seed)
  sim.args <- list(E=E, M=M, R=0, root=root)
  if(nsim == 1) {
    sim.args$R <- sim.set[[1]]
    res <- do.call(fast.sim.BM, sim.args)
  } else 
    res <- lapply(1:nsim, function(j){
      sim.args$R <- sim.set[[j]]
      do.call(fast.sim.BM, sim.args)
    })
  
  res
}


# makePD
# Alternative to nearPD in Matrix

makePD <- function (x) {
  eig.tol <- conv.tol <- 1e-06
  posd.tol <- 1e-08
  maxit <- 50
  n <- ncol(x)
  X <-x
  D <- matrix(0, n, n)
  iter <- 0
  converged <- FALSE
  conv <- Inf
  while (iter < maxit && !converged) {
    Y <- X
    R <- Y - D
    e <- eigen(R, symmetric = TRUE)
    Q <- e$vectors
    d <- e$values
    p <- d > eig.tol * d[1]
    if (!any(p)) 
      stop("Matrix seems negative semi-definite")
    Q <- Q[, p, drop = FALSE]
    X <- tcrossprod(Q * rep(d[p], each = nrow(Q)), Q)
    D <- X - R
    conv <- norm(Y - X, type = "I")/norm(Y, type = "I")
    iter <- iter + 1
    converged <- (conv <= conv.tol)
  }
  
  e <- eigen(X, symmetric = TRUE)
  d <- e$values
  Eps <- posd.tol * abs(d[1])
  if (d[n] < Eps) {
    d[d < Eps] <- Eps
    Q <- e$vectors
    o.diag <- diag(X)
    X <- Q %*% (d * t(Q))
    D <- sqrt(pmax(Eps, o.diag)/diag(X))
    X <- D * X * rep(D, each = n)
  }
  X
}

# ape replacement functions below ----------------------------------------------------

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
    preps <- pic.prep(phy.j, NROW(Y), NCOL(Y))
    preps$x <- Y
    preps$tip.label <- phy$tip.label
    out <- do.call(ace.pics, preps)
    out[n + 1,]
  }))
  
  if(length(out) == (n-1)) out <- t(out)
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
