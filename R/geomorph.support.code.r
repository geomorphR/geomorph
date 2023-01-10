#' @name geomorph-package
#' @docType package
#' @aliases geomorph
#' @title Geometric morphometric analyses for 2D/3D data
#' @author Dean C. Adams, Michael Collyer, Antigoni Kaliontzopoulou, and Erica Baken
#' @description Functions in this package allow one to read, manipulate, and digitize landmark data; generate shape
#'  variables via Procrustes analysis for points, curves and surface data, perform statistical analyses
#'  of shape variation and covariation, and provide graphical depictions of shapes and patterns of
#'  shape variation.
#'
#' @importFrom jpeg readJPEG
#' @importFrom ape multi2di.phylo
#' @importFrom ape root.phylo
#' @importFrom ape collapse.singles
#' @importFrom graphics abline arrows hist identify layout lines locator mtext par plot.new plot.window plot.xy points rasterImage segments text title
#' @importFrom stats na.omit .lm.fit anova as.dist as.formula cmdscale coef cor cov dist formula kmeans lm lm.fit model.matrix optimise pnorm prcomp predict quantile rnorm sd spline terms var
#' @importFrom utils combn object.size read.csv read.table setTxtProgressBar txtProgressBar write.csv write.table
#' @import RRPP
#' @import parallel
#' @import rgl
#' @import grDevices
#' @import ggplot2
#' @import Matrix
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

#' Head shape and phylogenetic relationships for several Plethodon salamander species
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

#' Dorsal head shape data of lizards
#'
#' @name lizards
#' @docType data
#' @author Antigoni Kaliontzopoulou
#' @keywords datasets
#' @description Superimposed landmark data of the dorsal view of lizard heads
#' @details
#' Dataset includes superimposed landmarks (coords), centroid size (cs), an index of individuals (ind) and digitizing repetitions (rep), and a table of symmetrical matching
#' landmarks (lm.pairs). The object is a \code{\link{geomorph.data.frame}}.
#' The dataset corresponds to the data for population "b" from Lazic et al. 2015.
#' @references Lazic, M., Carretero, M.A., Crnobrnja-Isailovic, J. & Kaliontzopoulou, A. 2015. Effects of
#' environmental disturbance on phenotypic variation: an integrated assessment of canalization, developmental 
#' stability, modularity and allometry in lizard head shape. American Naturalist 185: 44-58.
NULL

#' Estimate mean shape for a set of aligned specimens
#'
#' Estimate the mean shape for a set of aligned specimens
#'
#' The function estimates the average landmark coordinates for a set of aligned specimens. 
#' Three different methods are available for missing data (see Arguments and Examples).
#'  
#'  One can then use the generic function \code{\link{plot}} to produce a numbered plot of landmark 
#'  positions and potentially add links, in order to review landmark positions
#'
#' @param A Either a list (length n, p x k), A 3D array (p x k x n), or a matrix (n x pk) containing GPA-aligned coordinates for a set of specimens
#' @param na.action An index for how to treat missing values: 1 = stop analysis; 2 = return NA for coordinates
#' with missing values for any specimen; 3 = attempt to calculate means for coordinates for all non-missing values.
#' @keywords utilities
#' @export
#' @author Antigoni Kaliontzopoulou & Michael Collyer
#' @examples
#' data(plethodon)
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
#' A <- Y.gpa$coords
#' A[[1]] <- NA # make a missing value, just for example
#'
#' mshape(Y.gpa$coords)   # mean (consensus) configuration
#' # mshape(A, na.action = 1) # will return an error
#' mshape(A, na.action = 2) # returns NA in spot of missing value
#' mshape(A, na.action = 3) # finds mean values from all possible values
#' 
mshape <- function(A, na.action = 1){
  
  na.check <- na.action %in% 1:3
  if(!na.check) na.action <- 1
  
  if(!inherits(A, c("list", "array", "matrix")))
    stop("Data are neither a list nor array. mshape is not possible with data in the format used.\n",
         call. = FALSE)
  
  if(is.array(A)) {
    dims <- dim(A)
    if(length(dims) == 3) {
      p <- dims[[1]]
      k <- dims[[2]]
      n <- dims[[3]]
      L <- lapply(1:n, function(j) as.matrix(A[,,j]))
    }
    
    if(length(dims) == 2) {
      if(dims[[2]] == 2 || dims[[2]] == 3) {
        L <- list(A) 
        n <- 1
        p <- dims[[1]]
        k <- dims[[2]]
        } else {
          n <- dims[[1]]
          p <- dims[[2]]
          k <- 1
          warning(paste("\nWarning: It appears that data are in a matrix with specimens as rows.\n",
        "\nMeans are found for each column of the matrix.\n\n"), immediate. = TRUE)
        L <- lapply(1:n, function(j) matrix(A[j,], 1, ))
      }
    }
  }
    
    if(is.list(A)) {
      matrix.check <- sapply(A, is.matrix)
      if(any(!matrix.check)) stop("At least one specimen is not a data matrix.\n", call. = FALSE)
      dims <- dim(A[[1]])
      A.dims <- sapply(A, dim)
      A.unique <- apply(A.dims, 1, unique)
      if(!identical(dims, A.unique))
        stop("Not all specimens have the same number of landmarks or landmarks in the same dimensions.\n", 
             call. = FALSE)
      L <- A
      if(length(L) > 1) {
        p <- dims[[1]]
        k <- dims[[2]]
        n <- length(L)
      } else {
        n <- dims[[1]]
        p <- dims[[2]]
        k <- 1
      }
    }
    
    if(na.action == 1) {
      
      if(any(is.na(unlist(L))))
        stop("Data matrix contains missing values.\n 
             Estimate these first (see 'estimate.missing') or chamge the na.action (see Arguments).\n",
             call. = FALSE)
      
      res <- if(length(L) == 1) L else Reduce("+", L)/n
    }
    
    if(na.action == 2) {
      
      if(any(is.na(unlist(L)))) {
        warning(
          paste(
            "\nThis is not an error!  It is a friendly warning.\n",
            "\nMissing values detected.\n",
            "\nNA is returned for any coordinate where missing values were found.",
            "\nYou can estimate missing values (see 'estimate.missing')",
            "\nor change the na.action (see Arguments) to find the mean of just the remaining values.\n",
            "\nUse options(warn = -1) to turn off these warnings. \n\n", sep = " "),
          noBreaks. = TRUE, call. = FALSE, immediate. = TRUE) 
      }
      res <- if(length(L) == 1) L else Reduce("+", L)/n
      }
  
    
    if(na.action == 3) {
      
      if(any(is.na(unlist(L)))) {
        
        warning(
          paste(
            "\nThis is not an error!  It is a friendly warning.\n",
            "\nMissing values detected.\n",
            "\nMeans are calculated only for values that are found.",
            "\nYou can estimate missing values (see 'estimate.missing')",
            "\nor change the na.action (see Arguments) to return NA for coordinates that have missing values.\n",
            "\nUse options(warn = -1) to turn off these warnings. \n\n", sep = " "),
          noBreaks. = TRUE, call. = FALSE, immediate. = TRUE) 
           }
      
      mmean <- function(L) { 
        kp <- k * p
        U <- unlist(L)
        u <- length(U)
        starts <- seq.int(1, u, kp)
        M <- array(NA, p*k)
        for(i in 1:kp) {
          pts <- starts -1 + i
          M[i] <- mean(na.omit(U[pts]))
        }
        
        if(k == 1) M <- matrix(M, 1, p) else M <- matrix(M, p, k)
        M
      }
      
      res <- if(length(L) == 1) L else mmean(L)
    }
    
  if(inherits(res, "list")) res <- res[[1]]
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
    lm <- suppressWarnings(lapply(templm, as.numeric))
    p <- length(lm)
    k <- sapply(lm, length)
    if(length(unique(k)) == 1) k <- unique(k)
    scale <- as.numeric(unlist(strsplit(temp[scl], "SCALE="))[2])
    id <- unlist(strsplit(temp[idl], "ID="))[2]
    image <- unlist(strsplit(temp[iml], "IMAGE="))[2]
    image <- sub(".jpg", "", image, ignore.case)
    image <- sub(".tiff", "", image, ignore.case)
    image <- sub(".tif", "", image, ignore.case)
    image <- sub(".bmp", "", image, ignore.case)
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
tangents = function(s, x, scaled=FALSE){ # s = curves, x = landmarks
  if(nrow(s) > 1) { ts <- x[s[, 3], ] - x[s[, 1], ]} else { ts <- x[s[3]] - x[s[1]]}
  if(scaled==TRUE) {
    if(nrow(as.matrix(ts)) > 1) {
      ts.scale = sqrt(rowSums(ts^2))
    } else {ts.scale = sqrt(sum(ts^2))} 
    ts <- ts/ts.scale
  }
  y <- matrix(0, nrow(x), ncol(x))
  y[s[,2],] <- ts
  y
}

# nearest
# finds nearest points on surfaces for sliding semilandmakrs
# used in all functions associated with pPga.wCurves
nearest <- function(X, m, k = 4, n) {
  a <- X[m,]
  b <- vapply(1:n, function (j) sum((a-X[j,])^2), 1)
  match(sort(b)[2:(k + 1)], b)
}

# getU_s
# calculates U matrix for sliding semilandmarks
# creates sparseMatrix object, so must
# be adapted for Matrix calculations
getU_s <- function(y, tans = NULL, surf = NULL, 
                   BEp = NULL, 
                   m, k){
  
  if(!is.null(tans) && is.null(surf)) {
    U <- matrix(0, m * k, m)
    tag <- FALSE
  } else if(is.null(tans) && !is.null(surf)) {
    U <- matrix(0, m * k, 2 * m)
    tag <- FALSE
  } else {
    U <- matrix(0, m * k, 3 * m)
    tag <- TRUE
  }
  
  if(is.null(BEp)) BEp <- 1:m
  
  if(!is.null(tans)){
    tx <- tans[,1][BEp]
    ty <- tans[,2][BEp]
    tz <- if(k == 3) tans[,3 ][BEp] else NULL
    
    U[1:m, 1:m] <- diag(tx)
    U[(m+1):(2*m), 1:m] <- diag(ty)
    if(k == 3) U[(2*m+1):(3*m), 1:m] <- diag(tz)
    
  }
  
  if(!is.null(surf)) {
    tp <- if(tag) m else 0
    
    PC <- getSurfPCs(y, surf)
    p1x <- PC$p1x[BEp]
    p1y <- PC$p1y[BEp]
    p1z <- PC$p1z[BEp]
    p2x <- PC$p2x[BEp]
    p2y <- PC$p2y[BEp]
    p2z <- PC$p2z[BEp]
    
    U[1:m, (tp+1):(tp+m)] <- diag(p1x)
    U[(m+1):(2*m), (tp+1):(tp+m)] <- diag(p1y)
    if(k == 3) U[(2*m+1):(3*m), (tp+1):(tp+m)] <- diag(p1z)
    U[1:m, (tp+m+1):(tp+2*m)] <- diag(p2x)
    U[(m+1):(2*m), (tp+m+1):(tp+2*m)] <- diag(p2y)
    if(k == 3) U[(2*m+1):(3*m), (tp+m+1):(tp+2*m)] <- diag(p2z)
    
  }
  
  keep <- which(colSums(U^2) != 0)
  Matrix(U[, keep], sparse = TRUE)
  
}

# Ltemplate
# calculates inverse of bending energy matrix
# used in any function that calculates bending energy
# used in BE.slide
Ltemplate <- function(Mr, Mt=NULL){
  p <-nrow(Mr); k <- ncol(Mr)
  if(!is.null(Mt)) P <- as.matrix(dist(Mr-Mt)) else P <- as.matrix(dist(Mr))
  if(k==2) {P <-P^2*log(P); P[is.na(P)] <- 0}
  Q <- cbind(1, Mr)
  L <- matrix(0, p + k + 1, p + k + 1)
  L[1:p, 1:p] <- P
  L[1:p, (p + 1):(p + k + 1)] <- Q
  L[(p + 1):(p + k + 1), 1:p] <- t(Q)
  L <- forceSymmetric(L)
  Linv <- try(solve(L), silent = TRUE)
  if(inherits(Linv, "try-error")) {
    cat("\nSingular BE matrix; using generalized inverse")
    Linv <- -fast.ginv(as.matrix(L))
  }
  as.matrix(Linv)[1:p, 1:p]
}

# fast.solveSym
# same as fast.solve, but specifically for Ltemplate
# attempts to coerce solve
fast.solveSym <- function(x) { 
  res <- try(solve(x), silent = TRUE)
  if(inherits(res, "try-error")) {
    res <- try(as.matrix(solve(forceSymmetric(x))), silent = TRUE)
  }
  if(inherits(res, "try-error")) {
    res <- fast.ginv(x)
  }
  return(res)
}

# sparse.solve
# same as fast.solve, but specifically for sparse symmetric matrices
# used in any function requiring a generalized inverse of a known
# sparse symmetric matrix

sparse.solve <- function(X){
  keepx <- which(rowSums(X^2)!= 0)
  keepy <- which(rowSums(X^2)!= 0)
  Y <- X[keepx, keepy]
  Xn <- X
  Xn[keepx, keepy] <- fast.solveSym(Y)
  Xn
}

# fast.solveSym
# same as fast.solve, but specifically for Ltemplate
# attempts to coerce solve
fast.solveSym <- function(x) { 
  res <- try(solve(x), silent = TRUE)
  if(inherits(res, "try-error")) {
    res <- try(as.matrix(solve(forceSymmetric(x))), silent = TRUE)
  }
  if(inherits(res, "try-error")) {
    res <- fast.ginv(x)
  }
  return(res)
}

# For approximating Linv
# Should avoid ginv 

LtemplateApprox <- function(ref){ # only for sliding functions
  p <-nrow(ref); k <- ncol(ref)
  P <- as.matrix(dist(ref))
  if(k==2) {P <-P^2*log(P); P[is.na(P)] <- 0}
  Q <- cbind(1, ref)
  L <- matrix(0, p + k + 1, p + k + 1)
  L[1:p, 1:p] <- P
  L[1:p, (p + 1):(p + k + 1)] <- Q
  L[(p + 1):(p + k + 1), 1:p] <- t(Q)
  diag(L) <- 0.0001 * sd(L)
  L <- forceSymmetric(L)
  Linv <- try(solve(L), silent = TRUE)
  if(inherits(Linv, "try-error")) {
    cat("\nSingular BE matrix; using generalized inverse")
    Linv <- -fast.ginv(as.matrix(L))
  }
  as.matrix(Linv)[1:p, 1:p]
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
pGpa <- function(Y, PrinAxes = FALSE, Proj = FALSE, max.iter = 10, tol){
  max.iter <- min(max.iter, 4)
  
  iter <- 0
  pb <- txtProgressBar(min = 0, max = max.iter, initial = 0, style=3)
  setTxtProgressBar(pb,iter)
  n <- length(Y); dims <- dim(Y[[1]]); p <- dims[1]; k <- dims[2]
  Yc <- Map(function(y) center.scale(y), Y)
  CS <- sapply(Yc,"[[","CS")
  Ya <- lapply(Yc,"[[","coords")
  Ya <- apply.pPsup(Ya[[1]], Ya)
  M <- Reduce("+",Ya)/n
  Q <- ss <- n*(1-sum(M^2))
  M <- cs.scale(M)
  iter <- 1
  setTxtProgressBar(pb,iter)
  while(Q > tol){
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


.pGpa <- function(Y, PrinAxes = FALSE, Proj = FALSE, 
                  Parallel = TRUE, max.iter = 10, tol){
  
  max.iter <- min(max.iter, 4)
  
  ParCores <- NULL
  if (is.numeric(Parallel)) {
    ParCores <- Parallel
    Parallel <- TRUE
  }
  if (Parallel && is.null(ParCores)) {
    ParCores <- detectCores() - 1
  }
  
  if (is.numeric(ParCores)) {
    if(ParCores > detectCores() - 1) ParCores <- detectCores() - 1
  } 
  
  n <- length(Y); p <- nrow(Y[[1]]); k <- ncol(Y[[1]])
  Yc <- Map(function(y) center.scale(y), Y)
  CS <- sapply(Yc,"[[","CS")
  Ya <- lapply(Yc,"[[","coords")
  Ya <- apply.pPsup(Ya[[1]], Ya)
  M <- Reduce("+", Ya)/n
  Q <- ss <- n*(1-sum(M^2))
  M <- cs.scale(M)
  iter <- 1

  if(Parallel) {
    Unix <- .Platform$OS.type == "unix" 
    
    if(Unix) {
      
      while(Q > tol){
        iter <- iter+1
        Y <- Ya
        Ya <- mclapply(1:n, function(j) {
        y <- Y[[j]]
        MY <- crossprod(M, y)
        sv <- La.svd(MY, k, k)
        u <- sv$u; u[,k] <- u[,k]*  determinant(MY)$sign
        tcrossprod(y, u %*% sv$vt)
        }, mc.cores = ParCores)
        
        M <- Reduce("+",Ya)/n
        ss2 <- n*(1-sum(M^2))
        Q <- abs(ss-ss2)
        ss <- ss2
        M <- cs.scale(M)
        if(iter > max.iter) break
      }
      
    } else {
      
      cl <- makeCluster(ParCores)
      
      while(Q > tol){
        iter <- iter+1
        Y <- Ya
        Ya <- parLapply(cl = cl, 1:n, function(j) {
          y <- Y[[j]]
          MY <- crossprod(M, y)
          sv <- La.svd(MY, k, k)
          u <- sv$u; u[,k] <- u[,k]*  determinant(MY)$sign
          tcrossprod(y, u %*% sv$vt)
        })
        
          M <- Reduce("+",Ya)/n
          ss2 <- n*(1-sum(M^2))
          Q <- abs(ss-ss2)
          ss <- ss2
          M <- cs.scale(M)
          if(iter > max.iter) break
      }
      stopCluster(cl)
      }
    } else {
    
      while(Q > tol){
        iter <- iter+1
        Ya <- apply.pPsup(M, Ya)
        M <- Reduce("+",Ya)/n
        ss2 <- n*(1-sum(M^2))
        Q <- abs(ss-ss2)
        ss <- ss2
        M <- cs.scale(M)
        if(iter > max.iter) break
      }
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
  p <- nrow(y)
  pc.match <- 1:p
  pc.match[-surf] <- NA
  
  d <- as.matrix(dist(y))
  nearpts <- lapply(1:p, function(j) {
    nn <- pc.match[j]
    if(!is.na(nn)) order(d[,nn])[1:5]
  })
  
  pc.dir <- lapply(1:p, function(.) matrix(0, k, k)) 
  
  take <- which(vapply(nearpts, is.null, numeric(1)) == 0)
  
  pc.take <- lapply(1:length(take), function(j) {
    take.j <- nearpts[[take[j]]]
    x <- center(y[take.j, ])
    pc <- La.svd(x, nu=0)$vt
    s=sign(diag(crossprod(V,pc)))
    pc*s
  })
  
  pc.dir[take] <- pc.take

  p1x <- vapply(pc.dir, function(x) x[1, 1], 1)
  p1y <- vapply(pc.dir, function(x) x[1, 2], 1)
  p2x <- vapply(pc.dir, function(x) x[2, 1], 1)
  p2y <- vapply(pc.dir, function(x) x[2, 2], 1)
  if(k==3) {
    p1z <- vapply(pc.dir, function(x) x[1, 3], 1)
    p2z <- vapply(pc.dir, function(x) x[2, 3], 1)
  } else
  {p1z <- NULL; p2z <- NULL}
  
  list(p1x=p1x,p1y=p1y, p2x=p2x, p2y=p2y, p1z=p1z, p2z=p2z)
  }


# semilandmarks.slide.BE
# slides landmarks along tangents of curves and PC planes of surfaces using bending energy
# used in pGpa.wSliders

semilandmarks.slide.BE <- function(y, tans = NULL,
                                   surf = NULL, ref, Lk, appBE, BEp){
  yc <- y - ref
  p <- nrow(yc)
  k <- ncol(yc)
  
  m <- if(appBE) length(BEp) else p
  if(m == p) appBE <- FALSE
  if(!appBE) BEp <- 1:p
  yvec <- as.vector(yc[BEp,])
  U <- getU_s(y, tans, surf, BEp, m, k)
  tULk <- crossprod(U, Lk) 
  mid.part <- forceSymmetric(tULk %*% U) 
  tULky <- tULk %*% yvec
  tULk <- NULL # clear space
  res <- matrix(U %*% solve(mid.part, tULky), m, k)
  mid.part <- NULL # clear space
  if(appBE) {
    temp <- matrix(0, p, k)
    temp[BEp, ] <- res
  } else temp <- res
  y - temp  
  }

# semilandmarks.slide.procD
# slides landmarks along tangents of curves and within tangent planes on surfaces using minimized ProcD
# used in pGpa.wSliders

semilandmarks.slide.ProcD <- function(y, tans = NULL, surf = NULL, ref) {
  yc <- y - ref
  p <- nrow(yc)
  k <- ncol(yc)
  
  U <- as.matrix(getU_s(y, tans, surf = surf, BEp = NULL, p, k))
  yvec <- as.vector(yc)
  res <- matrix(U %*% crossprod(U, yvec), p, k)
  
  y - res  
  }

# BE.slide
# performs sliding iterations using bending energy
# used in pGpa.wSliders
BE.slide <- function(curves = NULL, surf = NULL, Ya, ref, appBE = TRUE, 
                     sen, max.iter=10, tol){
  n <- length(Ya); p <- nrow(Ya[[1]]); k <- ncol(Ya[[1]])
  
  iter <- 1 # from initial rotation of Ya
  pb <- txtProgressBar(min = 0, max = max.iter, initial = 0, style=3)
  slid0 <- Ya
  Q <- ss0 <- sum(Reduce("+",Ya)^2)/n
  
  if(!is.numeric(sen)) sen <- 0.5
  if(sen < 0.1) sen <- 0.1
  if(sen >= 1) appBE <- FALSE
  
  if(appBE) {
    BEp <- c(if(!is.null(curves)) unique(as.vector(curves)), 
             if(!is.null(surf)) surf)
    Up <- round(seq(1, p, length.out = p * sen))
    BEp <- sort(unique(c(BEp, Up)))
  }
  
  m <- if(appBE) length(BEp) else p
  
  setTxtProgressBar(pb,iter)
  
  doTans <- !is.null(curves)
  
  if(doTans) {
    doSlide <- function(j)
      semilandmarks.slide.BE(slid0[[j]], tans[[j]],
                       surf, ref, Lk, appBE, BEp)
  } else {
    doSlide <- function(j)
      semilandmarks.slide.BE(slid0[[j]], tans = NULL, 
                       surf, ref, Lk, appBE, BEp)
    }
  
  while(Q > tol){
    iter <- iter+1
    
    if(doTans) {
      tans <- Map(function(y) tangents(curves, y, scaled=TRUE), 
                    slid0)
    } else tans <- NULL
    
    
    L <- if(appBE) LtemplateApprox(ref[BEp,]) else Ltemplate(ref)
    Lk <- kronecker(diag(k), L)
    
    slid <- lapply(1:n, doSlide)
    
    ss <- sum(Reduce("+",slid)^2)/n
    slid0 <- apply.pPsup(ref,slid)
    slid <- NULL
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

# BE.slidePP
# performs sliding iterations using bending energy
# same as BE.slide, but with parallel processing
# used in pGpa.wSliders

BE.slidePP <- function(curves = NULL, surf = NULL, Ya, ref, max.iter=10, 
                       appBE = TRUE, sen, ParCores = TRUE, tol){ 
  
  if(is.logical(ParCores)) no_cores <- detectCores() - 1 else
    no_cores <- min(detectCores() - 1, ParCores)
  Unix <- .Platform$OS.type == "unix"
  
  n <- length(Ya); p <- nrow(Ya[[1]]); k <- ncol(Ya[[1]])
  iter <- 1 # from initial rotation of Ya
  slid0 <- Ya
  Q <- ss0 <- sum(Reduce("+",Ya)^2)/n
  
  if(!is.numeric(sen)) sen <- 0.5
  if(sen < 0.1) sen <- 0.1
  if(sen >= 1) appBE <- FALSE
  
  if(appBE) {
    BEp <- c(if(!is.null(curves)) unique(as.vector(curves)), 
             if(!is.null(surf)) surf)
    Up <- round(seq(1, p, length.out = p * sen))
    BEp <- sort(unique(c(BEp, Up)))
  }
  
  m <- if(appBE) length(BEp) else p
  
  doTans <- !is.null(curves)
  
  if(doTans) doSlide <- function(j)
    semilandmarks.slide.BE(slid0[[j]], tans[[j]], 
                     surf, ref, Lk, appBE, BEp) else
      doSlide <- function(j)
        semilandmarks.slide.BE(slid0[[j]], tans = NULL, 
                               surf, ref, Lk, appBE, BEp)
  

  while(Q > tol){
    iter <- iter+1
    
    if(doTans) {
      tans <- Map(function(y) tangents(curves, y, scaled=TRUE), 
                    slid0)
    } else tans <- NULL
    
    L <- if(appBE) LtemplateApprox(ref[BEp,]) else Ltemplate(ref)
    Lk <- kronecker(diag(k), L)
    
    if(Unix) slid <- mclapply( 1:n, doSlide, mc.cores = no_cores) else {
        
        cl <- makeCluster(no_cores)
        clusterExport(cl, c("doSlide"),
                      envir=environment())
        
        slid <- parLapply(cl = cl, 1:n, doSlide)
        
        }
    
      ss <- sum(Reduce("+",slid)^2)/n
      slid0 <- apply.pPsup(ref,slid)
      slid <- NULL
      ref = cs.scale(Reduce("+", slid0)/n)
      Q <- abs(ss0-ss)
      ss0 <- ss
      if(iter >= max.iter) break
  }
  if(!Unix) stopCluster(cl)
    
  list(coords=slid0, consensus=ref, iter=iter+1, Q=Q)
}

# .BE.slide
# same as BE.slide, but without progress bar option
# used in pGpa.wSliders
.BE.slide <- function(curves, surf, Ya, ref, appBE = TRUE, 
                      sen, max.iter=10, tol){
  n <- length(Ya); p <- nrow(Ya[[1]]); k <- ncol(Ya[[1]])
  
  iter <- 1 # from initial rotation of Ya
  slid0 <- Ya
  Q <- ss0 <- sum(Reduce("+",Ya)^2)/n
  
  if(!is.numeric(sen)) sen <- 0.5
  if(sen < 0.1) sen <- 0.1
  if(sen >= 1) appBE <- FALSE
  
  if(appBE) {
    BEp <- c(if(!is.null(curves)) unique(as.vector(curves)), 
             if(!is.null(surf)) surf)
    Up <- round(seq(1, p, length.out = p * sen))
    BEp <- sort(unique(c(BEp, Up)))
  }
  
  m <- if(appBE) length(BEp) else p
  
  doTans <- !is.null(curves)
  
  if(doTans) {
    doSlide <- function(j)
      semilandmarks.slide.BE(slid0[[j]], tans[[j]],
                       surf, ref, Lk, appBE, BEp)
  } else {
    doSlide <- function(j)
      semilandmarks.slide.BE(slid0[[j]], tans = NULL, 
                       surf, ref, Lk, appBE, BEp)
  }
  

  while(Q > tol){
    iter <- iter+1
    
    if(doTans) {
      tans <- Map(function(y) tangents(curves, y, scaled=TRUE), 
                  slid0)
    } else tans <- NULL
    
    L <- if(appBE) LtemplateApprox(ref[BEp,]) else Ltemplate(ref)
    Lk <- kronecker(diag(k), L)
    
    slid <- lapply(1:n, doSlide)
    
    ss <- sum(Reduce("+",slid)^2)/n
    slid0 <- apply.pPsup(ref,slid)
    slid <- NULL
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
procD.slide <- function(curves, surf, Ya, ref, max.iter=10, tol){
  n <- length(Ya); p <- nrow(Ya[[1]]); k <- ncol(Ya[[1]])
  
  iter <- 1 # from initial rotation of Ya
  pb <- txtProgressBar(min = 0, max = max.iter, initial = 0, style=3)
  slid0 <- Ya
  Q <- ss0 <- sum(Reduce("+",Ya)^2)/n

  setTxtProgressBar(pb,iter)
  
  doTans <- !is.null(curves)
  
  if(doTans) {
    doSlide <- function(j) {
      semilandmarks.slide.ProcD(slid0[[j]], tans[[j]],
                                surf, ref)
    }
      
  } else {
    doSlide <- function(j) {
      semilandmarks.slide.ProcD(slid0[[j]], tans = NULL, 
                                surf, ref)
    }
  }
  
  while(Q > tol){
    iter <- iter+1
    
    if(doTans) {
      tans <- Map(function(y) tangents(curves, y, scaled=TRUE), 
                  slid0)
    } else tans <- NULL
    
    slid <- lapply(1:n, doSlide)
    
    ss <- sum(Reduce("+",slid)^2)/n
    slid0 <- apply.pPsup(ref,slid)
    slid <- NULL
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


# procD.slidePP
# same as procD.slide, but without progress bar option
# Same as procD.slide but with parallel processing
# used in pGpa.wSliders
procD.slidePP <- function(curves, surf, Ya, ref, max.iter=10, 
                          tol, ParCores = TRUE){
  
  if(is.logical(ParCores)) no_cores <- detectCores() - 1 else
    no_cores <- min(detectCores() - 1, ParCores)
  Unix <- .Platform$OS.type == "unix"
  
  n <- length(Ya); p <- nrow(Ya[[1]]); k <- ncol(Ya[[1]])
  iter <- 1 # from initial rotation of Ya
  slid0 <- Ya
  Q <- ss0 <- sum(Reduce("+",Ya)^2)/n

  
  doTans <- !is.null(curves)
  
  if(doTans) doSlide <- function(j)
    semilandmarks.slide.ProcD(slid0[[j]], tans[[j]], surf, ref) else
      doSlide <- function(j) semilandmarks.slide.ProcD(slid0[[j]], 
                                                       tans = NULL, surf)
  
  while(Q > tol){
    iter <- iter+1
    
    if(doTans) {
      tans <- Map(function(y) tangents(curves, y, scaled=TRUE), 
                  slid0)
    } else tans <- NULL
    
    if(Unix) slid <- mclapply( 1:n, doSlide, mc.cores = no_cores) else {
      
      cl <- makeCluster(no_cores)
      clusterExport(cl, c("doSlide"),
                    envir=environment())
      
      slid <- parLapply(cl = cl, 1:n, doSlide)
      
    }
    
    ss <- sum(Reduce("+",slid)^2)/n
    slid0 <- apply.pPsup(ref,slid)
    slid <- NULL
    ref = cs.scale(Reduce("+", slid0)/n)
    Q <- abs(ss0-ss)
    ss0 <- ss
    if(iter >= max.iter) break
  }
  if(!Unix) stopCluster(cl)
  
  list(coords=slid0, consensus=ref, iter=iter+1, Q=Q)
}


# .procD.slide
# same as procD.slide, but without progress bar option
# used in pGpa.wSliders
.procD.slide <- function(curves, surf, Ya, ref, max.iter=10, tol){
    n <- length(Ya); p <- nrow(Ya[[1]]); k <- ncol(Ya[[1]])
    
    iter <- 1 # from initial rotation of Ya
    slid0 <- Ya
    Q <- ss0 <- sum(Reduce("+",Ya)^2)/n
    
    doTans <- !is.null(curves)
    
    if(doTans) {
      doSlide <- function(j)
        semilandmarks.slide.ProcD(slid0[[j]], tans[[j]],
                            surf, ref)
    } else {
      doSlide <- function(j)
        semilandmarks.slide.ProcD(slid0[[j]], tans = NULL, 
                            surf, ref)
    }
    
    while(Q > tol){
      iter <- iter+1
      
      if(doTans) {
        tans <- Map(function(y) tangents(curves, y, scaled=TRUE), 
                    slid0)
      } else tans <- NULL
      
      slid <- lapply(1:n, doSlide)
      
      ss <- sum(Reduce("+",slid)^2)/n
      slid0 <- apply.pPsup(ref,slid)
      slid <- NULL
      ref = cs.scale(Reduce("+", slid0)/n)
      Q <- abs(ss0-ss)
      ss0 <- ss
      if(iter >= max.iter) break
    }
    
    list(coords=slid0, consensus=ref, iter=iter+1, Q=Q)
}

# pGPA.wSliders
# GPA with partial Procrustes superimposition, incorporating semilandmarks
# used in gpagen
pGpa.wSliders <- function (Y, curves = NULL, surf = NULL, ProcD = TRUE, 
                           PrinAxes = FALSE, Proj = FALSE, appBE = TRUE, 
                           sen = 0.5, max.iter = 10, tol) {
  n <- length(Y)
  p <- nrow(Y[[1]])
  k <- ncol(Y[[1]])
  Yc <- Map(function(y) center.scale(y), Y)
  CS <- sapply(Yc, "[[", "CS")
  Ya <- lapply(Yc, "[[", "coords")
  Ya <- apply.pPsup(Ya[[1]], Ya)
  M <- Reduce("+", Ya)/n
  if (ProcD == FALSE) 
    gpa.slide <- BE.slide(curves = curves, surf = surf, Ya, 
                          ref = M, appBE = appBE, sen = 0.5, 
                          max.iter = max.iter, tol)
  else gpa.slide <- procD.slide(curves, surf, Ya, ref = M, 
                                max.iter = max.iter, tol)
  Ya <- gpa.slide$coords
  M <- gpa.slide$consensus
  iter <- gpa.slide$iter
  Q <- gpa.slide$Q
  if (PrinAxes == TRUE) {
    ref <- M
    rot <- prcomp(ref)$rotation
    for (i in 1:k) if (sign(rot[i, i]) != 1) 
      rot[1:k, i] = -rot[1:k, i]
    Ya <- Map(function(y) y %*% rot, Ya)
    M <- center.scale(Reduce("+", Ya)/n)$coords
  }
  list(coords = Ya, CS = CS, iter = iter, consensus = M, Q = Q, 
       nsliders = NULL)
}

# .pGPA.wSliders
# same as pGPA.wSliders, without option for progress bar
# used in gpagen
.pGpa.wSliders <- function(Y, curves = NULL, surf = NULL, ProcD = TRUE, 
                           PrinAxes = FALSE, Proj = FALSE, 
                           appBE = TRUE, sen = 0.5, max.iter = 10, tol, 
                           Parallel = FALSE){
  
  ParCores <- NULL
  if (is.numeric(Parallel)) {
    ParCores <- Parallel
    Parallel <- TRUE
  }
  if (Parallel && is.null(ParCores)) {
    ParCores <- detectCores() - 1
  }
  
  if (is.numeric(ParCores)) {
    if(ParCores > detectCores() - 1) ParCores <- detectCores() - 1
  } 
  
  n <- length(Y); p <- nrow(Y[[1]]); k <- ncol(Y[[1]])
  Yc <- Map(function(y) center.scale(y), Y)
  CS <- sapply(Yc,"[[","CS")
  Ya <- lapply(Yc,"[[","coords")
  Ya <- apply.pPsup(Ya[[1]], Ya)
  M <- Reduce("+", Ya)/n
  
  if(ProcD == FALSE) {
    gpa.slide <- if(Parallel) BE.slidePP(curves = curves, surf = surf, 
                                         Ya, ref=M, 
                                         appBE = appBE, sen = sen,
                                         max.iter = max.iter, 
                                         tol, ParCores = ParCores) else 
      .BE.slide(curves, surf, Ya, ref=M, appBE = appBE, sen = sen, 
                max.iter = max.iter, tol)
  } else {
    gpa.slide <- if(Parallel) procD.slidePP(curves = curves, surf = surf, Ya, ref=M, 
                               max.iter = max.iter, tol, ParCores = ParCores) else 
      .procD.slide(curves, surf, Ya, ref=M, max.iter = max.iter, tol)
    }
     
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

# perm.CR.index
# creates a permutation index for resampling, shuffling landmarks
# used in all functions utilizing CR (modularity)

perm.CR.index <- function(g, k, iter, seed=NULL){ # g is numeric partition.gp
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
  x <- center(as.matrix(x)); y <- center(as.matrix(y))
  px <- dim(x)[2]; py <- dim(y)[2]; pmin <- min(px,py)
  S12 <- crossprod(x, y)/(dim(x)[1] - 1)
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
# streamlines pls code
# used in all pls analyses
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

# pls.multi
# obtain average of pairwise PLS analyses for 3+modules
# used in all pls analyses
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
  if(length(CR.gp) > 1) {CR.mat <- dist(matrix(0, ngps,)) 
    for(i in 1:length(CR.mat)) CR.mat[[i]] <- CR.gp[i]
    CR.mat <- as.matrix(CR.mat) #added to specify which group is which (if out of numerical order)
    rownames(CR.mat) <- colnames(CR.mat) <- levels(factor(g, levels = unique(g)))
    CR.mat <- as.dist(CR.mat)
  }
  if(length(CR.gp)==1){  CR.mat <- NULL }
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

# phylo.mat
# estimate BM phylo.cov.matrix and transform from phylogeny
# used in: compare.evol.rates, compare.multi.evol.rates, phylo.integration, phylo.modularity, physignal
phylo.mat<-function(x,phy){
  C<-fast.phy.vcv(phy)
  C<-C[rownames(x),rownames(x)]
  invC <-fast.solve(C)
  eigC <- eigen(C)
  if(any(Im(eigC$values)==0)) eigC$values<-Re(eigC$values)
  if(any(Im(eigC$vectors)==0)) eigC$vectors<-Re(eigC$vectors)
  lambda <- zapsmall(abs(Re(eigC$values)))
  if(any(lambda == 0)){
    warning("Singular phylogenetic covariance matrix. Proceed with caution", immediate. = TRUE)
  }
  D.mat <- fast.solve(eigC$vectors%*% diag(sqrt(abs(eigC$values))) %*% t(eigC$vectors))
  rownames(D.mat) <- colnames(D.mat) <- colnames(C)
  list(invC = invC, D.mat = D.mat,C = C)
}

# pls.phylo
# phylogenetic pls
# used in: phylo.integration
pls.phylo <- function(x,y, Ptrans, verbose = FALSE){
  x <- as.matrix(x); y <- as.matrix(y)
  px <- ncol(x); py <- ncol(y); pmin <- min(px,py)
  x <- as.matrix(center(Ptrans %*% x))
  y <- as.matrix(center(Ptrans %*% y))
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
    warning("No specimens have scaling.  Unscaled landmarks imported, as a result\n", 
            immediate. = TRUE)
    scaled = FALSE
    } else {
      if(any(is.na(scaling))){
        sp.na <- which(is.na(scaling))
        warning(paste("\nWarning: Some specimens have no scaling:\n",
        sp.names[sp.na],
        "\nUnscaled landmarks imported, as a result.\n"), immediate. = TRUE)
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
# used in: readland.shapes and digit.curves
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

# get.VCV
# generates a VCV matrix, based on arguments
# used in K.modules

get.VCV <- function(A, phy = NULL, Cov = NULL,
                    transform. = TRUE) {
  x <- try(two.d.array(A), silent = TRUE)
  if(inherits(x, "try-error")) x <- try(as.matrix(A), silent = TRUE)
  if(inherits(x, "try-error"))
    stop("\nA is not a suitable data array for analysis. ", call. = FALSE)
  
  namesX <- rownames(x)
  if (is.null(namesX)) namesX <- 1:NROW(x)
  
  x <- as.matrix(x)
  n <- nrow(x)
  
  if(!is.null(phy) || !is.null(Cov)) {
    if(!is.null(phy) && is.null(Cov)) {
      Cov <- fast.phy.vcv(phy)
      Pcov <- Cov.proj(Cov, rownames(x))}
    
    if(is.null(phy) && !is.null(Cov)) {
      Pcov <- try(Cov.proj(Cov, rownames(x)), silent = TRUE)
      if(inherits(Pcov, "try-error"))
        stop("The names of Covariance matrix do not seem to match data names.\n",
             call. = FALSE)
    }
    
    if(!is.null(phy) && !is.null(Cov)) {
      Pcov <- try(Cov.proj(Cov, rownames(x)), silent = TRUE)
      if(inherits(Pcov, "try-error"))
        stop("The names of Covariance matrix do not seem to match data names.\n",
             call. = FALSE)
    }
    
    ones <- matrix(1, n)
    fit <- lm.fit(Pcov %*% ones, Pcov %*% x)
    B <- fit$coefficients
    R <- x - ones %*% B
    V <- if(transform.) crossprod(fit$resdiuals) / (n-1) else
      crossprod(R) / (n-1)
    
  } else V <- var(x)
  
  return(V)
}



### All functions below perform in a constrained way the same as functions in ape and geiger:
### -----------------------------------------------------------------------------------------


# A function to extract eigen values, as per mvrnorm
# But does not incorporate other traps and options, as in mvnorm.
Sig.eigen <- function(Sig, tol = 1e-06){
  E <- eigen(Sig)
  if(any(Im(E$values)==0)) E$values<-Re(E$values)
  if(any(Im(E$vectors)==0)) E$vectors<-Re(E$vectors)
  ev <- E$values
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
    if(any(Im(e$values)==0)) e$values<-Re(e$values)
    if(any(Im(e$vectors)==0)) e$vectors<-Re(e$vectors)
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
  if(any(Im(e$values)==0)) e$values<-Re(e$values)
  if(any(Im(e$vectors)==0)) e$vectors<-Re(e$vectors)
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

# physignal.z support functions

updateCov <- function(Cov, lambda) {
  D <- diag(diag(Cov))
  Cn <- (Cov - D) * lambda
  Cn + D
}

logLikh <- function(y, Cov = NULL, Pcov = NULL){
  y <- as.matrix(y)
  n <- nrow(y)
  p <- ncol(y)
  u <- matrix(1, n)
  gls <- !is.null(Pcov)
  yhat <- if(gls) u %*% lm.fit(Pcov %*% u, Pcov %*% y)$coefficients else fastFit(u, y, n, p)
  R <- as.matrix(y - yhat)
  Sig <- if(gls) crossprod(Pcov %*% R)/n else crossprod(R) / n
  logdetC <- if(gls) determinant(Cov, logarithm = TRUE)$modulus else 0
  logdetSig <- determinant(Sig, logarithm = TRUE)$modulus
  
  ll <- -0.5 * (n * p * log(2 * pi) + p * logdetC +
                  n * logdetSig + n) 
  
  as.numeric(ll)
}


lambda.opt<- function(y, tree) {
  y <- as.matrix(y)
  dims <- dim(y)
  n <- dims[1]
  p <- dims[2]
  Cov <- fast.phy.vcv(tree)
  
  LL <- function(lambda) {
    Cov.i <- updateCov(Cov, lambda)
    Cov.i <- Cov.i[rownames(y), rownames(y)]
    Pcov <- Cov.proj(Cov.i)
    logLikh(y, Cov = Cov.i, Pcov = Pcov)
  }
  opt <- optimise(LL, interval = c(1e-6, 1), maximum = TRUE, tol = 0.0001)$maximum
  
  if(opt < 1e-6) opt <- 1e-6
  if(opt > 0.999) opt <- 1
  
  opt
}

get.VCV <- function(A, phy = NULL, Cov = NULL,
                    transform. = TRUE) {
  x <- try(two.d.array(A), silent = TRUE)
  if(inherits(x, "try-error")) x <- try(as.matrix(A), silent = TRUE)
  if(inherits(x, "try-error"))
    stop("\nA is not a suitable data array for analysis. ", call. = FALSE)
  
  namesX <- rownames(x)
  if(is.null(namesX)) namesX <- 1:NROW(x)
  x <- as.matrix(x)
  n <- nrow(x)
  ones <- matrix(1, n)
  
  gls <- (!is.null(Cov) || !is.null(phy)) 
  
  if(gls) {
    if(!is.null(Cov)) {
      if(!is.null(phy)) {
        cat("\nBoth a phylogeny and covariance matrix were provided.")
        cat("\nOnly the covariance matrix will be used.\n")
      }
        
      Pcov <- try(Cov.proj(Cov, rownames(x)), silent = TRUE)
      if(inherits(Pcov, "try-error"))
        stop("The names of Covariance matrix do not seem to match data names.\n",
             call. = FALSE)
    } else {
      Cov <- fast.phy.vcv(phy)
      Pcov <- try(Cov.proj(Cov, rownames(x)), silent = TRUE)
      if(inherits(Pcov, "try-error")) {
        cat("\nUnable to match taxa names between tree and data.")
        cat("\nAssuming data have been sorted to match tree.\n")
        Pcov <- try(Cov.proj(Cov), silent = TRUE)
        if(inherits(Pcov, "try-error"))
          stop("Unable to create a covariance matrix from the tree.\n", 
               call. = FALSE)
      }
    }
    U <- qr.Q(qr(ones))
    if(transform.) {
      R <- fastLM(U, Pcov %*% x)$residuals
    } else {
      fit <- lm.fit(ones, as.matrix(Pcov%*% x))
      B <- coef(fit)
      R <- x - (ones %*% B)
    }
    
    V <-  crossprod(R) / (n-1)

    } else V <- var(x)
  
  return(V)
}
