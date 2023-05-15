#' Estimate locations of missing landmarks
#'
#' A function for estimating the locations of missing landmarks 
#' 
#' The function estimates the locations of missing landmarks for incomplete specimens  in a set of landmark
#' configurations, where missing landmarks in the incomplete specimens are designated by NA in place 
#' of the x,y,z coordinates.  Two distinct approaches are implemented.
#' 
#' The first approach (method="TPS") uses the thin-plate spline to interpolate landmarks on a reference specimen to estimate 
#' the locations of missing landmarks on a target specimen. Here, a reference specimen is obtained from 
#' the set of specimens for which all landmarks are present, Next, each incomplete specimen is aligned to 
#' the reference using the set of landmarks common to both. Finally, the thin-plate spline is used 
#' to estimate the locations of the missing landmarks in the target specimen (Gunz et al. 2009).
#' 
#' The second approach (method="Reg") is multivariate regression. Here each landmark with missing values is
#' regressed on all other landmarks for the set of complete specimens, and the missing landmark values are
#' then predicted by this linear regression model. Because the number of variables can exceed the number of
#' specimens, the regression is implemented on scores along the first set of PLS axes for the complete and 
#' incomplete blocks of landmarks (see Gunz et al. 2009). Note, however, that a minimum of k*m+k specimens 
#' are required to estimate m missing landmarks (of k-dimension) in any one specimen using the regression method.
#' More generally, if the number of missing landmarks approaches the number of reference specimens used to 
#' estimate them, estimation will become increasingly imprecise with the regression method. Additionally, 
#' the location of missing landmarks (contiguous versus disparate in location) can also influence the precision 
#' of estimation.  The user should be aware that the function will produce results but the results from the 
#' regression method might be influenced by the number of specimens, the number of total landmarks, and the 
#' number and location of missing landmarks in any one specimen.  It might be wise to compare multiple methods 
#' for specific cases, if uncertain about the precision of estimation.
#'  
#'  One can also exploit bilateral symmetry to estimate the locations of missing landmarks. Several
#'   possibilities exist for implementing this approach (see Gunz et al. 2009).  Example R code for one 
#'   implementation is found in Claude (2008).
#'  
#' NOTE: Because all geometric morphometric analyses and plotting functions implemented in geomorph 
#' require a full complement of landmark coordinates, the alternative to estimating the missing 
#' landmark coordinates is to proceed with subsequent analyses EXCLUDING
#' specimens with missing values. 
#' 
#' @param A An array (p x k x n) containing landmark coordinates for a set of specimens or a geomorphShapes object
#' @param method Method for estimating missing landmark locations
#' @author Dean Adams
#' @keywords utilities
#' @return Function returns an array (p x k x n) of the same dimensions as input A, including coordinates for the target specimens 
#' (the original landmarks plus the estimated coordinates for the missing landmarks).
#' If the input is a geomorphShapes object, this is returned with the original and estimated coordinates in $landmarks
#' In both cases, these data need to be Procrustes Superimposed prior to analysis, including sliding of semilandmarks (see \code{\link{gpagen}}).
#' @export
#' @references Claude, J. 2008. Morphometrics with R. Springer, New York.
#' @references  Bookstein, F. L., K. Schafer, H. Prossinger, H. Seidler, M. Fieder, G. Stringer, G. W. Weber, 
#' J.-L. Arsuaga, D. E. Slice, F. J. Rohlf, W. Recheis, A. J. Mariam, and L. F. Marcus. 1999. Comparing 
#' frontal cranial profiles in archaic and modern Homo by morphometric analysis. Anat. Rec. (New Anat.) 257:217-224.
#' @references Gunz, P., P. Mitteroecker, S. Neubauer, G. W. Weber, and F. L. Bookstein. 2009. Principles for 
#' the virtual reconstruction of hominin crania. J. Hum. Evol. 57:48-62.
#' @examples
#' data(plethodon)
#' plethland<-plethodon$land
#'   plethland[3,,2]<-plethland[8,,2]<-NA  #create missing landmarks
#'   plethland[3,,5]<-plethland[8,,5]<-plethland[9,,5]<-NA  
#'   plethland[3,,10]<-NA  
#'   
#' estimate.missing(plethland,method="TPS")
#' estimate.missing(plethland,method="Reg")
estimate.missing <- function(A, method=c("TPS","Reg")){
  if(inherits(A, "geomorphShapes")) {
    a <- A
    A <- simplify2array(a$landmarks)
  }
  
  dims <- dim(A)
  n <- dims[[3]]
  p <- dims[[1]]
  k <- dims[[2]]
  
  A.names <- dimnames(A)

  if(any(is.na(A))==FALSE)  {stop("No missing data.")}
  method <- match.arg(method)
  if (length(dim(A))!=3){
    stop("Data matrix not a 3D array (see 'arrayspecs').")  }
  if(method=="TPS"){
    A2 <- A
    spec.NA <- which(rowSums(is.na(two.d.array(A)))>0)
    Y.gpa <- gpagen(A[,,-spec.NA], PrinAxes=FALSE, print.progress = FALSE)
    p <- dim(Y.gpa$coords)[1]; k <- dim(Y.gpa$coords)[2]; n <- dim(Y.gpa$coords)[3]    
    ref <- mshape(arrayspecs(two.d.array(Y.gpa$coords)*Y.gpa$Csize, p, k))
    for (i in 1:length(spec.NA)){
      missing <- which(is.na(A2[,1,spec.NA[i]]))
      tmp <- tps2d3d(ref, ref[-missing,], A2[-missing,,spec.NA[i]], PB=FALSE)
      A2[,,spec.NA[i]] <- tmp
    }
  }
  if(method=="Reg"){
    spec.NA <- which(rowSums(is.na(two.d.array(A)))>0)    
    land.NA <- which(colSums(is.na(two.d.array(A)))>0)  
    p <- dim(A)[1]; k <- dim(A)[2]  
    A2 <- A
    complete <- A[,,-spec.NA]
    incomplete <- A[,,spec.NA]
    Y.gpa <- gpagen(complete, PrinAxes=FALSE, print.progress = FALSE)
    ref <- mshape(arrayspecs(two.d.array(Y.gpa$coords)*Y.gpa$Csize, p, k))
    complete <- arrayspecs(two.d.array(Y.gpa$coords)*Y.gpa$Csize, p, k)
    if(length(dim(incomplete))>2){
      for (i in 1:dim(incomplete)[3]){
        missing <- which(is.na(incomplete[,1,i]))
        lndmk <- which(!is.na(incomplete[,1,i]))
        tmp <- apply.pPsup(center(ref[-missing,]), list(center(incomplete[-missing,,i])))[[1]]
        incomplete[lndmk,,i] <- tmp
      }      
    }
    if(length(dim(incomplete))==2){
      missing <- which(is.na(incomplete[,1]))
      lndmk <- which(!is.na(incomplete[,1]))
      tmp <- apply.pPsup(center(ref[-missing,]), list(center(incomplete[-missing,])))[[1]]
      incomplete[lndmk,] <- tmp
    }
    A2[,,-spec.NA] <- complete
    A2[,,spec.NA] <- incomplete
    
    A2.list <- lapply(1:n, function(j) A2[,,j])
    A2.temp <- lapply(A2.list, na.omit)
    keep <- sapply(A2.temp, function(x) NROW(x) == p)
    A2.complete <- A2.list[keep]
    A2.2d.complete <- t(sapply(A2.complete, as.vector))
    V <- crossprod(A2.2d.complete)
    
    A2.updated <- lapply(A2.list, function(x){
      missing.coord <- attr(na.omit(as.vector(x)), "na.action")
      if(!is.null(missing.coord)) {
        Vr <- V[-missing.coord, missing.coord]
        pls <- svd(Vr)
        UU <- pls$u
        VV <- pls$v
        xx <- A2.2d.complete[, -missing.coord]
        yy <- A2.2d.complete[, missing.coord]
        XScores <- xx %*% UU
        YScores <- yy %*% VV
        beta <- coef(lm(YScores ~ XScores ))
        miss.xsc <- c(1,  as.vector(x)[-missing.coord] %*% UU)
        miss.ysc <- miss.xsc %*% beta
        pred.val <- miss.ysc %*% t(VV)
        x[missing.coord] <- pred.val
      } else x <- x
      x
      
    })
    
    names(A2.updated) <- A.names[[3]]
    NA.check <- sapply(A2.updated, function(x) any(is.na(x)))
    if(any(NA.check)) {
      
      oldwarn <- options()$warn
      options(warn = 1)
      warning("Some specimens could not have all landmarks estimated.\n")
      
      cat("This is probably because of too many missing landmarks.\n")
      cat("Consider using method = 'TPS' instead.\n")
      cat("Specimens with NA:\n")
      print(names(A2.updated[which(NA.check)]))
      options(warn = oldwarn)
    }

    A2 <- simplify2array(A2.updated)
    dimnames(A2) <- A.names
  }
  if("a" %in% ls()) {
    a$landmarks <- A2.updated
    names(A2.updated) <- A.names[[3]]
    A2 <- a
  }
   return(A2)
}
