## S3 GENERIC FUNCTIONS

## gpagen

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.gpagen <- function (x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call), fill=TRUE, "\n\n")
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

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.gpagen <- function(object, ...) {
  x <- object
  print.gpagen(x, ...)
}

#' Plot Function for geomorph
#' 
#' @param x plot object
#' @param ... other arguments passed to plot
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.gpagen <- function(x, ...){
  plotAllSpecimens(x$coords)
}


## procD.lm

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.procD.lm <- function (x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call), fill=TRUE, "\n\n")
  cat("\nType I (Sequential) Sums of Squares and Cross-products\n")
  if(x$perm.method == "RRPP") cat ("Randomized Residual Permutation Procedure Used\n") else
    cat("Randomization of Raw Values used\n")
  cat(paste(x$permutations, "Permutations"))
  if(!is.null(x$random.F)) cat("\n\n*** F values, Z scores, and P values updated for nested effects")
  cat("\n\n")
  print(x$aov.table)
  invisible(x)
}

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.procD.lm <- function(object, ...) {
  x <- object
  print.procD.lm(x, ...)
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

#' Plot Function for geomorph
#' 
#' @param x plot object
#' @param outliers Logical argument to include outliers plot
#' @param ... other arguments passed to plot
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.procD.lm <- function(x, outliers=FALSE, ...){
  r <- x$residuals
  f <- x$fitted
  if(!is.null(x$weights)) {r <- r*sqrt(x$weights); f <- f*sqrt(x$weights)}
  if(!is.null(x$Pcor)) {
    y <- x$Pcor%*%x$Y
    X <- x$Pcor%*%x$X
    fit <- lm.fit(X,y)
    f <- fit$fitted.values
    r <- fit$residuals
  }
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


## advanced.procD.lm

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.advanced.procD.lm <- function (x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call), fill=TRUE, "\n\n")
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

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.advanced.procD.lm <- function(object, ...) {
  x <- object
  print.advanced.procD.lm(x, ...)
}

#' Plot Function for geomorph
#' 
#' @param x plot object
#' @param outliers Logical argument to include outliers plot
#' @param ... other arguments passed to plot
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.advanced.procD.lm <- function(x, outliers = FALSE, ...) {
  plot.procD.lm(x,...)
}


## procD.allometry

# Two print options

printAllometry.HOS <- function(x){
  cat("\nCall:\n")
  cat(deparse(x$call), fill=TRUE, "\n\n")
  cat("\nHomogeneity of Slopes Test\n")
  print(x$HOS.test)
  if(x$HOS.test[2,7] > x$alpha) cat(paste("\nThe null hypothesis of parallel slopes is supported
  based on a signficance criteron of alpha =", x$alpha,"\n")) 
  if(x$HOS.test[2,7] <= x$alpha) cat(paste("\nThe null hypothesis of parallel slopes is rejected
  based on a signficance criteron of alpha =", x$alpha,"\n"))
  cat("\nBased on the results of this test, the following ANOVA table is most appropriate\n")
  cat("\nType I (Sequential) Sums of Squares and Cross-products\n")
  if(x$perm.method == "RRPP") cat ("Randomized Residual Permutation Procedure Used\n") else
    cat("Randomization of Raw Values used\n")
  cat(paste(x$permutations, "Permutations"))
  cat("\n\n")
  print(x$aov.table)
}

printAllometry.noHOS <- function(x){
  cat("\nCall:\n")
  cat(deparse(x$call), fill=TRUE, "\n\n")
  cat("\nType I (Sequential) Sums of Squares and Cross-products\n")
  if(x$perm.method == "RRPP") cat ("Randomized Residual Permutation Procedure Used\n") else
    cat("Randomization of Raw Values used\n")
  cat(paste(x$permutations, "Permutations"))
  cat("\n\n")
  print(x$aov.table)
}

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.procD.allometry <- function (x, ...) {
  if(!is.null(x$HOS.test)) printAllometry.HOS(x) else printAllometry.noHOS(x)
  invisible(x)
}

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.procD.allometry <- function(object, ...) {
  x <- object
  print.procD.allometry(x,...)
}

#' Plot Function for geomorph
#' 
#' @param x plot object
#' @param method Method for estimating allometric shape components; see \code{\link{procD.allometry}} for details
#' @param warpgrids A logical value indicating whether deformation grids for small and large shapes 
#'  should be displayed (note: if groups are provided no TPS grids are shown)
#' @param label An optional vector indicating labels for each specimen that are to be displayed
#' @param mesh A mesh3d object to be warped to represent shape deformation of the minimum and maximum size 
#' if {warpgrids=TRUE} (see \code{\link{warpRefMesh}}).
#' @param ... other arguments passed to plot
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.procD.allometry <- function(x, method=c("CAC","RegScore","PredLine"),warpgrids=TRUE,
                                 label=NULL, mesh=NULL, ...) {
  method <- match.arg(method)
  if(x$logsz) xlab <- "log(Size)" else xlab <- "Size"
  if(x$logsz) size <- log(x$size) else size <- x$size
  n <- length(size)
  if(!is.null(x$gps)) pt.cols <- as.numeric(x$gps) else pt.cols <- rep(1, length(size))
  if(method == "CAC"){
    layout(matrix(c(3,1,1,1,1,1,1,1,4,2,2,2,2,2,2,2,2,2),3,6))   
    plot(size,x$CAC,xlab=xlab, ylab="CAC",pch=21,bg=pt.cols,cex=1.25)
    if (length(label!=0)) {
      if(isTRUE(label)){text(size,x$CAC,seq(1, n),adj=c(-0.7,-0.7)) }
      else{text(size,x$CAC,label,adj=c(-0.1,-0.1))}
    }
    plot(x$CAC,x$RSC[,1], xlab="CAC",ylab="RSC 1", pch=21,bg=pt.cols,cex=1.25)
    if (!is.null(label)) {
      if(!is.null(label)){text(x$CAC,x$RSC,seq(1, n),adj=c(-0.7,-0.7)) }
      else{text(x$CAC,x$RSC,label,adj=c(-0.1,-0.1))}
    }
  }
  if(method=="PredLine"){
    layout(matrix(c(2,1,1,1,1,1,1,1,3),3,3))   
    plot(size,x$pred.val,xlab=xlab, ylab="Shape (Predicted)",pch=21,bg=pt.cols,cex=1.25)
    if (length(label!=0)) {
      if(isTRUE(label)){text(size,x$pred.val,seq(1, n),adj=c(-0.7,-0.7)) }
      else{text(size,x$pred.val,label,adj=c(-0.1,-0.1))}
    }
  }
  if(method=="RegScore"){
    layout(matrix(c(2,1,1,1,1,1,1,1,3),3,3))   
    plot(size,x$Reg.proj,xlab=xlab, ylab="Shape (Regression Score)",pch=21,bg=pt.cols,cex=1.25)
    if (length(label!=0)) {
      if(isTRUE(label)){text(size,x$Reg.proj,seq(1, n),adj=c(-0.7,-0.7)) }
      else{text(size,x$Reg.proj,label,adj=c(-0.1,-0.1))}
    }
  } 
  if(method=="CAC") y <- x$CAC else if (method=="RegScore") y <- x$Reg.proj else
    y <- x$pred.val
  
  if(is.null(x$gps)){
    if(warpgrids==TRUE){
      arrows(min(size), (0.7 * max(y)), min(size), 0, length = 0.1,lwd = 2)
      arrows(max(size), (0.7 * min(y)), max(size), 0, length = 0.1,lwd = 2)
    }
  }
  if(warpgrids==TRUE && x$k==3){
    if(is.null(mesh)){
      open3d()
      plot3d(x$Ahat[,,which.min(size)],type="s",col="gray",main="Shape at minimum size",size=1.25,aspect=FALSE)
      open3d()
      plot3d(x$Ahat[,,which.max(size)],type="s",col="gray",main="Shape at maximum size",size=1.25,aspect=FALSE)
      if(!is.null(mesh)){
        plotRefToTarget(x$ref, x$Ahat[,,which.min(size)], mesh, method = "surface")
        title3d(main="Shape at minimum size")
        plotRefToTarget(x$ref, x$Ahat[,,which.max(size)], mesh, method = "surface")
        title3d(main="Shape at maximum size")
      }}
  }
  if(warpgrids==TRUE && x$k==2){
    tps(x$ref,x$Ahat[,,which.min(size)],20)
    tps(x$ref,x$Ahat[,,which.max(size)],20)
  }
  layout(1) 
}



## morphol.disparity

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.morphol.disparity <- function (x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call), fill=TRUE, "\n\n")
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

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.morphol.disparity <- function(object, ...) {
  x <- object
  print.morphol.disparity(x, ...)
}


## pls

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.pls <- function (x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call), fill=TRUE, "\n\n")
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

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.pls <- function(object, ...) {
  x<- object
  print.pls(x, ...)
}

plotPLS <- function(p, label = NULL, warpgrids=TRUE){
  A1 <- p$A1; A2 <- p$A2
  XScores <- p$XScores; YScores <- p$YScores
  if(is.matrix(XScores)) XScores <- XScores[,1]
  if(is.matrix(YScores)) YScores <- YScores[,1]
  pc <- prcomp(cbind(XScores, YScores))$x[,1]
  px <- predict(lm(XScores~pc))
  py <- predict(lm(YScores~pc))
  pxmax <- max(px); pxmin <- min(px)
  pymax <- max(py); pymin <- min(py)
  
  if (length(dim(A1)) == 3) {
    A1.ref <- mshape(A1)
    A1.min <- arrayspecs(pxmin*p$left.pls.vectors[,1], 
                         nrow(A1.ref), ncol(A1.ref))[,,1]
    A1.max <- arrayspecs(pxmax*p$left.pls.vectors[,1], 
                         nrow(A1.ref), ncol(A1.ref))[,,1]
    pls1.min <- A1.ref + A1.min
    pls1.max <- A1.ref + A1.max
  }
  if (length(dim(A2)) == 3) {
    A2.ref <- mshape(A2)
    A2.min <- arrayspecs(pymin*p$right.pls.vectors[,1], 
                         nrow(A2.ref), ncol(A2.ref))[,,1]
    A2.max <- arrayspecs(pymax*p$right.pls.vectors[,1], 
                         nrow(A2.ref), ncol(A2.ref))[,,1]
    pls2.min <- A2.ref + A2.min
    pls2.max <- A2.ref + A2.max
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
                          0.8, 1), byrow = TRUE, ncol = 4))
    screen(1)
    plot(XScores, YScores, pch = 21, bg = "black", 
         main = "PLS1 Plot: Block 1 (X) vs. Block 2 (Y) ", 
         xlab = "PLS1 Block 1", ylab = "PLS1 Block 2")
    abline(lm(py~px), col="red")
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

#' Plot Function for geomorph
#' 
#' @param x plot object
#' @param label Optional vector to label points
#' @param warpgrids Logical argument whether to include warpgrids
#' @param ... other arguments passed to plot
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.pls <- function(x, label = NULL, warpgrids=TRUE, ...){
  plotPLS(x, label=label, warpgrids=warpgrids)
}


## bilat.symmetry

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.bilat.symmetry <- function (x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call), fill=TRUE, "\n\n")
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

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.bilat.symmetry <- function(object, ...) {
  x <- object
  print.bilat.symmetry(x, ...)
}

plotBilatSymmetry <- function(b, warpgrids = TRUE, mesh= NULL){
  k <- dim(b$symm.shape)[[2]]
  if(b$data.type == "Matching"){
    if(k==2){  
      par(mfrow=c(2,2),oma=c(1.5,0,1.5,0))
      plotAllSpecimens(b$symm.shape)
      plotAllSpecimens(b$asymm.shape)
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
        plotAllSpecimens(b$symm.shape)
        plotAllSpecimens(b$asymm.shape)
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

#' Plot Function for geomorph
#' 
#' @param x plot object
#' @param warpgrids Logical argument whether to include warpgrids
#' @param mesh Option to include mesh in warpgrids plots
#' @param ... other arguments passed to plot
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.bilat.symmetry <- function(x, warpgrids = TRUE, mesh= NULL, ...){
  plotBilatSymmetry(x, warpgrids = warpgrids, mesh = mesh)
}


## CR

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.CR <- function (x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call), fill=TRUE, "\n\n")
  cat(paste("\nCR:", round(x$CR, nchar(x$permutations))))
  cat(paste("\n\nP-value:", round(x$P.value, nchar(x$permutations))))
  cat(paste("\n\nBased on", x$permutations, "random permutations"))
  if(!is.null(x$CInterval)) cat(paste("\n\nConfidence Intervals", round(x$CInterval,nchar(x$permutations)))) 
  invisible(x)
}


#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.CR <- function(object, ...) {
  x <- object
  print.CR(x, ...)
}

#' Plot Function for geomorph
#' 
#' @param x plot object
#' @param ... other arguments passed to plot
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.CR <- function(x, ...){
  CR.val <- x$random.CR
  CR.obs <- x$CR
  p <- x$P.value
  ndec <- nchar(p)-2
  CR.obs <- round(CR.obs, ndec)
  main.txt <- paste("Observed CR =",CR.obs,";", "P-value =", p)
  hist(CR.val,30,freq=TRUE,col="gray",xlab="CR Coefficient",xlim=c(0,max(c(2,CR.val))),
       main=main.txt, cex.main=0.8)
  arrows(CR.obs,50,CR.obs,5,length=0.1,lwd=2)
}


#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Dean Adams
#' @keywords utilities
print.CR.phylo <- function (x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call), fill=TRUE, "\n\n")
  cat(paste("\nCR:", round(x$CR, nchar(x$permutations))))
  cat(paste("\n\nP-value:", round(x$P.value, nchar(x$permutations))))
  cat(paste("\n\nBased on", x$permutations, "random permutations"))
  if(!is.null(x$CInterval)) cat(paste("\n\nConfidence Intervals", round(x$CInterval,nchar(x$permutations)))) 
  invisible(x)
}

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Dean Adams
#' @keywords utilities
summary.CR.phylo <- function(object, ...) {
  x <- object
  print.CR.phylo(x, ...)
}

#' Plot Function for geomorph
#' 
#' @param x plot object
#' @param ... other arguments passed to plot
#' @export
#' @author Dean Adams
#' @keywords utilities
#' @keywords visualization
plot.CR.phylo <- function(x, ...){
  CR.val <- x$random.CR
  CR.obs <- x$CR
  p <- x$P.value
  ndec <- nchar(p)-2
  CR.obs <- round(CR.obs, ndec)
  main.txt <- paste("Observed CR =",CR.obs,";", "P-value =", p)
  hist(CR.val,30,freq=TRUE,col="gray",xlab="CR Coefficient",xlim=c(0,max(c(2,CR.val))),
       main=main.txt, cex.main=0.8)
  arrows(CR.obs,50,CR.obs,5,length=0.1,lwd=2)
}

## physignal

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.physignal <- function(x, ...){
  cat("\nCall:\n")
  cat(deparse(x$call), fill=TRUE, "\n\n")
  cat(paste("\nObserved Phylogenetic Signal (K):", round(x$phy.signal, nchar(x$permutations))))
  cat(paste("\n\nP-value:", round(x$pvalue, nchar(x$permutations))))
  cat(paste("\n\nBased on", x$permutations, "random permutations"))
  invisible(x)
}

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.physignal <- function(object, ...) {
  x <- object
  print.physignal(x, ...)
}

#' Plot Function for geomorph
#' 
#' @param x plot object
#' @param ... other arguments passed to plot
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.physignal <- function(x, ...){
  K.val <- x$random.K
  K.obs <- x$phy.signal
  p <- x$pvalue
  ndec <- nchar(p)-2
  K.obs <- round(K.obs, ndec)
  main.txt <- paste("Observed K =",K.obs,";", "P-value =", p)
  hist(K.val,30,freq=TRUE,col="gray",xlab="Phylogenetic Signal, K",
       main=main.txt, cex.main=0.8)
  arrows(K.obs,50,K.obs,5,length=0.1,lwd=2)
}


## evolrate

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.evolrate <- function (x, ...) {
  cat("\nCall:\n")
  cat(paste("\n\nObserved Rate Ratio:", round(x$sigma.d.ratio, nchar(x$permutations))))
  cat(paste("\n\nP-value:", round(x$P.value, nchar(x$permutations))))
  cat(paste("\n\nBased on", x$permutations, "random permutations"))
  cat(paste("\n\nThe rate for group",x$groups,"is",x$sigma.d.gp, ""))
  invisible(x)
}

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.evolrate <- function(object, ...) {
  x <- object
  print.evolrate(x, ...)
}

## evolrate1

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.evolrate1 <- function (x, ...) {
  cat("\nCall:\n")
  cat(paste("\n\nOne group only. Observed Rate:", x$sigma.d.all))
  invisible(x)
}

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.evolrate1 <- function(object, ...) {
  x <- object
  print.evolrate1(x, ...)
}

#' Plot Function for geomorph
#' 
#' @param x plot object
#' @param ... other arguments passed to plot
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.evolrate <- function(x, ...){
  Rate.val <- x$random.sigma
  Rate.obs <- x$random.sigma[1]
  p <- x$P.value
  ndec <- nchar(x$permutations)
  Rate.obs <- round(Rate.obs, ndec)
  p <- round(p, ndec)
  main.txt <- paste("Observed Rate Ratio =",Rate.obs,";", "P-value =", p)
  hist(Rate.val,30,freq=TRUE,col="gray",xlab="Rate Ratios",xlim=c(0,max(c(2,Rate.val))),
       main=main.txt, cex.main=0.8)
  arrows(Rate.obs,50,Rate.obs,5,length=0.1,lwd=2)
}


## trajectory.analysis

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param angle.type Choice between vector correlation or vector angles, in radians or degrees
#' for summarizing results ("r", "rad", "deg", respectively)
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.trajectory.analysis <- function(x, 
                angle.type = c("r", "rad", "deg"), ...) {
  angle.type = match.arg(angle.type)
  cat(deparse(x$call), fill=TRUE, "\n\n")
  cat("\nType I (Sequential) Sums of Squares and Cross-products\n")
  cat ("Randomized Residual Permutation Procedure Used\n") 
  cat(paste(x$permutations, "Permutations"))
  cat("\n\n")
  print(x$aov.table)
  cat("\n\n")
  cat("Pairwise statistical results:\n")
  cat("\n\n*** Path Distances\n")
  cat("\nObserved Path Distances\n")
  print(x$path.distances)
  cat("\nPairwise Absolute Differences Between Path Distances\n")
  print(x$magnitude.diff)
  cat("\nEffect Sizes\n")
  print(x$Z.magnitude.diff)
  cat("\nP-Values\n")
  print(x$P.magnitude.diff)
  cat("\n")
  if(angle.type == "r"){
    cat("\n*** Principal Vector Correlations\n")
    cat("\nPairwise Correlations\n")
    print(x$trajectory.cor)
    cat("\nEffect Sizes\n")
    print(x$Z.angle)
    cat("\nP-Values\n")
    print(x$P.angle)
    cat("\n")
  } else if(angle.type == "rad"){
    cat("\n*** Principal Vector Angles\n")
    cat("\nPairwise Angles (in radians)\n")
    print(x$trajectory.angle.rad)
    cat("\nEffect Sizes\n")
    print(x$Z.angle)
    cat("\nP-Values\n")
    print(x$P.angle)
    cat("\n")
  } else {
    cat("\n*** Principal Vector Angles\n")
    cat("\nPairwise Angles (in degrees)\n")
    print(x$trajectory.angle.deg)
    cat("\nEffect Sizes\n")
    print(x$Z.angle)
    cat("\nP-Values\n")
    print(x$P.angle)
    cat("\n")
  }
  tr1 <- x$random.trajectories[[1]][[1]]
  tp <- nrow(tr1)
  if(tp > 2){
    cat("\n*** Trajectory Shape Differences\n")
    cat("\nPairwise Shape Differences (Procrustes Distance)\n")
    print(x$trajectory.shape.dist)
    cat("\nEffect Sizes\n")
    print(x$Z.shape.diff)
    cat("\nP-Values\n")
    print(x$P.shape.diff)
  }
  invisible(x)
}

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param angle.type Choice between vector correlation or vector angles, in radians or degrees
#' for summarizing results ("r", "rad", "deg", respectively)
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.trajectory.analysis <- function(object,
                  angle.type = c("r", "rad", "deg"), ...) {
  x <- object
  angle.type <- match.arg(angle.type)
  print.trajectory.analysis(x, angle.type=angle.type, ...)
}

# general plotting functions for phenotypic trajectories
trajplot.w.int<-function(Data, M, TM, groups, group.cols = NULL, 
pattern = c("white", "gray", "black"), pt.scale = 1, ...){ # TM = trajectories from means
  n <- length(TM); tp<-dim(TM[[1]])[1]; p<-dim(TM[[1]])[2]
  if(length(pattern) != 3) stop("Point sequence color pattern must conatin three values")
  pmax <- max(Data[,1]); pmin <- min(Data[,1])
  plot(Data[,1:2],type="n",
       xlim = c(2*pmin, pmax),
       xlab="PC I", ylab="PC II",
       main="Two Dimensional View  of Phenotypic Trajectories",asp=1)
  
  if(is.null(group.cols)) gp.cols <- 1:n else gp.cols <- group.cols
  if(length(gp.cols) != nlevels(groups)) 
      stop("group.cols is not logical with respect to group levels") 
  
  points(Data[,1:2],pch=21,bg=pattern[2],cex=.75*pt.scale)
  # Sequence lines
  for(i in 1:n){
    y <- TM[[i]]
    for(ii in 1:(tp-1)) points(y[ii:(ii+1),1:2],  type="l", lwd=1.5*pt.scale, col=gp.cols[i])
  }
  # Sequence points
  points(M[,1:2], pch=21, bg=pattern[2], cex=1.5*pt.scale)
  for(i in 1:n){
    y <- TM[[i]]
    k <- nrow(y)
      points(y[1,1], y[1,2], pch=21, cex=1.5*pt.scale, bg=pattern[1])
      points(y[k,1], y[k,2], pch=21, cex=1.5*pt.scale, bg=pattern[3])
  }

  legend("topleft", levels(groups), lwd=2, col=gp.cols)
}

trajplot.by.groups<-function(Data, TM, groups, group.cols = NULL, 
            pattern = c("white", "gray", "black"), pt.scale = 1, ...) {
  n <- length(TM); tp <- nrow(TM[[1]]); p <- ncol(TM[[1]])
  if(length(pattern) != 3) stop("Point sequence color pattern must conatin three values")
  Data2 <- t(matrix(matrix(t(Data)),p,))
  pmax <- max(Data[,1]); pmin <- min(Data2[,1])
  plot(Data2[,1:2], type="n",
       xlim = c(2*pmin, pmax),
       xlab="PC I", ylab="PC II",
       main="Two Dimensional View  of Phenotypic Trajectories",asp=1)
  if(is.null(group.cols)) gp.cols <- as.numeric(groups) else gp.cols <- group.cols
  if(length(gp.cols) != length(groups)) {
    if(length(gp.cols) != nlevels(groups)) 
      stop("group.cols is not logical with respect to either groups or group levels") else
      {
         new.gp.cols <-array(,n)
         for(i in 1:n) new.gp.cols[i] <- gp.cols[match(groups[i], levels(groups))]
      } 
    gp.cols <- new.gp.cols
  }
  gp.index <- unique(gp.cols)
  point.seq <- function(x, p, tp, pt.col, pt.scale){
    for(i in 1:(tp-1)){
      y <- matrix(x[1:(2*p)],2,, byrow=TRUE)
      points(y, type="l", col=pt.col, lwd=1*pt.scale)
      x <- x[-(1:p)]
    }
  }
  for(i in 1:nrow(Data)) point.seq(Data[i,], p=p, tp=tp, pt.col = gp.cols[i], pt.scale=pt.scale)
  points(Data2, pch = 21, bg = pattern[2], cex = 1*pt.scale)
  for(i in 1:length(TM)){
    y <- TM[[i]]
    points(y[1,1], y[1,2], pch=21, bg = pattern[1], cex=1*pt.scale)
    points(y[nrow(y),1], y[nrow(y),2], pch=21, bg = pattern[3], cex=1*pt.scale)
  }
  legend("topleft", levels(groups), lwd=2, col=gp.index)
}

#' Plot Function for geomorph
#' 
#' @param x plot object
#' @param group.cols An optional vector of colors for group levels
#' @param pt.seq.pattern The sequence of colors for starting, middle, and end points of 
#' trajectories, respectivly.  E.g., c("green", "gray", "red") for gray points
#' but initial points with green color and end points with red color.
#' @param pt.scale An optional value to magnify or reduce points (1 = no change)
#' @param ... other arguments passed to plot
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.trajectory.analysis <- function(x, group.cols = NULL, 
            pt.seq.pattern  = c("white", "gray", "black"), pt.scale = 1,...){
  if(x$trajectory.type == 2)
  trajplot.w.int(Data=x$pc.data, M =x$pc.means,
           TM = x$pc.trajectories, groups = x$groups, 
           group.cols=group.cols, pattern = pt.seq.pattern, pt.scale=pt.scale)
  if(x$trajectory.type == 1)
    trajplot.by.groups(Data=x$pc.data, 
           TM = x$pc.trajectories, groups = x$groups, 
           group.cols=group.cols, pattern = pt.seq.pattern, pt.scale=pt.scale)
}


