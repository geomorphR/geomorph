## S3 GENERIC FUNCTIONS

## gpagen

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
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

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
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

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
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

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
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
plot.advanced.procD.lm <- function(x, outliers = FALSE, ...) {
  plot.procD.lm(x,...)
}

## morphol.disparity

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
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

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
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

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
summary.pls <- function(object, ...) {
  x<- object
  print.pls(x, ...)
}

plotPLS <- function(p, label = NULL, warpgrids=TRUE){
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

#' Plot Function for geomorph
#' 
#' @param x plot object
#' @param label Optional vector to label points
#' @param warpgrids Logical argument whether to include warpgrids
#' @param ... other arguments passed to plot
#' @export
#' @author Michael Collyer
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

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
summary.bilat.symmetry <- function(object, ...) {
  x <- object
  print.bilat.symmetry(x, ...)
}

plotBilatSymmetry <- function(b, warpgrids = TRUE, mesh= NULL){
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
#' Plot Function for geomorph
#' 
#' @param x plot object
#' @param warpgrids Logical argument whether to include warpgrids
#' @param mesh Option to include mesh in warpgrids plots
#' @param ... other arguments passed to plot
#' @export
#' @author Michael Collyer
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
print.CR <- function (x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call), "\n\n")
  cat(paste("\nCR:", round(x$CR, nchar(x$permutations))))
  cat(paste("\n\nP-value:", round(x$P.value, nchar(x$permutations))))
  cat(paste("\n\nBased on", x$permutations, "random permutations"))
  cat(paste("\n\nConfidence Intervals", round(x$CInterval,nchar(x$permutations))))
  invisible(x)
}


#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
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
print.CR.phylo <- function (x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call), "\n\n")
  cat(paste("\nCR:", round(x$CR, nchar(x$permutations))))
  cat(paste("\n\nP-value:", round(x$P.value, nchar(x$permutations))))
  cat(paste("\n\nBased on", x$permutations, "random permutations"))
  invisible(x)
}

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Dean Adams
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
print.physignal <- function(x, ...){
  cat("\nCall:\n")
  cat(deparse(x$call), "\n\n")
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
print.evolrate <- function (x, ...) {
  cat("\nCall:\n")
  cat(paste("\n\nObserved Rate Ratio:", round(x$sigma.d.ratio, nchar(x$permutations))))
  cat(paste("\n\nP-value:", round(x$P.value, nchar(x$permutations))))
  cat(paste("\n\nBased on", x$permutations, "random permutations"))
  cat(paste("\n\nThe Rates for each group:", x$sigma.d.gp, ""))
  invisible(x)
}

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
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





