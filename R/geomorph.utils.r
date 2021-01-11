## S3 GENERIC FUNCTIONS

## gpagen

#' Print/Summary Function for geomorph 
#' 
#' @param x print/summary object (from \code{\link{gpagen}})
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.gpagen <- function (x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call), fill=TRUE, "\n\n")
  cat("\nGeneralized Procrustes Analysis\n")
  cat("with Partial Procrustes Superimposition\n\n")
  cat(paste(x$p-x$nsliders-x$nsurf, "fixed landmarks\n"))
  cat(paste(x$nsliders+x$nsurf, "semilandmarks (sliders)\n"))
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
#' @param object print/summary object (from \code{\link{gpagen}})
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
#' @param x plot object (from \code{\link{gpagen}})
#' @param ... other arguments passed to plotAllSpecimens
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.gpagen <- function(x, ...){
  plotAllSpecimens(x$coords, ...)
}


## procD.lm

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object (from \code{\link{procD.lm}})
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords utilities
print.procD.lm<- function(x,...) {
  class(x) <- "lm.rrpp"
  print.lm.rrpp(x, ...)
}

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object (from \code{\link{procD.lm}})
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.procD.lm <- function(object,...){
  class(object) <- "lm.rrpp"
  effect.type <- object$ANOVA$effect.type
  anova.lm.rrpp(object, effect.type=effect.type, ...)
}

#' Plot Function for geomorph
#' 
#' @param x plot object (from \code{\link{procD.lm}})
#' @param type Indicates which type of plot, choosing among diagnostics,
#' regression, or principal component plots.  Diagnostic plots are similar to 
#' \code{\link{lm}} diagnostic plots, but for multivariate data.  Regression plots
#' plot multivariate dispersion in some fashion against predictor values. PC plots
#' project data onto the eigenvectors of the covariance matrix for fitted values.
#' @param outliers Logical argument to include outliers plot, if diagnostics
#' are performed
#' @param predictor An optional vector if "regression" plot type is chosen, 
#' and is a variable likely used in \code{\link{procD.lm}}.
#' This vector is a vector of covariate values equal to the number of observations.
#' @param reg.type If "regression" is chosen for plot type, this argument
#' indicates whether a prediction line 
#' (Predline) plot, or regression score (RegScore) plotting is performed.  
#' @param ... other arguments passed to plot (helpful to employ
#' different colors or symbols for different groups).  See
#' \code{\link{plot.default}} and \code{\link{par}}
#' @return An object of class "plot.procD.lm" is a list with components
#'  that can be used in other plot functions, such as the type of plot, points, 
#'  a group factor, and other information depending on the plot parameters used.
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.procD.lm <- function(x, type = c("diagnostics", "regression",
                                      "PC"), outliers=FALSE, predictor = NULL,
                          reg.type = c("PredLine", "RegScore"), ...){
  type <- match.arg(type)
  reg.type <- match.arg(reg.type)
  out <- plot.lm.rrpp(x, type = type,  predictor = predictor,
               reg.type = reg.type, ...)
  if(type == "diagnostics" && outliers) plotOutliers(x$GM$residuals)
  if(!is.null(x$GM)) out$GM <- x$GM
  class(out) <- "plot.procD.lm"
  invisible(out)
}

## morphol.disparity

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object (from \code{\link{morphol.disparity}})
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.morphol.disparity <- function (x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call), fill=TRUE, "\n\n")
  cat("\nRandomized Residual Permutation Procedure Used\n")
  cat(paste(x$permutations, "Permutations\n"))
  cat("\nProcrustes variances for defined groups")
  if(x$partial) cat("; partial variances (disparities) were calculated.\n") else cat("\n")
  print(x$Procrustes.var)
  if(x$partial) {
    propor <- x$Procrustes.var / sum(x$Procrustes.var)
    cat("\nProportion of total disparity for each group:\n")
    print(propor)
  }
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
#' @param object print/summary object (from \code{\link{morphol.disparity}})
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
#' @param x print/summary object (from \code{\link{phylo.integration}} or \code{\link{two.b.pls}})
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.pls <- function (x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call), fill=TRUE, "\n\n")
  
  cat("\nr-PLS:", round(x$r.pls, nchar(x$permutations)-1))
  cat("\n\nEffect Size (Z):", round(x$Z, nchar(x$permutations)))
  cat("\n\nP-value:", x$P.value)
  cat("\n\nBased on", x$permutations, "random permutations\n")
  if(!is.null(x$pairwise.Z)) {
    Z <- x$pairwise.Z
    P <- x$pairwise.P.values
    nms <- unique(unlist(strsplit(names(Z), "-")))
    m <- matrix(0, length(nms), length(nms))
    dimnames(m) <- list(nms, nms)
    dz <- dp <- as.dist(m)
    dz[1:length(Z)] <- Z
    dp[1:length(P)] <- P
    
    cat("\nPairwise statistics\n")
    cat("\nr-PLS:\n")
    print(x$r.pls.mat)
    cat("\nEffect Sizes (Z):\n")
    print(dz)
    cat("\nP-values::\n")
    print(dp)
  }
  
  invisible(x)
}

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object (from \code{\link{phylo.integration}} or \code{\link{two.b.pls}})
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.pls <- function(object, ...) {
  x<- object
  print.pls(x, ...)
}

#' Plot Function for geomorph
#' 
#' @param x plot object (from \code{\link{phylo.integration}} or \code{\link{two.b.pls}})
#' @param label Optional vector to label points
#' @param ... other arguments passed to plot. The function returns values that can be used with 
#' \code{\link{picknplot.shape}} (in a limited capacity). In most cases, greater 
#' flexibility can be attained with using \code{\link{plotRefToTarget}} and \code{\link{shape.predictor}}.
#' @return If shapes = TRUE, function returns a list containing the shape coordinates of the extreme ends of axis1 and axis2 
#' if 3D arrays were originally provided for each
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
#' 
plot.pls <- function(x, label = NULL, ...) {
  XScores <- x$XScores
  YScores <- x$YScores
  if(is.matrix(XScores)) XScores <- XScores[,1]
  if(is.matrix(YScores)) YScores <- YScores[,1]
  plot.args <- list(x = XScores, y = YScores,
                    main = "PLS1 Plot: Block 1 (X) vs. Block 2 (Y)", xlab = "PLS1 Block 1", 
                    ylab = "PLS1 Block 2", ...)
  do.call(plot, plot.args)
  pc <- prcomp(cbind(XScores, YScores))$x[,1]
  px <- predict(lm(XScores~pc))
  py <- predict(lm(YScores~pc))
  abline(lm(py ~ px), col = "red")
  if (length(label != 0)) {
    text(XScores, YScores, label, adj = c(-0.7, -0.7))
  }
  out <- list()
  out$plot.args <- plot.args
  out$A1 <- x$A1
  out$A2 <- x$A2
  class(out) <- "plot.pls"
  invisible(out)
}


## bilat.symmetry

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object (from \code{\link{bilat.symmetry}})
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
#' @param object print/summary object (from \code{\link{bilat.symmetry}})
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.bilat.symmetry <- function(object, ...) {
  x <- object
  print.bilat.symmetry(x, ...)
}

#' Plot Function for geomorph
#' 
#' @param x plot object (from \code{\link{bilat.symmetry}})
#' @param warpgrids Logical argument whether to include warpgrids
#' @param mesh Option to include mesh in warpgrids plots
#' @param ... other arguments passed to plot
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.bilat.symmetry <- function(x, warpgrids = TRUE, mesh= NULL, ...){
    k <- dim(x$symm.shape)[[2]]
    if(x$data.type == "Matching"){
      if(k==2){  
        par(mfrow=c(2,2),oma=c(1.5,0,1.5,0))
        plotAllSpecimens(x$symm.shape)
        plotAllSpecimens(x$asymm.shape)
        plotRefToTarget(x$DA.component[,,1],x$DA.component[,,2],method="TPS",main="Directional Asymmetry")
        plotRefToTarget(x$FA.component[,,1],x$FA.component[,,2],method="TPS",main="Fluctuating Asymmetry")
        mtext("Symmetric Shape Component (left) and Asymmetric Shape Component (right)",outer = TRUE,side=3)
        mtext("Mean directional (left) and fluctuating (right) asymmetry",side = 1, outer = TRUE)
        par(mfrow=c(1,1))
      }
      if(k==3){
        if (is.null(mesh)){
          open3d() ; mfrow3d(1, 2) 
          plotRefToTarget(x$DA.component[,,1],x$DA.component[,,2],method="points",main="Directional Asymmetry",box=FALSE, axes=FALSE)
          next3d()
          plotRefToTarget(x$FA.component[,,1],x$FA.component[,,2],method="points",main="Fluctuating Asymmetry",box=FALSE, axes=FALSE)
        } 
        if(!is.null(mesh)){
          open3d() ; mfrow3d(1, 2) 
          cat("\nWarping mesh\n")
          plotRefToTarget(x$DA.component[,,1],x$DA.component[,,2],mesh,method="surface")
          title3d(main="Directional Asymmetry")
          next3d()
          cat("\nWarping mesh\n")
          plotRefToTarget(x$FA.component[,,1],x$FA.component[,,2],mesh,method="surface")
          title3d(main="Fluctuating Asymmetry")
        }
      }
      layout(1) 
    }
    if(x$data.type == "Object"){
      if(warpgrids==TRUE){
        if(k==2){  
          par(mfrow=c(2,2),oma=c(1.5,0,1.5,0))
          plotAllSpecimens(x$symm.shape)
          plotAllSpecimens(x$asymm.shape)
          plotRefToTarget(x$DA.component[,,1],x$DA.component[,,2],method="TPS",main="Directional Asymmetry")
          plotRefToTarget(x$FA.component[,,1],x$FA.component[,,2],method="TPS",main="Fluctuating Asymmetry")
          mtext("Symmetric Shape Component (left) and Asymmetric Shape Component (right)",outer = TRUE,side=3)
          mtext("Mean directional (left) and fluctuating (right) asymmetry",side = 1, outer = TRUE)
        }
        if(k==3){
          if(is.null(mesh)) {
            open3d() ; mfrow3d(1, 2) 
            plotRefToTarget(x$DA.component[,,1],x$DA.component[,,2],method="points",main="Directional Asymmetry",box=FALSE, axes=FALSE)
            next3d()
            plotRefToTarget(x$FA.component[,,1],x$FA.component[,,2],method="points",main="Fluctuating Asymmetry",box=FALSE, axes=FALSE)
          } 
          if(!is.null(mesh)){
            open3d() ; mfrow3d(1, 2) 
            cat("\nWarping mesh\n")
            plotRefToTarget(x$DA.component[,,1],x$DA.component[,,2],mesh,method="surface")
            title3d(main="Directional Asymmetry")
            next3d()
            cat("\nWarping mesh\n")
            plotRefToTarget(x$FA.component[,,1],x$FA.component[,,2],mesh,method="surface")
            title3d(main="Fluctuating Asymmetry")
          }  
        }
        layout(1) 
      } 
    }
  }

## CR

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object (from \code{\link{phylo.modularity}})
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.CR <- function (x, ...) {
  cat("\nCall:\n")
  cat(deparse(x$call), fill=TRUE, "\n\n")
  cat(paste("\nCR:", round(x$CR, nchar(x$permutations))))
  cat(paste("\n\nP-value:", round(x$P.value, nchar(x$permutations))))
  cat(paste("\n\nEffect Size:", round(x$Z, nchar(x$permutations))))
  cat(paste("\n\nBased on", x$permutations, "random permutations"))
  if(!is.null(x$CInterval)) cat(paste("\n\nConfidence Intervals", round(x$CInterval,nchar(x$permutations)))) 
  invisible(x)
}


#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object (from \code{\link{phylo.modularity}})
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
#' @param x plot object (from \code{\link{phylo.modularity}})
#' @param ... other arguments passed to plot
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.CR <- function(x, ...){
  CR.val <- x$random.CR
  CR.obs <- x$CR
  options(scipen = 999)
  p <- x$P.value
  ndec <- nchar(x$permutations)
  CR.obs <- round(CR.obs, ndec)
  p <- round(p,ndec)
  main.txt <- paste("Observed CR =",CR.obs,";", "P-value =", p)
  hist(CR.val,30,freq=TRUE,col="gray",xlab="CR Coefficient",xlim=c(0,max(c(2,CR.val))),
       main=main.txt, cex.main=0.8)
  arrows(CR.obs,50,CR.obs,5,length=0.1,lwd=2)
  options(scipen = 0)
}

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object (from \code{\link{phylo.modularity}})
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
#' @param object print/summary object (from \code{\link{phylo.modularity}})
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
#' @param x plot object (from \code{\link{phylo.modularity}})
#' @param ... other arguments passed to plot
#' @export
#' @author Dean Adams
#' @keywords utilities
#' @keywords visualization
plot.CR.phylo <- function(x, ...){
  CR.val <- x$random.CR
  CR.obs <- x$CR
  options(scipen = 999)
  p <- x$P.value
  ndec <- nchar(x$permutations)
  p <- round(p,ndec)
  CR.obs <- round(CR.obs, ndec)
  main.txt <- paste("Observed CR =",CR.obs,";", "P-value =", p)
  hist(CR.val,30,freq=TRUE,col="gray",xlab="CR Coefficient",xlim=c(0,max(c(2,CR.val))),
       main=main.txt, cex.main=0.8)
  arrows(CR.obs,50,CR.obs,5,length=0.1,lwd=2)
  options(scipen = 0)
}

## physignal

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object (from \code{\link{physignal}})
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.physignal <- function(x, ...){
  cat("\nCall:\n")
  cat(deparse(x$call), fill=TRUE, "\n\n")
  cat(paste("\nObserved Phylogenetic Signal (K):", round(x$phy.signal, nchar(x$permutations))))
  cat(paste("\n\nP-value:", round(x$pvalue, nchar(x$permutations))))
  cat(paste("\n\nEffect Size:", round(x$Z, nchar(x$permutations))))
  cat(paste("\n\nBased on", x$permutations, "random permutations"))
  invisible(x)
}

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object (from \code{\link{physignal}})
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
#' @param x plot object (from \code{\link{physignal}})
#' @param ... other arguments passed to plot
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.physignal <- function(x, ...){
  K.val <- x$random.K
  K.obs <- x$phy.signal
  p <- x$pvalue
  ndec <- nchar(1 / x$permutations) - 2
  K.obs <- round(K.obs, ndec)
  p <- round(p, ndec)
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
  cat(paste("\n\nEffect Size:", round(x$Z, nchar(x$permutations))))
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

# compare.pls

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.compare.pls <- function(x,...){
  z <- x$sample.z
  z.pw <- x$pairwise.z
  p <- x$pairwise.P
  cat("\nEffect sizes\n\n")
  print(z)
  cat("\nEffect sizes for pairwise differences in PLS effect size\n\n")
  print(z.pw)
  cat("\nP-values\n\n")
  print(p)
  invisible(x)
}

# compare.CR

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.compare.CR<- function(object, ...) print.compare.CR(object,...)

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.compare.CR<- function(x,...){
  cat("\n", x$comment, "\n\n")
  z <- x$sample.z
  z.pw <- x$pairwise.z
  p <- x$pairwise.P
  cat("\nEffect sizes\n\n")
  print(z)
  cat("\nEffect sizes for pairwise differences in CR effect size\n\n")
  print(z.pw)
  cat("\nP-values\n\n")
  print(p)
  invisible(x)
}

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.compare.pls <- function(object, ...) print.compare.pls(object,...)

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' 
print.combined.set <- function(x,...){
  p <- x$points.by.set
  g <- length(p)
  cs <- x$CS
  css <- colMeans(cs)
  rcs <- css/sum(css)
  cat(paste("\nA total of", g, "subsets were combined\n\n"))
  y <- matrix(0, 3, g)
  if(!is.null(names(p))) colnames(y) <- names(p)
  y[1, ] <- p
  y[2, ] <- css
  y[3, ] <- rcs
  y <- as.data.frame(y)
  rownames(y) <- c("Number of points in subset", "Mean centroid size", "Mean relative size")
  if(x$norm.CS) rownames(y)[2] <- "Mean normalized centroid size"
  if(x$weighted) rownames(y)[2] <- "Mean weighted centroid size"
  print(y)
}

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
#'
summary.combined.set <- function(object, ...) print.combined.set(object, ...)

# mshape

#' Plot Function for geomorph
#' 
#' @param x plot object (from \code{\link{mshape}})
#' @param links An optional matrix defining for links between landmarks
#' @param ... other arguments passed to plot
#' @export
#' @author Antigoni Kaliontzopoulou
#' @keywords utilities
#' @keywords visualization
#' @seealso  \code{\link{define.links}}

plot.mshape <- function(x, links=NULL,...){
  if(any(is.na(x))) stop("Consensus configuration contains missing values. Please choose a different na.action in 'mshape'.")
  
  x <- as.matrix(x)
  class(x) <- "matrix"
  if(ncol(x)==2){
    x <- xy.coords(x)
    par(xpd=T)
    plot.new()
    plot.window(1.05*range(x$x), 1.05*range(x$y), asp = 1,
                xlab="", ylab="", xaxt="n", yaxt="n", bty="n",...)
    if(!is.null(links)){
      for (i in 1:nrow(links)){
        segments(x$x[links[i,1]], x$y[links[i,1]], 
                 x$x[links[i,2]], x$y[links[i,2]])
      }
    }
    plot.xy(x, type="p", cex=3, pch=21, bg="white")
    text(x, labels=1:length(x$x))
    par(xpd=F)
  } else {
  if(ncol(x)==3){
    plot3d(x, type="n", aspect=FALSE, xlab="", ylab="", zlab="", axes=F,...)
    if(!is.null(links)){
      for(i in 1:nrow(links)){
        segments3d(c(x[links[i,1], 1], x[links[i,2], 1]),
                   c(x[links[i,1], 2], x[links[i,2], 2]), 
                   c(x[links[i,1], 3], x[links[i,2], 3]))
      }
    }
    plot3d(x, add=T, type="s", col="white", alpha=0.25, shininess=2, fog=F)
    text3d(x, texts=1:nrow(x), cex=0.7, font=2)
  }
  }
}

# gm.prcomp

#' Print/Summary function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Antigoni Kaliontzopoulou
#' @keywords utilities
print.gm.prcomp <- function (x, ...) {
  class(x) <- "ordinate"
  xx <- summary(x)
  
  if(!is.null(x$ancestors)) {
    
    cat("\n\n")
    cat("Dispersion (variance) of points, after projection:\n")
      
    vars <- apply(x$x, 2, var)
    p <- vars/sum(vars)
    cp <- cumsum(vars)/sum(vars)
    tab <- rbind(vars, p, cp)
    rownames(tab) <- c("Tips Dispersion", "Proportion Tips Dispersion", "Cumulative Tips Dispersion")
    
    if(!is.null(x$anc.x)) {
      vars <- apply(x$anc.x, 2, var)
      p <- vars/sum(vars)
      cp <- cumsum(vars)/sum(vars)
      ancs <- rbind(vars, p, cp)
      rownames(ancs) <- c("Ancestors Dispersion", "Proportion Ancestors Dispersion", "Cumulative Ancestors Dispersion")
      tab <- rbind(tab, ancs)
    }
    print(tab)
    cat("\n\n")
    
  } else tab <- NULL
  
  out <- list(PC.summary = xx, Ancestor.summary = tab)
  invisible(out)
}

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Antigoni Kaliontzopoulou
#' @keywords utilities
summary.gm.prcomp <- function (object, ...) {
  print.gm.prcomp(object, ...)
}

#' Plot Function for geomorph
#' 
#' @param x An object of class \code{\link{gm.prcomp}}
#' @param axis1 A value indicating which PC axis should be displayed as the X-axis (default = PC1)
#' @param axis2 A value indicating which PC axis should be displayed as the Y-axis (default = PC2)
#' @param flip An argument that if not NULL can be used to flip components in the plot.  
#' The values need to match axis1 or axis2.  For example, if axis1 = 3 and axis2 = 4, flip = 1 will not
#' change either axis; flip = 3 will flip only the horizontal axis; flip = c(3, 4) will flip both axes.
#' @param phylo A logical value indicating whether the phylogeny should be projected to PC space
#' @param time.plot A logical value indicating if a 3D plot with the phylogeny and time as the 
#' z-axis is desired
#' @param phylo.par A list of plotting parameters for the inclusion of a phylogeny, including: logicals for 
#' whether features should be included (tip.labels, nodel.labels, anc.states), toggled as TRUE/FALSE; 
#' edge parameters (edge.color, edge.width, edge.lty); node parameters (node.bg, node.pch, node.cex);
#' and label parameters (tip.txt.cex, tip.txt.col, tip.txt.adj, node.txt.cex, node.txt.col, node.txt.adj).
#' @param ... other arguments passed to plot.  For plots with a phylogeny, these parameters pertain to 
#' the tip values.
#' @return An object of class "plot.gm.prcomp" is a list with components
#'  that can be used in other plot functions, such as the type of plot, points, 
#'  a group factor, and other information depending on the plot parameters used.  A time plot
#'  is an addendum to the normal 2D plot, and does not add additional output.
#'  
#'  NOTE: To visualize shape variation across PC axes in 2d plots, use \code{\link{picknplot.shape}}.
#'  
#'  
#' @export
#' @author Antigoni Kaliontzopoulou, Michael Collyer
#' @keywords utilities
#' @keywords visualization
#' @seealso  \code{\link{plotRefToTarget}} \code{\link{picknplot.shape}}


plot.gm.prcomp <- function(x, axis1 = 1, axis2 = 2, flip = NULL, phylo = FALSE, 
                           time.plot = FALSE, 
                           phylo.par = list(tip.labels = TRUE, 
                                            node.labels = TRUE, 
                                            anc.states = TRUE,
                                            node.pch = 21, 
                                            node.bg = "grey", 
                                            node.cex = 1, 
                                            edge.color = "black", 
                                            edge.width = 1,
                                            tip.txt.cex = 1, 
                                            tip.txt.col = "black", 
                                            tip.txt.adj = c(-0.1,-0.1),
                                            node.txt.cex = 1, 
                                            node.txt.col = "grey",
                                            node.txt.adj = c(-0.1, -0.1)), 
                           ...){

  class(x) <- "ordinate"
  pcdata <- as.matrix(x$x[, c(axis1, axis2)])
  Pcov <- x$Pcov
  xx <- plot(x, axis1 = axis1, axis2 = axis2, flip = flip, ...)
  plot.args <- xx$plot.args
  if(!is.null(plot.args$axes)) axes <- plot.args$axes else axes <- TRUE

  if(axes){
    abline(h = 0, lty=2)
    abline(v = 0, lty=2)
  }
  
  if(phylo || time.plot) {
    
    if(is.null(x$phy))
      stop("\nx must include a phylogeny for plotting\n", 
           call. = FALSE)
    
    if(!phylo && time.plot) {
      phylo = TRUE
    }
  }
  
  if(phylo) {
    
    p.p <- list(tip.labels = TRUE, node.labels = TRUE, anc.states = TRUE,
                            node.bg = "grey", node.pch = 21, node.cex = 1,
                            edge.color = "black", edge.width = 1,
                            tip.txt.cex = 1, tip.txt.col = "black", 
                            tip.txt.adj = c(-0.1, -0.1),
                            node.txt.cex = 1, node.txt.col = "grey",
                            node.txt.adj = c(-0.1, -0.1))
    
    m.p <- match(names(phylo.par), names(p.p))
    if(any(is.na(m.p)))
      stop("Some of the arguments in phylo.pars are different than those that are possible (see Arguments).\n",
           call. = FALSE)
    p.p[m.p] <- phylo.par
    
    phy <- x$phy
    tp <- add.tree(xx, phy, edge.col = p.p$edge.color,
                   edge.lwd = p.p$edge.width,
                   edge.lty = p.p$edge.lty, 
                   anc.pts = p.p$anc.states,
                   cex = p.p$node.cex,
                   pch = p.p$node.pch,
                   col = p.p$edge.color, 
                   bg = p.p$node.bg,
                   return.ancs = TRUE)
    
    if(p.p$tip.labels) {
      text(xx$points, rownames(xx$points), adj = p.p$tip.txt.adj,
           cex = p.p$tip.txt.cex, col = p.p$tip.txt.col)
    }
    
    if(p.p$node.labels) {
      text(tp, rownames(tp), adj = p.p$node.txt.adj,
           cex = p.p$node.txt.cex, col = p.p$node.txt.col)
    }
    
    phy.pcdata <- rbind(xx$points[phy$tip.label,], tp)
    N.tips <- length(phy$tip.label)
    
  }
  
  if(time.plot) {
    
    zaxis <- getNodeDepth(phy)
    zaxis <- abs(zaxis - max(zaxis))
    
    limits = function(x,s){ 
      r = range(x)
      rc = scale(r, scale=F)
      l = mean(r) + s*rc 
      l
    }
    
    t.p <- plot.args
    t.p$pch <- 19
    if(!is.null(p.p$bg)) t.p$col <- p.p$bg
    if(!is.null(t.p$bg)) t.p$col <- t.p$bg
    if(is.null(t.p$col)) t.p$col <- "black"
    if(is.null(t.p$cex)) t.p$cex <- 1
    
    view3d(phi=90, fov=30)
    plot3d(xx$points[,1], xx$points[,2], zaxis, type="n", 
           xlim = limits(phy.pcdata[,1], 1.5),
           ylim = limits(phy.pcdata[,2], 1.5),
           zlim = c(max(zaxis), min(zaxis)),
           size = 0.1,
           asp = c(1,1,1),
           xlab = xx$plot.args$xlab, 
           ylab = xx$plot.args$ylab, zlab = "Time")
    
    for (i in 1:nrow(phy$edge)) {
      lines3d(phy.pcdata[(phy$edge[i, ]), 1], phy.pcdata[(phy$edge[i, ]), 2], zaxis[(phy$edge[i, ])], 
              col = p.p$edge.color, lwd = p.p$edge.width)}
    
    points3d(xx$points[,1], xx$points[,2], zaxis[1:N.tips],
             col = t.p$col, size = t.p$cex*4)
    
    if(p.p$anc.states){
      points3d(tp[, 1], 
               tp[, 2], 
               zaxis[(N.tips + 1):nrow(phy.pcdata)], 
               col = p.p$node.bg, size = p.p$node.cex*4)
    }
    
    if(p.p$tip.labels){
      text3d(xx$points[,1], xx$points[,2], zaxis[1:N.tips], 
             rownames(xx$points),
             col = p.p$tip.txt.col, 
             cex = p.p$tip.txt.cex, 
             adj = p.p$tip.txt.adj) }
    
    if(p.p$node.labels){
      text3d(tp[, 1], tp[, 2],
             zaxis[(N.tips + 1):nrow(phy.pcdata)], 
             rownames(tp),
             col = p.p$node.txt.col, 
             cex = p.p$node.txt.cex, 
             adj = p.p$node.txt.adj)}
    
  }
  
  
  out <- list(PC.points = pcdata,   
              call = match.call())
  out$GM <- list()
  out$GM$A <- x$A
  
  out$plot.args <- plot.args
  out$Pcov <- Pcov
  if(phylo) {
    out$phylo <- list()
    out$phylo$phy <- phy
    out$phylo$phylo.par <- p.p
    out$phylo$phy.pcdata <- phy.pcdata
  }
  class(out) <- "plot.gm.prcomp" 
  
  invisible(out)
  
}


# geomorphShapes

#' Print/Summary function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.geomorphShapes <- function (x, ...) {
  cat("\nLandmark data for", x$n, "specimens")
  cat("\nNumber of fixed landmarks:", length(x$fixed))
  cat("\nNumber of (sliding) semilandmarks:", length(x$sliders))
  cat("\nNumber of total landmarks:", x$p)
  if(x$scaled) cat("\nlandmarks have been scaled.") else
    cat("\nLandmarks have not been scaled.")
  cat("\n\nThis information is based on information available from a class 'shapes' object,")
  cat("\nas it was originally read.  (It might have been since edited but this summary will not change.)")
  cat("\nThe curves matrix (for use in gpagen) is also based on the same information.")
  cat("\nThis matrix can be modified to alter which landmarks are semilandmarks.\n")
  invisible(x)
}

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.geomorphShapes <- function (object, ...) {
  print.geomorphShapes(object, ...)
}
