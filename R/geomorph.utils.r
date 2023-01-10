## S3 GENERIC FUNCTIONS

## gpagen

#' Print/Summary Function for geomorph 
#' 
#' @param x print/summary object (from \code{\link{gpagen}})
#' @param ... other arguments passed to print/summary
#' @method print gpagen
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
#' @method summary gpagen
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
#' @method plot gpagen
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
#' @method print procD.lm
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
#' @method summary procD.lm
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
#'@method plot procD.lm
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
#' @method print morphol.disparity
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
#' @method summary morphol.disparity
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
#' @method print pls
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
#' @method summary pls
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
#' @method plot pls
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
  out$Pcov <- x$Pcov
  class(out) <- "plot.pls"
  invisible(out)
}


## bilat.symmetry

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object (from \code{\link{bilat.symmetry}})
#' @param ... other arguments passed to print/summary
#' @method print bilat.symmetry
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
#' @method summary bilat.symmetry
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
#' @method plot bilat.symmetry
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
#' @method print CR
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
#' @method summary CR
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
#' @method plot CR
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
#' @method print CR.phylo
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
#' @method summary CR.phylo
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
#' @method plot CR.phylo
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

## physignal abd physignal.z

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object (from \code{\link{physignal}})
#' @param ... other arguments passed to print/summary
#' @method print physignal
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.physignal <- function(x, ...){
  cat("\nCall:\n")
  cat(deparse(x$call), fill=TRUE, "\n\n")
  cat(paste("\nObserved Phylogenetic Signal (K):", round(x$phy.signal, nchar(x$permutations))))
  cat(paste("\n\nP-value:", round(x$pvalue, nchar(x$permutations))))
  cat(paste("\n\nBased on", x$permutations, "random permutations"))
  cat("\n\n Use physignal.z to estimate effect size.")
  invisible(x)
}

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object (from \code{\link{physignal.z}})
#' @param ... other arguments passed to print/summary
#' @method print physignal.z
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.physignal.z <- function(x, ...){
  cat("\nCall:\n")
  cat(deparse(x$call), fill=TRUE, "\n")
  
  cat("Evaluation of phylogenetic signal effect size\n")
  cat("Optimization method:", x$opt.method)
  cat("\nOptimization performed in", x$opt.dim, "data dimensions.\n")
  
  if(is.na(x$Z)) {
    cat("The scaling parameter, lambda, was optimized at 0.")
    cat("\nThis means that the log-likehood was invariant across permutations")
    cat("\nand there is no phylogenetic signal in the data.\n\n")
  } else {
    cat(paste("\nObserved phylogenetic signal effect size (Z):", round(x$Z, nchar(x$permutations))))
    cat(paste("\n\nP-value:", round(x$pvalue, nchar(x$permutations))), "based on", x$permutations, "random permutations")
    cat("\n\nFor a model with a log-likelihood of", round(x$rand.logL[[1]], nchar(x$permutations)))
    cat("\na branch-scaling (lambda) of", round(x$lambda, nchar(x$permutations)))
    cat("\nand a ratio of Bownian Motion fit (K) of", round(x$K, nchar(x$permutations)), "\n")
    
    if(!is.null(x$K.by.p)) {
      cat("\nK measured across phylogenetically-aligned components (1, 1:2, 1:3, ...\n")
      print(x$K.by.p)
      cat("\nLambda measured across phylogenetically-aligned components (1, 1:2, 1:3, ...\n")
      print(x$lambda.by.p)
      cat("\nLog-likelihood measured across phylogenetically-aligned components (1, 1:2, 1:3, ...\n")
      print(x$logL.by.p)
    }
  }
    
  invisible(x)
}

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object (from \code{\link{physignal}})
#' @param ... other arguments passed to print/summary
#' @method summary physignal
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.physignal <- function(object, ...) {
  x <- object
  print.physignal(x, ...)
}

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object (from \code{\link{physignal.z}})
#' @param ... other arguments passed to print/summary
#' @method summary physignal.z
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.physignal.z <- function(object, ...) {
  x <- object
  print.physignal.z(x, ...)
}

#' Plot Function for geomorph
#' 
#' @param x plot object (from \code{\link{physignal}})
#' @param ... other arguments passed to plot
#' @method plot physignal
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

#' Plot Function for geomorph
#' 
#' @param x plot object (from \code{\link{physignal.z}})
#' @param ... other arguments passed to plot
#' @method plot physignal.z
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.physignal.z <- function(x, ...){
  
  if(is.na(x$Z)) {
    cat("The scaling parameter, lambda, was optimized at 0.")
    cat("\nThis means that the log-likehood was invariant across permutations,")
    cat("\nthere is no phylogenetic signal in the data, and there is nothing to plot.\n\n")
  } else {
    
    ll <- x$rand.logL
    z <- box.cox(ll)$transformed
    z[1] <- x$Z
    
    opar <- par()$mfrow
    par(mfrow = c(1,2))
    p <- x$pvalue
    ndec <- nchar(1 / x$permutations) - 2
    ll.obs <- round(ll[[1]], ndec)
    Z.obs <- round(z[[1]], ndec)
    p <- round(p, ndec)
    main.txt <- paste("Observed log-likelihood =", ll.obs,";", "P-value =", p)
    hist(ll, 30, freq=TRUE, col = "gray", 
         xlab = paste("log-liklihoods:", x$permutations, "RRPP premutations"),
         main = main.txt, cex.main = 0.8)
    arrows(ll.obs, 50, ll.obs, 5, length = 0.1, lwd=2)
    
    main.txt <- paste("Standardized effect size =", Z.obs)
    hist(z, 30, freq=TRUE, col = "gray", 
         xlab = "Standardized values",
         main = main.txt, cex.main = 0.8)
    arrows(Z.obs, 50, Z.obs, 5, length = 0.1, lwd=2)
    
    par(mfrow = opar)
  }
}


## evolrate

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @method print evolrate
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
#' @method summary evolrate
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
#' @method print evolrate1
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
#' @method summary evolrate1
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
#' @method plot evolrate
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
#' @method print compare.pls
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
#' @method summary compare.CR
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.compare.CR<- function(object, ...) print.compare.CR(object,...)

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @method print compare.CR
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
#' @method summary compare.pls
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.compare.pls <- function(object, ...) print.compare.pls(object,...)

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @method print compare.physignal.z
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.compare.physignal.z<- function(x,...){
  z <- x$sample.z
  z.pw <- x$pairwise.z
  p <- x$pairwise.P
  cat("\nEffect sizes (Z-scores)\n\n")
  print(z)
  cat("\nEffect sizes for pairwise differences in phylogenetic signal.\n\n")
  print(z.pw)
  cat("\nP-values\n\n")
  print(p)
  invisible(x)
}

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @method summary compare.physignal.z
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.compare.physignal.z <- function(object, ...) 
  print.compare.physignal.z(object,...)


#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @method print combined.set
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
#' @method summary combined.set
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
#' @method plot mshape
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
#' @method print gm.prcomp
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
#' @method summary gm.prcomp
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
#' @method plot gm.prcomp
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
#' @method print geomorphShapes
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
#' @method summary geomorphShapes
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.geomorphShapes <- function (object, ...) {
  print.geomorphShapes(object, ...)
}

# module.eigen

#' Print/Summary function for geomorph
#' 
#' @param x An object of class \code{\link{module.eigen}}
#' @param ... other arguments passed to print/summary
#' @method print module.eigen
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.module.eigen <- function(x, ...) {
  cat("\nEigen-analysis of covariance matrices for", x$n.modules, "modules\n")
}


#' Print/Summary function for geomorph
#' 
#' @param object An object of class \code{\link{module.eigen}}
#' @param use.rel.dims A logical argument for whether to use only the relevant dimensions in the summary, 
#' based on a broken stick model.
#' @param PC A numeric value to indicate for which principal component to perform vector correlations
#' @param ... other arguments passed to print/summary
#' @method summary module.eigen
#' @export
#' @author Michael Collyer
#' @return An object of class "summary.module.eigen" is a list with components
#' including table of eigenvalues, table of proportion of variation explained,
#' a pairwise table of vector correlations (for indicated principal component), and a pairwise
#' table of Krzanowski squared correlations.
#' @keywords utilities
summary.module.eigen <- function(object, use.rel.dims = TRUE, PC = 1) {
  
  m <- if(use.rel.dims) 1:max(unlist(object$rel.dims)) else 
    1:length(object$evals$total)
  
  if(PC > max(m)) {
    cat("\nPC is a number greater than the number of relevant dimensions.  Setting it to PC = 1\n")
    PC <- 1
  }
  
  if(!is.null(object$evecs)) {
    v0 <- object$evecs$ind[,m]
    vM <- object$evecs$mod[,m]
    vI <- object$evecs$int[,m]
    vT <- object$evecs$total[,m]
    
    rp <- abs(vec.cor.matrix(rbind(Independence= v0[, PC],
                                   Modularity = vM[, PC], 
                                   Integration = vI[, PC], 
                                   Total = vT[, PC])))
    
    krzCor <- function(v1, v2){
      m <- min(ncol(v1), ncol(v2))
      sum(crossprod(v1[, 1:m], v2[, 1:m])^2) / m
    }
    
    rK <- as.dist(diag(4))
    rK[1:6] <- c(krzCor(vM, v0), krzCor(v0, vI), krzCor(v0, vT),
                 krzCor(vM, vI), krzCor(vM, vT), krzCor(vI, vT))
    rK <- as.matrix(rK)
    diag(rK) <- 1
    dimnames(rK) <- dimnames(rp)
    
  } else rp <- rK <- NULL
  
  
  rel.dims <- if(use.rel.dims) object$rel.dims[c("mod", "int", "total")] else lapply(1:3, function(.) m)
  rel.eigs <- object$evals[c("mod", "int", "total")]
  rel.eigs <- Map(function(e, d) e[d], rel.eigs, rel.dims)
  rel.eigs <- lapply(rel.eigs, function(x){
    y <- matrix(NA, max(m))
    y[1:length(x)] <- x
    y
  })
  rel.eigs <- c(list(ind = matrix(object$evals$ind[m])), rel.eigs)
  eig.tab <- as.data.frame(do.call(cbind, rel.eigs))
  colnames(eig.tab) <- c("Independence", "Modularity", "Integration", "Total")
  rownames(eig.tab) <- paste("Comp", m)
  
  prop.tab <- eig.tab / sum(object$evals$total)
  

  out <- list(rPC.mat = rp, rKrz.mat = rK,
              eig.tab = eig.tab, prop.tab = prop.tab,
              use.rel.dims = use.rel.dims,
              PC = PC, total.dims = length(object$eval$total),
              prop.mod.cells = object$prop.mod.cells,
              prop.int.cells = object$prop.int.cells)
  
  class(out) <- "summary.module.eigen"
  out
}


#' Print/Summary function for geomorph
#' 
#' @param x An object of class \code{\link{summary.module.eigen}}
#' @param ... other arguments passed to print/summary
#' @method print summary.module.eigen
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.summary.module.eigen <- function(x, ...) {
  
  cat("\nPercentage of covariance matrix represented by modularity covariances:", 
      round(x$prop.mod.cells * 100, 2), "percent.\n")
  cat("\nPercentage of covariance matrix represented by integration covariances:", 
      round(x$prop.int.cells * 100, 2), "percent.\n")
  
  if(x$use.rel.dims) {
    cat("\nDisplaying results for only", nrow(x$eig.tab), "relevant dimensions, of", x$total.dims, "dimensions, total,\n")
    cat("based on a broken stick model (from the distribution of 'independence' eigen values).\n")
  } else {
    cat("\nDisplaying results for all components.\n")
  }
  
  cat("\nEigen values:\n")
  print(x$eig.tab)
  cat("\nProportion of variance explained:\n")
  print(x$prop.tab)
  
  if(!is.null(x$rPC.mat)) {
    
    cat("\nPairwise vector correlations between PC", x$PC, "\n")
    print(as.dist(x$rPC.mat))
    
    cat("\nKrzanowski (squared) vector correlations, using defined components\n")
    print(as.dist(x$rKrz.mat))
    
  }
    
}

#' Plot function for geomorph
#' 
#' @param x An object of class \code{\link{module.eigen}}
#' @param use.rel.dims A logical argument for whether to use only the relevant dimensions in the plot, 
#' based on a broken stick model.
#' @param ... other arguments passed to plot.  
#' 
#' This function vectorizes a summary table (for proportions of variance) 
#' from \code{\link{summary.module.eigen}}.  For example, if a 5 x 4 table is produced, a vector of length 20 is produced
#' by the columns of the table.  Some plot parameters can be manipulated, but one has to do it with respect to all values.  
#' For example, for the 5 x 4 table, bg = c(rep(1:4, each = 5)) produces a vector of length 20 to describe point 
#' background colors.
#' @method plot module.eigen
#' @export
#' @author Michael Collyer
#' @keywords utilities
plot.module.eigen <- function(x, use.rel.dims = TRUE, 
                             add.lines = TRUE, add.legend = TRUE,
                             ...) {
  s <- summary.module.eigen(x, use.rel.dims = use.rel.dims)
  df <- s$prop.tab
  y <- as.vector(as.matrix(df))
  x <- rep(1:nrow(df), 4)
  
  plot.args <- list(...)
  plot.args$x <- x
  plot.args$y <- y
  if(is.null(plot.args$pch)) plot.args$pch <- rep(19, length(y))
  if(is.null(plot.args$col)) plot.args$col <- rep(c(8, 4, 6, 1), each = nrow(df))
  if(is.null(plot.args$bg)) plot.args$bg <- rep(c(8, 4, 6, 1), each = nrow(df))
  
  if(is.null(plot.args$xlab)) plot.args$xlab = "Component"
  if(is.null(plot.args$ylab)) plot.args$ylab = "Proportion of total variance"

  do.call(plot, plot.args)
  
  if(add.lines) {
    pt.lty <- plot.args$lty
    if(is.null(pt.lty)) pt.lty <- 1
    pt.lwd <- plot.args$lwd
    if(is.null(pt.lwd)) pt.lwd<- 1
    points(1:nrow(df), df[,1], type = "l", lty = pt.lty, lwd = pt.lwd,
           col = plot.args$col[1])
    points(1:nrow(df), df[,2], type = "l", lty = pt.lty, lwd = pt.lwd,
           col = plot.args$col[nrow(df) + 1])
    points(1:nrow(df), df[,3], type = "l", lty = pt.lty, lwd = pt.lwd,
           col = plot.args$col[2* nrow(df) + 1])
    points(1:nrow(df), df[,4], type = "l", lty = pt.lty, lwd = pt.lwd,
           col = plot.args$col[3* nrow(df) + 1])
    
    plot.args$type = "p"
    do.call(points, plot.args)
  } else pt.lty <- pt.lwd <- NULL
  
  if(add.legend) {
    pt.pch <- plot.args$pch[c(1, nrow(df) + 1, 2* nrow(df) + 1, 3* nrow(df) + 1)]
    pt.col<- plot.args$col[c(1, nrow(df) + 1, 2* nrow(df) + 1, 3* nrow(df) + 1)]
    pt.bg <- plot.args$bg[c(1, nrow(df) + 1, 2* nrow(df) + 1, 3* nrow(df) + 1)]
    
    if(add.lines) 
      legend("topright",
           c("Independence", "Modularity", "Integration", "Total"),
           lty = pt.lty, lwd= pt.lwd, 
           pch = pt.pch, col = pt.col, pt.bg = pt.bg) else 
             legend("topright",
                    c("Independence", "Modularity", "Integration", "Total"),
                    pch = pt.pch, col = pt.col, pt.bg = pt.bg)
  }
    
}

# K.modules

#' Print/Summary function for geomorph
#' 
#' @param x An object of class \code{\link{K.modules}}
#' @param ... other arguments passed to print/summary
#' @method print K.modules
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.K.modules <- function(x, ...) {
  M <- x$modules[[1]]
  nlev <- nlevels(as.factor(M))
  nsims <- length(x$modules)
  
  cat("\nK.modules analysis performed for", nlev, "modules, with", nsims, "simulations\n")
  
}

#' Print/Summary function for geomorph
#' 
#' @param object An object of class \code{\link{K.modules}}
#' @param ... other arguments passed to print/summary
#' @method summary K.modules
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.K.modules <- function(object, ...) {
  print.K.modules(object)
  
  if(object$hypothesis){
    
    if(object$rel.dims > 1){
      
      cat("An a priori hypothesis was defined as:\n")
      cat(object$modules[[object$hypothesis.rank]])
      cat("\n\nA total of", object$rel.dims, "relevant dimensions were used for evaluation.\n")
      
      cat("This hypothesis ranked", object$hypothesis.rank, "among the simulations\n")
      cat("(A percentile of", round(object$hypothesis.rank/length(object$modules), 3) * 100, "percent.)\n")
      cat("\n")
      cat("The sum of the eigenvalues in the", object$rel.dims, "relevant dimensions for this hypothesis was:",
          object$eig.sums[1], ",")
      cat("which is", round(object$eig.sums[[object$hypothesis.rank]] / object$eig1.ref *100, 2), 
          "% of the sum of the corresponding eigenvalues from the total covariance.")
      cat("\n\nA summary of all simulated eigenvalue sums is presented below (see plot.K.modules for plotting tools)\n\n")
      summary(object$eig.sums)
      
    } else {
      
      cat("An a priori hypothesis was defined as:\n")
      cat(object$modules[[object$hypothesis.rank]])
      cat("\n\nOne relevant dimension was used for evaluation.\n")
      
      cat("This hypothesis ranked", object$hypothesis.rank, "among the simulations\n")
      cat("(A percentile of", round(object$hypothesis.rank/length(object$modules), 3) * 100, "percent.)\n")
      cat("\n")
      cat("The eigenvalue was:", object$eig.sums[1], ", which is", round(object$eig.sums[[object$hypothesis.rank]] / object$eig1.ref *100, 2), 
          "% of the maximum value from the total covariance matrix.")
      cat("\n\nA summary of all simulated eigenvalues is presented below (see plot.K.modules for plotting tools)\n\n")
      summary(object$eig.sums)
    }
    
   
  }
  
}


#' Plot function for geomorph
#' 
#' @param x An object of class \code{\link{K.modules}}
#' @param modules A vector or integer for which modules to consider in plot (should generally be small in length
#' if results are to be tractable).  If a profile plot is chosen, all eigenvalues will be plotted.  (The default value
#' is the highest ranking 6 modules.)
#' @param type One of "profile", config", or "pcoa" for plotting a profile of eigenvalues, 
#' landmark configurations with modules formatted differently, 
#' or a principal coordinate plot, respectively.
#' @param relativize A logical argument for whether to make eigenvalues a fraction of the maximum possible
#' (for profile plots).
#' @param ... other arguments passed to plot, especially for changing points in a configuration.
#' 
#' This function provides plotting support for \code{\link{K.modules}}.  One of two plots can be produced: (1) a 
#' canvas of landmark configuration with points formatted differently for modules or (2) a principal coordinate plot, based
#' on Riemannian distances among modular covariance matrices.  Based on the modules requested (output from a K.modules
#' object), several configuration plots can be generated or one single principal coordinate plot with points equal in 
#' number to the vector length for "modules".  The length of this vector should be reasonable to avoid long computation of 
#' Riemannian distance among modular covariance matrices, or to avoid generating many configurations.  See 
#' \code{\link{K.modules}} for plotting examples.
#' @method plot K.modules
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' 
plot.K.modules <- function(x, modules = 1:6, 
                           type = c("profile", "config", "pcoa"),
                           relativize = TRUE, ...) {
  type <- match.arg(type)
  M <- x$modules[modules]
  plot.args <- list(...)
  
  if(type == "profile") {
    y <- if(relativize) as.vector(x$eig.sums) / x$eig1.ref else 
      as.vector(x$eig.sums)
    plot.args$x <- 1:length(y)
    plot.args$y<- y
    plot.args$type <- "l"
    if(is.null(plot.args$xlab)) {
      plot.args$xlab <- if(x$rel.dims == 1) "Rank of first eigenvalue" else "Rank of summed eigenvalues"
    }
    if(is.null(plot.args$ylab)) {
      if(x$rel.dims == 1) {
        
        plot.args$ylab <- if(relativize)  "Eigenvalue (relativized)" else 
          "Eigenvalue"
        
      } else {
        
        plot.args$ylab <- if(relativize)  "Summed eigenvalues (relativized)" else 
          "Summed eigenvalues"
      }
    }
    
    do.call(plot, plot.args)
    
    if(x$hypothesis) {
      point.args <- plot.args
      point.args$x <- x$hypothesis.rank
      point.args$y <- y[x$hypothesis.rank]
      point.args$type <- "p"
      if(is.null(point.args$pch)) point.args$pch <- 19
      if(is.null(point.args$cex)) point.args$cex <- 1.2
      do.call(points, point.args)
      if(is.null(plot.args$main)) {
        if(x$rel.dims == 1) 
          title("Modular eigenvalue profile, with described hypothesis indicated", 
              cex.main = 0.7)
      if(x$rel.dims > 1) 
        title("Modular summed eigenvalues profile, with described hypothesis indicated", 
              cex.main = 0.7)
      }
      
    } else {
      if(is.null(plot.args$main)) {
        if(x$rel.dims == 1) 
          title("Modular eigenvalue profile")
        if(x$rel.dims > 1) 
          title("Modular summed eigenvalues profile")
        
      }
        
    }
  }
  
  if(type == "config"){
    
    opar <- par()
    
    if(is.null(x$mean))
      stop("Coordinate data not used in K.modules, so no configuration is possible.\n",
           call. = FALSE)
    mn <- x$mean
    class(mn) <- "matrix"
    
    if(length(modules) > 16){
      warning("Too many modules were chosen for effective plotting. Only the first 16 will shown.\n",
           immediate. = TRUE)
      M <- modules[1:16]
    }
    
    m <- length(M)
    crow <- ceiling(sqrt(m))
    ccol <- ceiling(m / crow)
    cdims <- c(crow, ccol)
    
    dims <- dim(mn)
    plot.args$x <- mn[,1]
    plot.args$y <- mn[,2]
    plot.args$xlab <- "x"
    plot.args$ylab <- "y"
    plot.args$asp = 1
    
    if(dims[2] == 3) {
      plot.args$z <- mn[,3]
      plot.args$zlab <- "z"
    }
    
    if(dims[2] == 3) mfrow3d(crow, ccol) else
      par(mfrow = cdims)
    
    for(i in 1:m){
      plot.args$col <- plot.args$bg <- M[[i]]
      if(dims[2] == 3) {
        do.call(plot3d, plot.args)
        aspect3d("iso")
        title3d(names(M)[[i]])
      } else {
        do.call(plot, plot.args)
        title(names(M)[[i]])
      }
    }
    par(opar)
  }
  
  if(type == "pcoa"){
    
    cat("\nThis can be a computationally intensive procedure.\n")
    if(length(modules) > 100) {
      modules <- modules[1:100]
      cat("Truncating the number of modules to the first 100 requested,\n")
      cat("in order to have a reasonable result.\n")
    }
    
    M <- if(length(x$modules) < length(modules)) x$modules else
       x$modules[modules]
    
    Dist <- function(C){
      d <-lapply(1:(length(C) - 1), function(j){
        xx <- C[[j]]
        yy <- C[-(1:j)]
        sapply(1:length(yy), function(jj){
          d <- Re(eigen(fast.solve(xx) %*% yy[[jj]])$values)
          d <- d[d>0]
          sqrt(sum(log(d)^2))
        })
      })
      d <- unlist(d)
      Dmat <- dist(matrix(0, length(C), length(C)))
      for(i in 1:length(Dmat)) Dmat[[i]] <- d[i]
      as.matrix(Dmat)
    }
    
    K <- length(unique(M[[1]]))
    Covs <- lapply(M, function(y){
      KM <- K.modules(x$A, K = K, hyp = y, nsims = 1)
      cov <- KM$VCV
      k <- ncol(cov) / length(y)
      keep <- rep(y, each = k)
      gps <- levels(as.factor(y))
      res <- matrix(0, nrow(cov), ncol(cov))
      for(i in 1:length(gps)) {
        found <- which(keep == gps[i])
        res[found, found] <- cov[found, found]
      }
     res 
    })
    
    names(Covs) <- names(M)
    d <- Dist(Covs)
    pcoa <- as.data.frame(cmdscale(d))
    names(pcoa) <- paste("C", 1:ncol(pcoa), sep = "")
    plot(pcoa, ...)
    text(pcoa, names(M), pos=3, ...)
  }
  
}


#' Print/Summary function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @method print compare.ZVrel
#' @export
#' @author Dean Adams
#' @keywords utilities
print.compare.ZVrel <- function(x,...){
  z <- x$sample.Z.obs
  z.pw <- x$pairwise.z
  p <- x$pairwise.P
  cat("\nEffect sizes\n\n")
  print(z)
  cat("\nEffect sizes for pairwise differences in rel.eig effect size\n\n")
  print(z.pw)
  cat("\nP-values\n\n")
  print(p)
  invisible(x)
}

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @method summary compare.ZVrel
#' @export
#' @author Dean Adams
#' @keywords utilities
summary.compare.ZVrel <- function(object, ...) {
  print.compare.ZVrel(object,...)
}


## geomorph.data.frame

#' Handle missing values in rrpp.data.frame objects
#'
#' @param object object (from \code{\link{geomorph.data.frame}})
#' @param ... further arguments (currently not used)
#' @method na.omit geomorph.data.frame
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @examples
#' data(plethspecies)
#' Y.gpa <- gpagen(plethspecies$land, verbose = TRUE)
#' gdf <- geomorph.data.frame(Y.gpa)
#' gdf$d <- Y.gpa$procD
#' gdf$group <- c(rep(1, 4), rep(2, 4), NA) # one unknown group designation
#' gdf
#' ndf <- na.omit(gdf)
#' ndf

na.omit.geomorph.data.frame <- function(object, ...) {
  class(object) <- "rrpp.data.frame"
  out <- na.omit(object)
  class(out) <- "geomorph.data.frame"
  out
}

