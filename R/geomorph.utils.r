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
  if(any(is.na(match(effect.type, 
                     c("SS", "MS", "Rsq", "F", "cohen")))))
    effect.type = "F"
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
  if(x$method=="RV") {
    cat(paste("\nRV:", round(x$RV, nchar(x$permutations)-1)))
    cat(paste("\n\nP-value:", round(x$P.value, nchar(x$permutations)-1)))
    cat(paste("\n\nBased on", x$permutations, "random permutations\n"))
  }
  if(x$method=="PLS") {
    cat(paste("\nr-PLS:", round(x$r.pls, nchar(x$permutations)-1)))
    cat(paste("\n\nP-value:", round(x$P.value, nchar(x$permutations)-1)))
    cat(paste("\n\nBased on", x$permutations, "random permutations\n"))
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
        plotRefToTarget(x$DA.mns[,,1],x$DA.mns[,,2],method="TPS",main="Directional Asymmetry")
        plotRefToTarget(x$FA.mns[,,1],x$FA.mns[,,2],method="TPS",main="Fluctuating Asymmetry")
        mtext("Symmetric Shape Component (left) and Asymmetric Shape Component (right)",outer = TRUE,side=3)
        mtext("Mean directional (left) and fluctuating (right) asymmetry",side = 1, outer = TRUE)
        par(mfrow=c(1,1))
      }
      if (k==3){
        if (is.null(mesh)){
          open3d() ; mfrow3d(1, 2) 
          plotRefToTarget(x$DA.mns[,,1],x$DA.mns[,,2],method="points",main="Directional Asymmetry",box=FALSE, axes=FALSE)
          next3d()
          plotRefToTarget(x$FA.mns[,,1],x$FA.mns[,,2],method="points",main="Fluctuating Asymmetry",box=FALSE, axes=FALSE)
        } 
        if(!is.null(mesh)){
          open3d() ; mfrow3d(1, 2) 
          cat("\nWarping mesh\n")
          plotRefToTarget(x$DA.mns[,,1],x$DA.mns[,,2],mesh,method="surface")
          title3d(main="Directional Asymmetry")
          next3d()
          cat("\nWarping mesh\n")
          plotRefToTarget(x$FA.mns[,,1],x$FA.mns[,,2],mesh,method="surface")
          title3d(main="Fluctuating Asymmetry")
        }
      }
      layout(1) 
    }
    if(x$data.typ == "Object"){
      if(warpgrids==TRUE){
        if(k==2){  
          par(mfrow=c(2,2),oma=c(1.5,0,1.5,0))
          plotAllSpecimens(x$symm.shape)
          plotAllSpecimens(x$asymm.shape)
          plotRefToTarget(x$DA.mns[,,1],x$DA.mns[,,2],method="TPS",main="Directional Asymmetry")
          plotRefToTarget(x$FA.mns[,,1],x$FA.mns[,,2],method="TPS",main="Fluctuating Asymmetry")
          mtext("Symmetric Shape Component (left) and Asymmetric Shape Component (right)",outer = TRUE,side=3)
          mtext("Mean directional (left) and fluctuating (right) asymmetry",side = 1, outer = TRUE)
        }
        if (k==3){
          if(is.null(mesh)) {
            open3d() ; mfrow3d(1, 2) 
            plotRefToTarget(x$DA.mns[,,1],x$DA.mns[,,2],method="points",main="Directional Asymmetry",box=FALSE, axes=FALSE)
            next3d()
            plotRefToTarget(x$FA.mns[,,1],x$FA.mns[,,2],method="points",main="Fluctuating Asymmetry",box=FALSE, axes=FALSE)
          } 
          if(!is.null(mesh)){
            open3d() ; mfrow3d(1, 2) 
            cat("\nWarping mesh\n")
            plotRefToTarget(x$DA.mns[,,1],x$DA.mns[,,2],mesh,method="surface")
            title3d(main="Directional Asymmetry")
            next3d()
            cat("\nWarping mesh\n")
            plotRefToTarget(x$FA.mns[,,1],x$FA.mns[,,2],mesh,method="surface")
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
  p <- x$P.value
  ndec <- nchar(x$permutations)
  CR.obs <- round(CR.obs, ndec)
  main.txt <- paste("Observed CR =",CR.obs,";", "P-value =", p)
  hist(CR.val,30,freq=TRUE,col="gray",xlab="CR Coefficient",xlim=c(0,max(c(2,CR.val))),
       main=main.txt, cex.main=0.8)
  arrows(CR.obs,50,CR.obs,5,length=0.1,lwd=2)
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
  p <- x$P.value
  ndec <- nchar(x$permutations)
  CR.obs <- round(CR.obs, ndec)
  main.txt <- paste("Observed CR =",CR.obs,";", "P-value =", p)
  hist(CR.val,30,freq=TRUE,col="gray",xlab="CR Coefficient",xlim=c(0,max(c(2,CR.val))),
       main=main.txt, cex.main=0.8)
  arrows(CR.obs,50,CR.obs,5,length=0.1,lwd=2)
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


# plotTangentSpace

#' Print/Summary Function for geomorph
#' 
#' @param x print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.plotTangentSpace <- function (x, ...) {
  cat("\nPC Summary\n\n")
  print(x$pc.summary)
  invisible(x)
}

#' Print/Summary Function for geomorph
#' 
#' @param object print/summary object
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.plotTangentSpace <- function (object, ...) {
  print.plotTangentSpace(object, ...)
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
  sum.tab <- function(x) {
    var.PCs <- apply(x$x, 2, var)
    y <- rbind(x[["d"]], var.PCs/sum(var.PCs), cumsum(var.PCs/sum(var.PCs)))
    
    colnames(y) <- paste("PC", 1:ncol(y), sep="")
    rownames(y) <- c("Eigenvalues", "Proportion of variance", "Cumulative Proportion")
    y
  }
    tab.list <- sum.tab(x)
    cat("Importance of components:", "\n")
    print(tab.list); cat("\n")
    invisible(tab.list)
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
#' @param phylo A logical value indicating whether the phylogeny should be projected to PC space
#' @param phylo.par A list of plotting parameters for the phylogeny edges (edge.color, edge.width, edge.lty)
#' and nodes (node.bg, node.pch, node.cex)
#' @param ... other arguments passed to plot
#' @return An object of class "plot.gm.prcomp" is a list with components
#'  that can be used in other plot functions, such as the type of plot, points, 
#'  a group factor, and other information depending on the plot parameters used.
#'  
#'  NOTE: To visualize shape variation across PC axes, use \code{\link{picknplot.shape}}.
#' @export
#' @author Antigoni Kaliontzopoulou
#' @keywords utilities
#' @keywords visualization
#' @seealso  \code{\link{plotRefToTarget}}


plot.gm.prcomp <- function(x, axis1 = 1, axis2 = 2, phylo = FALSE, 
                           phylo.par = list(edge.color = "black", edge.width = 1, edge.lty = 1,
                                            node.bg = "black", node.pch = 21, node.cex = 1), ...) {
  options(warn = -1)
  if(NCOL(x$x) == 1) stop("Only one PC.  No plotting capability with this function.\n", 
                          call. = FALSE)
  v <- x$d/sum(x$d)
  xlabel <- paste("PC ", axis1, ": ", round(v[axis1] * 100, 2), "%", sep = "")
  ylabel <- paste("PC ", axis2, ": ", round(v[axis2] * 100, 2), "%", sep = "")
  plot.args <- list(x = x$x[, axis1], y = x$x[, axis2], xlab = xlabel, ylab = ylabel, ...)
  pcdata <- as.matrix(x$x[, c(axis1, axis2)])
  if(!is.null(plot.args$axes)) axes <- plot.args$axes else axes <- TRUE
  if(!is.logical(axes)) axes <- as.logical(axes)
  plot.args$xlim <- 1.05*range(plot.args$x)
  plot.args$ylim <- 1.05*range(plot.args$y)
  if(is.null(plot.args$asp)) plot.args$asp <- 1
  
  if(phylo) {
    phy <- x$phy
    phy.pcdata <- rbind(x$x, x$anc.x)
    phy.pcdata <- as.matrix(phy.pcdata[, c(axis1, axis2)])
    plot.args$x <- pcdata[,1]
    plot.args$y <- pcdata[,2]

  }
  
  do.call(plot, plot.args)
  
  if(phylo) {
    for (i in 1:nrow(phy$edge)) {
      dt.xy <- xy.coords(phy.pcdata[phy$edge[i,], ])
      plot.xy(dt.xy, type="l", col = phylo.par$edge.color, 
              lwd = phylo.par$edge.width, lty = phylo.par$edge.lty)
    }
    plot.xy(xy.coords(phy.pcdata[1:length(phy$tip),]), type="p",...)
    plot.xy(xy.coords(phy.pcdata[(length(phy$tip)+1):nrow(phy.pcdata),]), type="p",
            pch = phylo.par$node.pch, cex = phylo.par$node.cex, bg = phylo.par$node.bg)
  }

  if(axes){
    abline(h = 0, lty=2, ...)
    abline(v = 0, lty=2, ...)
  }

  options(warn = 0)
  out <- list(PC.points = pcdata,   
              call = match.call())
  out$GM <- list()
  out$GM$A <- x$A
  
  out$plot.args <- plot.args
  if(phylo) {
    out$phylo <- list()
    out$phylo$phy <- phy
    out$phylo$phylo.par <- phylo.par
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
