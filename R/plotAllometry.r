#' Plotting to assist visualization of shape-size covariation (allometry)
#'
#' Function performs plotting for a procD.lm fit and a vector of size measures.
#' 
#' Prior to geomorph 3.0.0, the function, plotAllometry, was used to perform linear regression
#' of shape variables and size, and produce plots to visualize shape allometries.  This function was deprecated
#' when procD.allometry was launched with geomorph 3.0.0, which performed homogeneity of slopes tests to determine
#' if a common allometry or unique group allometries were more appropriate as a model.  The S3 generic, plot.procD.allometry
#' provided the same plotting as plotAllometry before it.  In geomorph 3.0.8, procD.allometry has been deprectaed in favor of using
#' \code{\link{procD.lm}} and \code{\link{pairwise}} for analyses, which can include additional variables, 
#' thus eliminating plot.procD.allometry.
#' 
#' Typical allometric plots can still be generated using procD.lm (for PredLine and RegScore plots) or two.b.pls (for CAC).  This
#' function is merely a wrapper for condensing such different allometry plot options, following shape analysis with \code{\link{procD.lm}}.  
#' 
#' There are fundamentally two different kinds of allometry plots: those based on linear models and those that do not have a linear
#' model basis.  The common allometric component (CAC) and size-shape PCA (Mitteroecker et al. 2004) are plotting strategies that do not have results
#' that vary with linear model parameters.  By contrast, predicition lines (PredLine, Adams and Nistri 2010) and regression scores 
#' (RegScore, Drake and Klingenberg 2008) use fitted values and regression coefficients, respectively, to visualize allometric patterns.  The plotAllometry
#' function will extract necessary components from a \code{\link{procD.lm}} fit to calculate these various statistics (although the variables
#' used in the \code{\link{procD.lm}} fit are inconsequntial for CAC and size-shape PCA; only the shape variables are used).
#' 
#' If the variable for size is used in the \code{\link{procD.lm}} fit, the plot options will resemble past allometry plots found in
#' geomorph.  However, with this updated function philosophy, the model fit does not have to necessarily contain size.  This might be useful if 
#' one wishes to visualize whether shape, size, and some other variable covary in some way (by first performing a \code{\link{procD.lm}} fit 
#' between shape and another covariate, then performing plotAllometry with that fit and size).
#' 
#' The following are brief descriptions of the different plotting methods, with references.
#' 
#'\itemize{
#' \item {If "method=CAC" (the default) the function calculates the 
#'   common allometric component of the shape data, which is an estimate of the average allometric trend 
#'   for group-mean centered data (Mitteroecker et al. 2004). The function also calculates the residual shape component (RSC) for 
#'   the data.}
#'   \item {If "method = RegScore" the function calculates shape scores 
#'   from the regression of shape on size, and plots these versus size (Drake and Klingenberg 2008). 
#'   For a single allometry, these shape scores are mathematically identical to the CAC (Adams et al. 2013).  This method can
#'   allow for other model variable to be incorporated.}
#'   \item {If "method = PredLine" the function calculates predicted values from a regression of shape on size, and 
#'   plots the first principal component of the predicted values versus size as a stylized graphic of the 
#'   allometric trend (Adams and Nistri 2010). This method can
#'   allow for other model variable to be incorporated.}
#'   \item {If "method = size.shape" the function perform principal components analysis on a data space containing both shape 
#'   and size (sensu Mitteroecker et al. 2004).}
#'   }
#'   
#' The function returns values that can be used with \code{\link{picknplot.shape}} to visualize shape changes in the plot.
#' 
#' @param fit A procD.lm fit.
#' @param size A vector of the same length as the numner of observations in the fit.
#' @param logsz A logical value to indicate whether to first find the logarithm of size.
#' @param method The method of allometric visualization; choice among CAC, PredLine, RegScore, and size.shape (PCA)
#' @param ... Other arguments passed on to plot.default
#' @keywords utilities
#' @export
#' @author Michael Collyer
#' @return An object of class plotAllometry returns CAC values, the residual shape component (RSC, associated with CAC approach),
#' PredLine values, RegScore values, the size variable,
#' PC points for the size-shape PCA, and PCA statistics.
#' @references Adams, D. C., and A. Nistri. 2010. Ontogenetic convergence and evolution of foot morphology 
#'   in European cave salamanders (Family: Plethodontidae). BMC Evol. Biol. 10:1-10.
#' @references Adams, D.C., F.J. Rohlf, and D.E. Slice. 2013. A field comes of age: geometric morphometrics 
#'   in the 21st century. Hystrix. 24:7-14.
#' @references Drake, A. G., and C. P. Klingenberg. 2008. The pace of morphological change: Historical 
#'   transformation of skull shape in St Bernard dogs. Proc. R. Soc. B. 275:71-76.
#' @references Mitteroecker, P., P. Gunz, M. Bernhard, K. Schaefer, and F. L. Bookstein. 2004. 
#'   Comparison of cranial ontogenetic trajectories among great apes and humans. J. Hum. Evol. 46:679-698.
#' 
#' @examples 
#' 
#' # Simple allometry
#' data(plethodon) 
#' Y.gpa <- gpagen(plethodon$land, print.progress = FALSE)    #GPA-alignment  
#' 
#' gdf <- geomorph.data.frame(Y.gpa, site = plethodon$site, 
#' species = plethodon$species) 
#' fit <- procD.lm(coords ~ log(Csize), data=gdf, iter=0, print.progress = FALSE)
#' 
#' # Predline
#' plotAllometry(fit, size = gdf$Csize, logsz = TRUE, method = "PredLine", pch = 19)
#' 
#' # same as
#' logSize <- log(gdf$Csize)
#' plot(fit, type = "regression", reg.type = "PredLine", predictor = logSize, pch = 19)
#' 
#' # RegScore
#' plotAllometry(fit, size = gdf$Csize, logsz = TRUE, method = "RegScore", pch = 19)
#' 
#' # same as
#' plot(fit, type = "regression", reg.type = "RegScore", predictor = logSize, pch = 19)
#' 
#' # CAC
#' plotAllometry(fit, size = gdf$Csize, logsz = TRUE, method = "CAC", pch = 19)
#' 
#' # same (first plot) as
#' PLS <- two.b.pls(log(gdf$Csize), gdf$coords, print.progress = FALSE)
#' plot(PLS, warpgrids = FALSE)
#' 
#' # Group Allometries
#' fit2 <- procD.lm(coords ~ Csize * species * site, data=gdf, iter=0, print.progress = FALSE)
#' 
#' # CAC (should not change from last time; model change has no effect)
#' plotAllometry(fit2, size = gdf$Csize, logsz = TRUE, method = "CAC", pch = 19)
#' 
#' # Predline
#' plotAllometry(fit2, size = gdf$Csize, logsz = TRUE, method = "PredLine", 
#' pch = 19, col = as.numeric(interaction(gdf$species, gdf$site)))
#' 
#' # RegScore
#' plotAllometry(fit2, size = gdf$Csize, logsz = TRUE, method = "RegScore", 
#' pch = 19, col = as.numeric(interaction(gdf$species, gdf$site)))
#' 
#' # Size-Shape PCA
#' 
#' pc.plot <- plotAllometry(fit2, size = gdf$Csize, logsz = TRUE, method = "size.shape", 
#' pch = 19, col = as.numeric(interaction(gdf$species, gdf$site)))
#' summary(pc.plot$size.shape.PCA)
#' 
#' # Are species' shape differences just a manifestation of shape allometry?
#' 
#' fit3 <- procD.lm(coords ~ species, data=gdf, iter=0, print.progress = FALSE)
#' plotAllometry(fit3, size = gdf$Csize, logsz = TRUE, method = "RegScore", 
#' pch = 19, col = as.numeric(gdf$species))
#' 
#' # No evidence this is the case
#' 
plotAllometry <- function(fit, size, logsz = TRUE, 
        method = c("CAC", "PredLine", "RegScore", "size.shape"), ...) {
  type <- match.arg(method)
  f <- if(fit$LM$gls) fit$LM$gls.fitted else fit$LM$wFitted
  y <- fit$LM$Y
  X <- fit$LM$X
  n <- NROW(X)
  if(logsz) xc <- log(size) else xc <- size
  if(length(xc) != n) stop("The number of size observations does not match the number of shape observations.\n", call. = FALSE)
  b <- as.matrix(lm.fit(cbind(xc, X), y)$coefficients)[1,]
  PL <- prcomp(f)$x[,1]
  a <- crossprod(center(y), xc)/sum(xc^2)
  a <- a/sqrt(sum(a^2))
  r <- center(y)
  CAC <- r%*%a  
  p <- fit$LM$p
  resid <- r%*%(diag(p) - matrix(crossprod(a),p,p))
  RSC <- prcomp(resid)$x
  Reg.proj <- r %*% b %*% sqrt(1/crossprod(b))
  PCA <- prcomp(cbind(y, xc))
  PC.points <-PCA$x
  
  if(type == "CAC") {
    par(mfcol = c(1,2))
    plot(xc, CAC, 
         xlab = if(logsz) "log(Size)" else "Size",
         ...)
    plot(CAC, RSC[,1],
         ylab = "RSC1", ...)
    par(mfcol = c(1,1))
    
  } else if(type == "RegScore") {
    plot(xc, Reg.proj,
         xlab = if(logsz) "log(Size)" else "Size",
         ylab = "Regression Score",
         ...)
  } else if(type == "size.shape") {
    v <- round(PCA$sdev^2/sum(PCA$sdev^2) * 100, 2)
    plot(PC.points, 
         main = "Size-Shape PC plot",
         xlab = paste("PC 1: ", v[1], "%", sep = ""),
         ylab = paste("PC 2: ", v[2], "%", sep = ""),
         ...)
    } else {
    plot(xc, PL,
         xlab = if(logsz) "log(Size)" else "Size",
         ylab = "PC 1 for fitted values",
         ...)
  }
  
  out <- list(CAC = CAC, PredLine = PL, RegScore = Reg.proj, PC.points = PC.points,
              size.var = xc, size.shape.PCA = PCA, logsz = logsz)
  class(out) <- "plotAllometry"
  invisible(out)
  
}