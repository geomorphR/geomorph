#' Plotting to assist visualization of shape-size covariation (allometry)
#'
#' Function performs plotting for a procD.lm fit and a vector of size measures.
#' 
#' Prior to geomorph 3.0.0, the function, plotAllometry, was used to perform linear regression
#' of shape variables and size, and produce plots to visualize shape allometries.  This function was deprecated
#' when procD.allometry was launched with geomorph 3.0.0, which performed homogeneity of slopes tests to determine
#' if a common allometry or unique group allometries were more appropriate as a model.  The S3 generic, plot.procD.allometry
#' provided the same plotting as plotAllometry before it.  In geomorph 3.1.0, procD.allometry has been deprecated in favor of using
#' \code{\link{procD.lm}} and \code{\link{pairwise}} for analyses, which can include additional variables, 
#' thus eliminating plot.procD.allometry.  This function coalesces a few plotting options found in other functions, 
#' as a wrapper, for the purpose of retaining the plot.procD.allometry options in one place.  
#' 
#'
#' There are fundamentally two different kinds of allometry plots: those based on linear models and those that do not have a linear
#' model basis (more detail below).  The common allometric component (CAC) and size-shape PCA (Mitteroecker et al. 2004) are plotting 
#' strategies that do not have results
#' that vary with linear model parameters.  By contrast, prediction lines (PredLine, Adams and Nistri 2010) and regression scores 
#' (RegScore, Drake and Klingenberg 2008) are based on fitted values and regression coefficients, respectively, to visualize allometric patterns.  
#' The plotAllometry function will extract necessary components from a \code{\link{procD.lm}} fit to calculate these various statistics 
#' (although the variables used in the \code{\link{procD.lm}} fit are inconsequential for CAC and size-shape PCA; only the shape variables are used).
#' 
#' There are multiple ways to visualize allometry.  One way is to simply append a size variable to shape variables and perform a principal component analysis
#' (PCA).  In the event that size and shape strongly covary, the first PC scores might reflect this (Mitteroecker et al. 2004).  Alternatively, the major
#' axis of covariation between size and shape can be found by a singular value decomposition of their cross-products, a process known as two-block partial 
#' least squares (PLS; Rohlf and Corti 2000).  This major axis of variation is often referred to as the common allometric component 
#' (CAC; Mitteroecker et al. 2004).  Neither of these methods is associated with a model of allometric shape change, especially as such change might vary
#' for different groups.  As such, these methods have limited appeal for comparing group allometries (although color-coding groups in plots might reveal
#' different trends in the plot scatter).
#' 
#' By contrast, describing a linear model (with \code{\link{procD.lm}}) that has an explicit definition of how shape allometries vary by group can be 
#' more informative.  The following are the three most general models:
#' 
#' simple allometry: shape ~ size
#' 
#' common allometry, different means: shape ~ size + groups
#' 
#' unique allometries: shape ~ size * groups
#' 
#' However, other covariates can be added to these models.  One could define these models with \code{\link{procD.lm}}
#' and use \code{\link{anova.lm.rrpp}} to explicitly test which model
#' is most appropriate.  The function, \code{\link{pairwise}} can also be used to test pairwise differences among least-squares means or slopes.
#' To visualize different allometric patterns, wither prediction lines (PredLine; Adams and Nistri 2010) or regression scores 
#' (RegScore; Drake and Klingenberg 2008) can be used.  The former plots first PCs of fitted values against size; the latter calculates a regression score
#' as a projection of data on normalized vector that expresses the covariation between shape and the regression coefficients for size, conditioned
#' on other model effects.  For a simple allometry model, CAC and RegScore are the same (Adams et al. 2013) but RegScore, like PredLine but unlike CAC,
#'  generalizes to complex models.
#' Either PredLine or RegScore can help elucidate divergence in allometry vectors among groups.
#' 
#' If the variable for size is used in the \code{\link{procD.lm}} fit, the plot options will resemble past allometry plots found in
#' geomorph.  However, with this updated function philosophy, the model fit does not have to necessarily contain size.  This might be useful if 
#' one wishes to visualize whether shape, size, and some other variable covary in some way (by first performing a \code{\link{procD.lm}} fit 
#' between shape and another covariate, then performing plotAllometry with that fit and size).  For example, one can entertain the question,
#' "Are species differences in shape merely a manifestation of shape allometry, when species differ in size?"  By fitting a model, shape ~ species,
#'  then using plotAllometry for the model fit (with either PredLine or RegScore), the plot will help reveal if allometry and species effects are confounded.
#'
#' 
#' The following are brief descriptions of the different plotting methods, with references.
#' 
#'\itemize{
#'   \item {If "method = PredLine" (the default) the function calculates fitted values from a \code{\link{procD.lm}} fit, and 
#'   plots the first principal component of the "predicted" values versus size as a stylized graphic of the 
#'   allometric trend (Adams and Nistri 2010). This method is based on linear models and 
#'   can allow for other model variable to be incorporated.}
#'   \item {If "method = RegScore" the function calculates standardized shape scores 
#'   from the regression of shape on size, and plots these versus size (Drake and Klingenberg 2008). 
#'   For a single allometry, these shape scores are mathematically identical to the CAC (Adams et al. 2013).  
#'   This method is based on linear models and can allow for other model variable to be incorporated.}
#'   \item {If "method = size.shape" the function perform principal components analysis on a data space containing both shape 
#'   and size (sensu Mitteroecker et al. 2004).  This method is not based on linear models and results will not be changed by 
#'   changing the allometry model.}
#'   \item {If "method = CAC"  the function calculates the 
#'   common allometric component of the shape data, which is an estimate of the average allometric trend 
#'   for group-mean centered data (Mitteroecker et al. 2004). The function also calculates the residual shape component (RSC) for 
#'   the data.  This method is not based on linear models and results will not be changed by 
#'   changing the allometry model.}
#'   }
#'   
#' The function returns values that can be used with \code{\link{picknplot.shape}} or  
#' \code{\link{shape.predictor}} and \code{\link{plotRefToTarget}} to visualize shape changes in the plot.
#' 
#' @param fit A procD.lm fit.
#' @param size A vector of the same length as the number of observations in the fit.
#' @param logsz A logical value to indicate whether to first find the logarithm of size.
#' @param method The method of allometric visualization, which includes 
#' CAC, PredLine, RegScore, and size.shape (PCA)
#' @param ... Other arguments passed on to plot.default
#' @keywords utilities
#' @export
#' @author Michael Collyer
#' @return An object of class plotAllometry returns some combination of 
#' CAC values, the residual shape component (RSC, associated with CAC approach),
#' PredLine values, RegScore values, PC points for the size-shape PCA, and PCA statistics,
#' depending on arguments used.  The size variable and GM statistics from the original model fit
#' are also returned.
#' .
#' @references Adams, D. C., and A. Nistri. 2010. Ontogenetic convergence and evolution of foot morphology 
#'   in European cave salamanders (Family: Plethodontidae). BMC Evol. Biol. 10:1-10.
#' @references Adams, D.C., F.J. Rohlf, and D.E. Slice. 2013. A field comes of age: geometric morphometrics 
#'   in the 21st century. Hystrix. 24:7-14.
#' @references Drake, A. G., and C. P. Klingenberg. 2008. The pace of morphological change: Historical 
#'   transformation of skull shape in St Bernard dogs. Proc. R. Soc. B. 275:71-76.
#' @references Mitteroecker, P., P. Gunz, M. Bernhard, K. Schaefer, and F. L. Bookstein. 2004. 
#'   Comparison of cranial ontogenetic trajectories among great apes and humans. J. Hum. Evol. 46:679-698.
#' @references  Rohlf, F.J., and M. Corti. 2000. The use of partial least-squares to study covariation in shape. 
#' Systematic Biology 49: 740-753.
#' 
#' @examples 
#' 
#' # Simple allometry
#' data(plethodon) 
#' Y.gpa <- gpagen(plethodon$land, print.progress = FALSE)    #GPA-alignment  
#' 
#' gdf <- geomorph.data.frame(Y.gpa, site = plethodon$site, 
#' species = plethodon$species) 
#' fit <- procD.lm(coords ~ log(Csize), data=gdf, iter=0, 
#' print.progress = FALSE)
#' 
#' # Predline
#' plotAllometry(fit, size = gdf$Csize, logsz = TRUE, 
#' method = "PredLine", pch = 19)
#' 
#' # same as
#' logSize <- log(gdf$Csize)
#' plot(fit, type = "regression", reg.type = "PredLine", 
#' predictor = logSize, pch = 19)
#' 
#' # RegScore
#' plotAllometry(fit, size = gdf$Csize, logsz = TRUE, 
#' method = "RegScore", pch = 19)
#' 
#' # same as
#' plot(fit, type = "regression", reg.type = "RegScore", 
#' predictor = logSize, pch = 19)
#' 
#' # CAC
#' plotAllometry(fit, size = gdf$Csize, logsz = TRUE, 
#' method = "CAC", pch = 19)
#' 
#' # same (first plot) as
#' PLS <- two.b.pls(log(gdf$Csize), gdf$coords, print.progress = FALSE)
#' plot(PLS)
#' 
#' # Group Allometries
#' fit <- procD.lm(coords ~ Csize * species * site, data=gdf, iter=0, 
#' print.progress = FALSE)
#' 
#' # CAC (should not change from last time; model change has no effect)
#' plotAllometry(fit, size = gdf$Csize, logsz = TRUE, method = "CAC", 
#' pch = 19)
#' 
#' # Predline
#' plotAllometry(fit, size = gdf$Csize, logsz = TRUE, method = "PredLine", 
#' pch = 19, col = as.numeric(interaction(gdf$species, gdf$site)))
#' 
#' # RegScore
#' plotAllometry(fit, size = gdf$Csize, logsz = TRUE, method = "RegScore", 
#' pch = 19, col = as.numeric(interaction(gdf$species, gdf$site)))
#' 
#' # Size-Shape PCA
#' 
#' pc.plot <- plotAllometry(fit, size = gdf$Csize, logsz = TRUE, 
#' method = "size.shape", 
#' pch = 19, col = as.numeric(interaction(gdf$species, gdf$site)))
#' summary(pc.plot$size.shape.PCA)
#' 
#' # Are species' shape differences just a manifestation of shape allometry?
#' 
#' fit3 <- procD.lm(coords ~ species, data = gdf, iter = 0, 
#' print.progress = FALSE)
#' plotAllometry(fit3, size = gdf$Csize, logsz = TRUE, method = "RegScore", 
#' pch = 19, col = as.numeric(gdf$species))
#' 
#' # No evidence this is the case
#' 
plotAllometry <- function(fit, size, logsz = TRUE, 
        method = c("PredLine", "RegScore", "size.shape", "CAC"), ...){
  method <- match.arg(method)
  n <- length(size)
  if(n != fit$LM$n) 
    stop("Different number of observations between model fit and size.\n", call. = FALSE)
  if(!is.numeric(size))
    stop("The argument, size, must be a numeric vector.\n", call. = FALSE)
  if(!is.vector(size))
    stop("The argument, size, must be a numeric vector.\n", call. = FALSE)
  
  dat <- fit$LM$data
  Y <- fit$LM$Y
  form <- if(logsz) as.formula(Y ~ log(size)) else as.formula(Y ~ size)
  dat$size <- size
  GM <- fit$GM
  xc <- if(logsz) log(size) else size
  
  if(method == "PredLine") {
    
    if(logsz) out <- plot(fit, type = "regression", 
                  reg.type = "PredLine", 
                  predictor = log(size), ...) else
                    out <- plot(fit, type = "regression", 
                                reg.type = "PredLine", 
                                predictor = size, ...)          
    
  } else if(method == "RegScore") {
    
    if(logsz) out <- plot(fit, type = "regression", 
                          reg.type = "RegScore", 
                          predictor = log(size), ...) else
                            out <- plot(fit, type = "regression", 
                                        reg.type = "RegScore", 
                                        predictor = size, ...) 
    
  } else 
      
      {
        if(fit$LM$gls){
          if(!is.null(fit$LM$weights))
            fit <- lm.rrpp(form, data = dat, weights = fit$LM$weights, 
                            iter = 0, print.progress = FALSE)
          if(!is.null(fit$LM$Cov))
            fit <- lm.rrpp(form, data = dat, Cov = fit$LM$Cov, 
                            iter = 0, print.progress = FALSE)
        } else {
          fit <- lm.rrpp(form, data = dat, iter = 0, print.progress = FALSE)
        }
      
        f <- if(fit$LM$gls) fit$LM$gls.fitted else fit$LM$fitted
        b <- as.matrix(lm(f ~ xc)$coefficients)[2,]
        y <- fit$LM$Y
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
        
        if(method == "CAC") {
          
          par(mfcol = c(1,2))
          plot.args <- list(x = xc, y = CAC,
                            xlab = if(logsz) "log(Size)" else "Size",
                            ylab = "CAC", ...)
          plot2.args <- list(x = CAC, y = RSC[, 1],
                             xlab = "CAC", ylab = "RSC1", ...)
          do.call(plot, plot.args)
          do.call(plot, plot2.args)
          par(mfcol = c(1,1))
          out <- list(CAC = CAC, RSC = RSC, plot.args = plot.args,
                      all.plot.args = list(CAC=plot.args, RSC = plot2.args))
        }
        
        if(method == "size.shape") {
          
          v <- round(PCA$sdev^2/sum(PCA$sdev^2) * 100, 2)
          plot.args <- list(x = PC.points[,1], y = PC.points[,2],
                            xlab = paste("PC 1: ", v[1], "%", sep = ""), 
                            ylab = paste("PC 2: ", v[2], "%", sep = ""),
                            ...)
          do.call(plot, plot.args)
          title("Size-Shape PC plot")
          
          out <- list(size.shape.PCA = PCA,
                      PC.points = PC.points, plot.args = plot.args)
        }
      }
  
  out$size.var <- xc
  out$logsz <- logsz
  out$GM <- GM
  
  class(out) <- "plotAllometry"
  invisible(out)
  
}