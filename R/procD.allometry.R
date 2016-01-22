#' Procrustes ANOVA/regression, specifically for shape-size covariation (allometry)
#'
#' Function performs Procrustes ANOVA with permutation procedures to assess statistical hypotheses describing 
#'   patterns of shape covariation with size for a set of Procrustes-aligned coordinates.  Other factors or
#'   covariates can also be included in the analysis.  This function also provides results for plotting allometric
#'   curves.
#'
#' The function quantifies the relative amount of shape variation attributable to covariation with organism size (allometry)
#' plus other factors in a linear model, plus estimates the probability of this variation ("significance") for a null model, 
#' via distributions generated from resampling permutations. Data input is specified by formulae (e.g., 
#'   Y ~ X), where 'Y' specifies the response variables (shape data), and 'X' contains one or more independent 
#'   variables (discrete or continuous). The response matrix 'Y' can be either in the form of a two-dimensional data 
#'   matrix of dimension (n x [p x k]), or a 3D array (p x n x k).  It is assumed that  -if the data are based
#'   on landmark coordinates - the landmarks have previously been aligned using Generalized Procrustes Analysis (GPA) 
#'   [e.g., with \code{\link{gpagen}}]. 
#'   
#'   There are three formulae that need to be input (see Arguments). The first must contain variables for shape and size,
#'   e.g., Y ~ X, where Y (dependent variable) is shape and X (independent variable) is size.  The other two formulae
#'   are optional to indicate (1) groups for separate allometric curves and (2) additional model variables to consider in
#'   the ANOVA.  The groups input must be a single factor or multiple factors; e.g., ~ group, or ~ a*b.
#'   The resulting ANOVA uses sequential (Type I) sums of squares and cross-products with variables in this order:
#'   size, groups (if provided), size*groups (if warranted), other variables (if provided).  If a factor for groups is provided,
#'   ANOVA for a "homogeneity of slopes" test will also be performed.
#'   
#'   It is assumed that the order of the specimens in the shape matrix matches the order of values in the independent variables.  
#'   Linear model fits (using the  \code{\link{lm}} function) can also be input in place of formulae.  
#'   Arguments for \code{\link{lm}} can also be passed on via this function.  For further information about ANOVA in geomorph, resampling
#'   procedures used, and output, see \code{\link{procD.lm}} or \code{\link{advanced.procD.lm}}.
#'   If greater flexibility is required for variable order, \code{\link{advanced.procD.lm}} should be used.
#'   
#'   It is recommended that \code{\link{geomorph.data.frame}} is used to create and input a data frame.  This will reduce problems caused
#'   by conflicts between the global and function environments.  In the absence of a specified data frame, \code{\link{procD.allometry}} 
#'   will attempt to coerce input data into a data frame, but success is not guaranteed.
#'
#'   The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{procD.allometry}}.
#'   The generic function, \code{\link{plot}}, produces plots of allometric curves, using one of  three methods input (see below).
#'   If diagnostic plots on model residuals are desired, \code{\link{procD.lm}} should be used with the resulting model formula.  
#'   This, along with the data frame resulting from analysis with \code{\link{procD.allometry}} can be used directly in \code{\link{procD.lm}},
#'   which might be useful for extracting ANOVA components (as \code{\link{procD.allometry}} 
#'   is far more basic than \code{\link{procD.lm}}, in terms of output).
#'   
#' \subsection{Notes for geomorph 3.0 and making allometry plots}{ 
#' Former versions of geomorph had a "plotAllometry" function that performed ANOVA and produced
#' plots of allometry curves.  In geomorph 3.0, the \code{\link{plot}} function is used with 
#' \code{\link{procD.allometry}} objects to produce such plots.  The following arguments can be used in 
#' \code{\link{plot}} to achieve desired results.
#' \itemize{
#' \item{method = ("CAC, "RegScore, "PredLine").  Choose the desired plot method.}
#' \item{warpgrids: default = TRUE.  Logical value to indicate whether warpgrids should be plotted.} 
#' (Only workds with 3D array data)
#' \item{label: can be logical to label points (1:n) - e.g., label = TRUE - or a vector indicating
#' text to use as labels.}
#' \item{mesh: A mesh3d object to be warped to represent shape deformation of the minimum and maximum size 
#' if {warpgrids=TRUE} (see \code{\link{warpRefMesh}}).}
#' }
#' Use ?\code{\link{plot.procD.allometry}} to understand the arguments used.  The following are brief 
#' descriptions of the different plotting methods using \code{\link{plot}}, with references.
#'\itemize{
#' \item {If "method=CAC" (the default) the function calculates the 
#'   common allometric component of the shape data, which is an estimate of the average allometric trend 
#'   within groups (Mitteroecker et al. 2004). The function also calculates the residual shape component (RSC) for 
#'   the data.}
#'   \item {If "method=RegScore" the function calculates shape scores 
#'   from the regression of shape on size, and plots these versus size (Drake and Klingenberg 2008). 
#'   For a single group, these shape scores are mathematically identical to the CAC (Adams et al. 2013).}
#'   \item {If "method=PredLine" the function calculates predicted values from a regression of shape on size, and 
#'   plots the first principal component of the predicted values versus size as a stylized graphic of the 
#'   allometric trend (Adams and Nistri 2010). }
#'   }
#'   }
#'   
#' @param f1 A formula for the relationship of shape and size; e.g., Y ~ X.
#' @param f2 An optional right-hand formula for the inclusion of groups; e.g., ~ groups.
#' @param f3 A optional right-hand formula for the inclusion of additional variables; e.g., ~ a + b + c + ...
#' @param logsz A logical argument to indicate if the variable for size should be log-transformed.
#' @param iter Number of iterations for significance testing
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.  
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param alpha The significance level for the homegeneity of slopes test
#' @param RRPP A logical value indicating whether residual randomization should be used for significance testing
#' @param data A data frame for the function environment, see \code{\link{geomorph.data.frame}} 
#' @param ... Arguments passed on to procD.fit (typically associated with the lm function)
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @return An object of class "procD.allometry" is a list containing the following:
#' \item{HOS.test}{ANOVA for a homogeneity of slopes test (if groups are provided).}
#' \item{aov.table}{An analysis of variance table, based on inputs and the homogenetiy of slopes test.}
#' \item{alpha}{The significance level criterion for the homogeneity of slopes test.}
#' \item{perm.method}{A value indicating whether "RRPP" or randomization of "raw" vales was used.}
#' \item{permutations}{The number of random permutations used in the resampling procedure.}
#' \item{call}{The matched call.}
#' \item{formula}{The resulting formula, which can be used in follow-up analyses.  Irrespective of input, shape = Y
#' in the formula, and the variable used for size is called "size".}
#' \item{CAC}{The 
#'   common allometric component of the shape data, which is an estimate of the average allometric trend 
#'   within groups (Mitteroecker et al. 2004). The function also calculates the residual shape component (RSC) for 
#'   the data.}
#' \item{RSC}{The residual shape component (associated with CAC approach)}
#' \item{Reg.proj}{The projected regression scores on the regression of shape on size. 
#'   For a single group, these shape scores are mathematically identical to the CAC (Adams et al. 2013).}
#' \item{pred.val}{Principal component scores (first PC) of predicted values.}
#' \item{ref}{the reference configuration (if input coordinates are in a 3D array).}
#' \item{gps}{A vector of group names.}
#' \item{size}{A vector of size scores.}
#' \item{logsz}{A logical value to indicate if size values were log=transformed for analysis.}
#' \item{A}{Procrustes (aligned) residuals.}
#' \item{Ahat}{Predicted Procrustes residuals(if input coordinates are in a 3D array).}
#' \item{p}{landmark number}
#' \item{k}{landmark dimensions}
#' 
#' @references Adams, D.C., F.J. Rohlf, and D.E. Slice. 2013. A field comes of age: geometric morphometrics 
#'   in the 21st century. Hystrix. 24:7-14. 
#' @references Adams, D. C., and A. Nistri. 2010. Ontogenetic convergence and evolution of foot morphology 
#'   in European cave salamanders (Family: Plethodontidae). BMC Evol. Biol. 10:1-10.
#' @references Drake, A. G., and C. P. Klingenberg. 2008. The pace of morphological change: Historical 
#'   transformation of skull shape in St Bernard dogs. Proc. R. Soc. B. 275:71-76.
#' @references Mitteroecker, P., P. Gunz, M. Bernhard, K. Schaefer, and F. L. Bookstein. 2004. 
#'   Comparison of cranial ontogenetic trajectories among great apes and humans. J. Hum. Evol. 46:679-698.
#' @seealso \code{\link{procD.lm}} and \code{\link{advanced.procD.lm}} within geomorph;
#' \code{\link[stats]{lm}} for more on linear model fits
#' @examples
#' # Simple allometry
#' data(plethodon) 
#' Y.gpa <- gpagen(plethodon$land)    #GPA-alignment  
#' gdf <- geomorph.data.frame(Y.gpa, site = plethodon$site, 
#' species = plethodon$species) # geomorph data frame
#' plethAllometry <- procD.allometry(coords~Csize, f2 = NULL, f3=NULL, 
#' logsz = TRUE, data=gdf, iter=499)
#' summary(plethAllometry)
#' plot(plethAllometry, method = "PredLine")
#' plot(plethAllometry, method = "RegScore")
#' 
#' # Group Allometries
#' plethAllometry <- procD.allometry(coords~Csize, ~species*site, 
#' logsz = TRUE, data=gdf, iter=499, RRPP=TRUE)
#' summary(plethAllometry)
#' plot(plethAllometry, method = "PredLine")
#' 
#' # Using procD.lm to perform diagnostic residual plots
#' plethANOVA <- procD.lm(plethAllometry$formula, 
#' data = plethAllometry$data, iter = 499, RRPP=TRUE)
#' summary(plethANOVA) # Same ANOVA
#' plot(plethANOVA) # diagnostic plot instead of allometry plot
procD.allometry<- function(f1, f2 = NULL, f3 = NULL, logsz = TRUE,
                   iter = 999, seed=NULL, alpha = 0.05, RRPP = TRUE, data=NULL, ...){
  pfit <- procD.fit(f1, data=data, pca=FALSE)
  dat <- pfit$data
  Y <- pfit$Y
  if(!is.null(seed) && seed=="random") seed = sample(1:iter, 1)
  if((ncol(dat) - ncol(Y)) != 1) stop("Only a single covariate for size is permitted") 
  dat <- data.frame(Y=Y, size = dat[,ncol(dat)])
  size <- dat$size
  if(logsz) {
    dat <- model.frame(Y ~ log(size), data = dat) 
    form1 <- Y ~ log(size)
  } else {
    dat <- model.frame(Y ~ size, data = dat)
    form1 <- Y ~ size
  }
  
  if(!is.null(f2) || !is.null(f3)){
    if(!is.null(data)) {
      data.types <- lapply(data, class)
      keep = sapply(data.types, function(x) x != "array" & x != "phylo")
      dat2 <- as.data.frame(data[keep])
      } else dat2 <- NULL
      
    if(!is.null(f2)) {
      if(length(f2) > 2) f2 <- f2[-2]
      dat.g <- model.frame(f2, data=dat2) 
      dat <- data.frame(dat, dat.g)
      g.Terms <- terms(dat.g)
      if(any(attr(g.Terms, "dataClasses") == "numeric")) stop("groups formula (f2) must contain only factors")
      if(ncol(dat.g) > 1) gps <- factor(apply(dat.g, 1,function(x) paste(x, collapse=":"))) else 
        gps <- as.factor(unlist(dat.g))
    }  else {
      dat.g <- NULL
      g.Terms <- NULL
      gps <- NULL
    }
  
    if(!is.null(f3)) {
      if(length(f3) > 2) f3 <- f3[-2]
      dat.o <- model.frame(f3, data=dat2) 
      dat <- data.frame(dat, dat.o)
      o.Terms <- terms(dat.o)
    }  else {
      dat.o <- NULL
      o.Terms <- NULL
    }
  }

    if(!is.null(f2) & !is.null(f3)) {
      if(!logsz){
        form4 <- update(f3, ~. + size + gps)
        form5 <- update(f3, ~. + size * gps)
      }
      if(logsz){
        form4 <- update(f3, ~. + log(size) + gps)
        form5 <- update(f3, ~. + log(size) * gps)
      }
    }
      if(!is.null(f2) & is.null(f3)) {
        form4 <- form1
        form5 <- update(form1, ~.  * gps)
      }
  
  if(!is.null(f2) & !is.null(f3)) {
    formfull <-as.formula(c("~",paste(unique(
      c(c("size", attr(g.Terms, "term.labels"), paste("size", attr(g.Terms, "term.labels"), sep=":")),
        c("size", attr(o.Terms, "term.labels"), paste("size", attr(o.Terms, "term.labels"), sep=":")))),
      collapse="+")))
    form.type <- "go"
  } else if(!is.null(f2) & is.null(f3)) {
    if(!logsz) formfull <-as.formula(c("~",paste(unique(
      c(c("size", attr(g.Terms, "term.labels"), paste("size", attr(g.Terms, "term.labels"), sep=":")))),
      collapse="+"))) else
        formfull <-as.formula(c("~",paste(unique(
          c(c("log(size)", attr(g.Terms, "term.labels"), paste("log(size)", attr(g.Terms, "term.labels"), sep=":")))),
          collapse="+")))
    form.type <- "g"
  } else if(is.null(f2) & !is.null(f3)) {
    if(!logsz) formfull <-as.formula(c("~",paste(unique(
      c(c("size", attr(o.Terms, "term.labels"), paste("size", attr(o.Terms, "term.labels"), sep=":")))),
      collapse="+"))) else
        formfull <-as.formula(c("~",paste(unique(
          c(c("log(size)", attr(o.Terms, "term.labels"), paste("log(size)", attr(o.Terms, "term.labels"), sep=":")))),
          collapse="+")))
    form.type <- "o"
  } else {
    formfull <- form1
    form.type <- NULL}
  
  if(!is.null(f2)){
    form4 <- update(form4, Y ~.)
    form5 <- update(form5, Y ~.)
    datHOS <- data.frame(dat, size=size, gps=gps)
    HOS <- advanced.procD.lm(form4, form5, data=datHOS, iter=iter, seed=seed)$anova.table
    rownames(HOS) = c("Common Allometry", "Group Allometries")
    hos.pval <- HOS[2,7]
    if(hos.pval > alpha){
      if(form.type == "go") {
        if(!logsz) rhs.formfull <- paste(c("size", attr(g.Terms, "term.labels"), 
                                           attr(o.Terms, "term.labels")), collapse="+") else
                   rhs.formfull <- paste(c("log(size)", attr(g.Terms, "term.labels"), 
                                           attr(o.Terms, "term.labels")), collapse="+")  
        formfull <- as.formula(c("Y ~", rhs.formfull))
      }
      if(form.type == "g") {
        if(!logsz) rhs.formfull <- paste(c("size", attr(g.Terms, "term.labels")),  collapse="+") else
          rhs.formfull <- paste(c("log(size)", attr(g.Terms, "term.labels")),  collapse="+")
        formfull <- as.formula(c("Y ~", rhs.formfull))
      }
    } 
  } else HOS <- NULL
  
  formfull <- update(formfull, Y~.)
  fitf <- procD.fit(formfull, data=dat, pca=FALSE)
  anovafull <- procD.lm(formfull, data=dat, iter=iter, seed=seed, RRPP=RRPP)$aov.table
  if(RRPP) perm.method = "RRPP" else perm.method = "raw"
  
  # Plot set-up
  yhat <- fitf$wFitted[[length(fitf$wFitted)]]
  B <- fitf$wCoefficients[[length(fitf$wCoefficients)]]
  y.cent <- fitf$wResiduals[[length(fitf$wResiduals)]]
  if(logsz) sz <- log(size) else sz = size
  a<-(t(y.cent)%*%sz)%*%(1/(t(sz)%*%sz)); a<-a%*%(1/sqrt(t(a)%*%a))
  CAC<-y.cent%*%a  
  resid<-y.cent%*%(diag(dim(y.cent)[2])-a%*%t(a))
  RSC<-prcomp(resid)$x
  Reg.proj<-Y%*%B[2,]%*%sqrt(solve(t(B[2,])%*%B[2,])) 
  pred.val<-prcomp(yhat)$x[,1] 
  if(!is.null(data)) lm.dim <- dim(data[[match(as.character(f1[[2]]), names(data))]]) else {
    Z <- eval(f1[[2]], parent.frame())
    lm.dim <- dim(Z)
  }
  Ahat <- arrayspecs(yhat, lm.dim[[1]], lm.dim[[2]])
  A <- arrayspecs(Y, lm.dim[[1]], lm.dim[[2]])
  ref<-mshape(A)
  if(is.null(f2)) gps <- NULL
  out <- list(HOS.test = HOS, aov.table =anovafull, call = match.call(),
              alpha = alpha, perm.method = perm.method, permutations=iter+1,
              formula = formfull, data=dat,
              CAC = CAC, RSC=RSC, Reg.proj = Reg.proj,
              pred.val=pred.val,
              ref=ref, gps=gps, size=size, logsz=logsz, 
              A=A, Ahat=Ahat, p=lm.dim[[1]], k= lm.dim[[2]])
  class(out) <- "procD.allometry"
  out
}