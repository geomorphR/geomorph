#' Shape prediction from numeric predictors
#'
#' Function estimates one or more configurations based on one or more linear predictors, such as PC scores
#' allometric relationships, or any other least squares or partial least squares regression.  These configurations
#' can be used with \code{\link{plotRefToTarget}} to generate graphical representations of shape change, based on prediction criteria. 
#' 
#' @param A A 3D array (p x k x n) containing Procrustes shape variables either from GPA or fitted values from a previous
#' analytical procedure. 
#' @param x Linear (numeric) predictors.  Can be a vector or a matrix, or a list containing vectors or matrices.  Values must
#' be numeric.  If a factor is desired, one should use \code{\link[stats]{model.matrix}} to obtain a design matrix.  This will
#' impact how prediction criteria need to be provided (see below).
#' @param Intercept Logical value to indicate whether an intercept should be used in the linear equation for predictions.  Generally, 
#' this value will be FALSE for shape predictions made in ordination plots.  It should be TRUE in cases where the expected
#' shape at the point the predictor has a value of 0 is not the mean shape.
#' @param method A choice between least squares (LS) or partial least squares (PLS) regression for prediction.  The function defaults
#' to LS prediction.  PLS might be chosen in cases where correlation is preferred over linear regression.  If PLS is chosen, a two-block
#' PLS analysis using  \code{\link{two.b.pls}} should be performed first, as only the first singular vector for predictors will be used
#' for defining prediction criteria (see below).
#' @param ... Any number of prediction criteria.  Criteria should be presented as either a scalar (if one predictor is provided) or 
#' a vector (if more than one predictor or a prediction matrix is provided); e.g., pred1 = c(0.1, -0.5), pred2 = c(-0.2, -0.1) (which
#' would be the case if two predictors were provided).  It is essential that the number of elements in any prediction criterion matches
#' the number of predictors.  Caution should be used when providing a design matrix to ensure that correct dummy variables are used in
#' prediction criteria, and that either 1) an intercept is not included in the design and 2) is TRUE in the Intercept argument; or
#' or 1) an intercept is included in the design and 2) is FALSE in the Intercept argument; or 1) an intercept is not included in the design
#' and 2) is FALSE in the Intercept argument, if no intercept is desired.
#' @export
#' @author Michael Collyer
#' @return A list of predicted shapes matching the number of vectors of prediction criteria provides.  The predictions 
#' also have names matching those of the prediction criteria.
#' @examples
#' # Examples using Plethodon data
#' # NOT RUN
#' # data("plethodon")
#'
#' # Y.gpa <- gpagen(plethodon$land)    #GPA-alignment    
#' # plot(gm.prcomp(Y.gpa$coords))
#'
#' # preds <- shape.predictor(Y.gpa$coords, x= NULL, Intercept = FALSE, 
#' # pred1 = -0.1, pred2 = 0.1) # PC 1 extremes, sort of
#' # M <- mshape(Y.gpa$coords)
#' # plotRefToTarget(M, preds$pred1)
#' # plotRefToTarget(M, preds[[1]]) # same result
#' # plotRefToTarget(M, preds$pred2)
#'
#' # PCA <- gm.prcomp(Y.gpa$coords)
#' # PC <- PCA$x[,1]
#' # preds <- shape.predictor(Y.gpa$coords, x= PC, Intercept = FALSE, 
#' # pred1 = min(PC), pred2 = max(PC)) # PC 1 extremes, more technically
#' # plotRefToTarget(M, preds$pred1)
#' # plotRefToTarget(M, preds$pred2)
#'
#' # PC <- PCA$x[,1:2]
#'# user-picked spots can be anything, but it in this case, apparent groups
#' # preds <- shape.predictor(Y.gpa$coords, x= PC, Intercept = FALSE, 
#' #                          pred1 = c(0.045,-0.02), 
#' #                          pred2 = c(-0.025,0.06), 
#' #                          pred3 = c(-0.06,-0.04)) 
#' # plotRefToTarget(M, preds$pred1)
#' # plotRefToTarget(M, preds$pred2)
#' # plotRefToTarget(M, preds$pred3)
#'
#'# allometry example - straight-up allometry
#'
#' # preds <- shape.predictor(Y.gpa$coords, x= log(Y.gpa$Csize), 
#' #                          Intercept = TRUE, 
#' #                          predmin = min(log(Y.gpa$Csize)), 
#' #                          predmax = max(log(Y.gpa$Csize))) 
#'
#' # plotRefToTarget(M, preds$predmin, mag=3)
#' # plotRefToTarget(M, preds$predmax, mag=3)
#'
#' # allometry example - using RegScore or PredLine via procD.lm
#'
#' # gdf <- geomorph.data.frame(Y.gpa)
#' # plethAllometry <- procD.lm(coords ~ log(Csize), data=gdf)
#' # allom.plot <- plot(plethAllometry, 
#' # type = "regression", 
#' # predictor = log(gdf$Csize),
#' # reg.type ="RegScore") # make sure to have a predictor 
#'
#' # preds <- shape.predictor(plethAllometry$GM$fitted, 
#' #                          x= allom.plot$RegScore, Intercept = FALSE, 
#' #                          predmin = min(allom.plot$RegScore), 
#' #                          predmax = max(allom.plot$RegScore)) 
#' # plotRefToTarget(M, preds$predmin, mag=3)
#' # plotRefToTarget(M, preds$predmax, mag=3)
#'
#' # allom.plot <- plot(plethAllometry, 
#' # type = "regression", 
#' # predictor = log(gdf$Csize),
#' # reg.type ="PredLine")
#' # preds <- shape.predictor(plethAllometry$GM$fitted, 
#' #                          x= allom.plot$PredLine, Intercept = FALSE, 
#' #                          predmin = min(allom.plot$PredLine), 
#' #                          predmax = max(allom.plot$PredLine)) 
#' # plotRefToTarget(M, preds$predmin, mag=3)
#' # plotRefToTarget(M, preds$predmax, mag=3)
#'
#'# using factors via PCA
#'
#' # gdf <- geomorph.data.frame(Y.gpa, species = plethodon$species, 
#' #         site = plethodon$site)
#' # pleth <- procD.lm(coords ~ species*site, data=gdf)
#' # PCA <- prcomp(pleth$fitted)
#' # plot(PCA$x, asp=1, pch=19)
#'
#' # means <- unique(round(PCA$x,3))
#' # means # note: suggests 3 PCs useful enough
#'
#' # preds <- shape.predictor(arrayspecs(pleth$fitted, 12,2), x= PCA$x[,1:3],
#' #                          Intercept = FALSE,
#' #                          pred1 = means[1,1:3], 
#' #                          pred2 = means[2,1:3],
#' #                          pred3 = means[3,1:3], 
#' #                          pred4 = means[4,1:3])                   
#' # plotRefToTarget(M, preds$pred1, mag=2)
#' # plotRefToTarget(M, preds$pred2, mag=2)
#' # plotRefToTarget(M, preds$pred3, mag=2)
#' # plotRefToTarget(M, preds$pred4, mag=2)
#'
#'# Using a design matrix for factors
#'
#' # X <- pleth$X
#' # X # includes intercept; remove for better functioning 
#' # X <- X[,-1]
#' # symJord <- c(0,1,0) # design for P. Jordani in sympatry
#' # alloJord <- c(0,0,0) # design for P. Jordani in allopatry
#' # preds <- shape.predictor(arrayspecs(pleth$fitted, 12,2), x = X, 
#' #                          Intercept = TRUE, 
#' #                          symJord=symJord, alloJord=alloJord)
#' # plotRefToTarget(M, preds$symJord, mag=2)
#' # plotRefToTarget(M, preds$alloJord, mag=2)
#'
#'# PLS Example
#'
#' # data(plethShapeFood) 
#' # Y.gpa<-gpagen(plethShapeFood$land)    #GPA-alignment    
#'
#'# 2B-PLS between head shape and food use data
#' # PLS <-two.b.pls(A1 = plethShapeFood$food, A2 = Y.gpa$coords, iter=999) 
#' # summary(PLS)
#' # plot(PLS)
#'
#' # preds <- shape.predictor(Y.gpa$coords, plethShapeFood$food, 
#' #                          Intercept = FALSE,
#' #                          method = "PLS",
#' #                          pred1 = 2, pred2 = -4, pred3 = 2.5) 
#'                         # using PLS plot as a guide
#' # M <- mshape(Y.gpa$coords)
#' # plotRefToTarget(M, preds$pred1, mag=2)
#' # plotRefToTarget(M, preds$pred2, mag=2)
#' # plotRefToTarget(M, preds$pred3, mag=2)
#'
shape.predictor <- function(A, 
                            x = NULL, 
                            Intercept = FALSE, 
                            method = c("LS", "PLS"), ...){
  if(length(dim(A)) != 3) stop("Shape data must be in the form of a p x k x n array")
  dims <- dim(A); p <- dims[1]; k <- dims[2]; n <- dims[3]
  Y <- two.d.array(A)
  if(is.list(x)){
    if(any(sapply(x, is.factor)) == TRUE) stop("Predictors must be numeric.  Consider using model.matrix first.")
    if(any(sapply(x, is.character)) == TRUE) stop("Predictors must be numeric.  Consider using model.matrix first.")
    if(any(sapply(x, is.logical)) == TRUE) stop("Predictors must be numeric.  Consider using model.matrix first.")
    Ns <- sapply(x, NROW)
    if(length(unique(Ns)) > 1) stop("Predictors have different numbers of observations")
    N <- Ns[1]
    if(N != n) (stop("Predictors and Shape variables do not match in length"))
    X <- matrix(unlist(x),N)
  }
  if(!is.list(x)){
    if(is.matrix(x) || is.vector(x)) {
      if(is.factor(x)) stop("Predictors must be numeric.  Consider using model.matrix first.")
      if(is.character(x)) stop("Predictors must be numeric.  Consider using model.matrix first.")
      if(is.logical(x)) stop("Predictors must be numeric.  Consider using model.matrix first.")
      X <-x    
    }
  }
  if(is.null(x)) {
    pca <- prcomp(Y)
    X <- pca$x[,1]
  }
  X <- as.matrix(X)
  if(is.null(method)) method <- "LS" else
    method = match.arg(method, c("LS", "PLS"))
  if(method == "PLS"){
    plsRaw <- pls(X, Y, verbose=TRUE)
    XScores <- as.matrix(plsRaw$XScores); YScores <- as.matrix(plsRaw$YScores)
    X <- as.matrix(XScores[,1])
    Y <- predict(lm(Y~YScores+0))
  }
  colnames(X)<-paste("P",1:NCOL(X), sep="")
  form <- Y~X
  if(method == "PLS") Intercept <- FALSE
  if(!Intercept) form <- update(form, ~.+0)
  dat <- data.frame(Y=Y, X=X)
  fit <- lm(form, data=dat)
  B <- fit$coefficients
  dots <- list(...)
  if(is.null(dots)) stop("No prediction criteria given")
  np <- length(dots)
  if(length(unlist(dots)) != np*NCOL(X)) 
    stop("\nPrediction crtieria not consistent with the number of predictors.  \nCheck number of values input")
  dots <- lapply(dots, as.vector)
  if(Intercept){
    add.int <- function(x) c(1,x)
    dots <- lapply(dots, add.int)
  }
  pr.shape <- function(x) crossprod(x,B)
  preds <- lapply(dots, pr.shape)
  M <- mshape(A)
  if(is.null(rownames(M))) rownames(M) <- 1:nrow(M)
  if(is.null(colnames(M))) colnames(M) <- 1:ncol(M)
  M0 <- M; M0[M0!=0] <- 0
  if(!Intercept) reshape <- function(x) M + arrayspecs(x,p,k)[,,1] else
    reshape <- function(x) M0 + arrayspecs(x,p,k)[,,1] 
  preds <- lapply(preds, reshape)
  return(preds)
}
