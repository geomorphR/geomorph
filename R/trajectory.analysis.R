#'  Quantify and compare shape change trajectories
#'
#'  Function estimates attributes of shape change trajectories or motion trajectories for a set of 
#'  Procrustes-aligned specimens and compares them statistically
#'
#'  The function quantifies phenotypic shape change trajectories from a set of specimens, and assesses variation 
#'  in these parameters via permutation. A shape change trajectory is defined by a sequence 
#'  of shapes in tangent space. These trajectories can be quantified for various attributes (their size, orientation, 
#'  and shape), and comparisons of these attribute enable the statistical comparison of shape change 
#'  trajectories (see Collyer and Adams 2013; Collyer and Adams 2007; Adams and Collyer 2007; Adams and Collyer 2009). 
#'
#'  Data input is specified by a two formulae (e.g., Y ~ X), where 'Y' specifies the response variables (trajectory data), 
#'  and 'X' contains one or more independent variables (discrete or continuous). The response matrix 'Y' can be either in the form of a two-dimensional data 
#'   matrix of dimension (n x [p x k]), or a 3D array (p x n x k).. The function
#'  \code{\link{two.d.array}} can be used to obtain a two-dimensional data matrix from a 3D array of landmark
#'  coordinates. It is assumed that the order of the specimens 'Y' matches the order of specimens in 'X'. 
#'  It is  also assumed that the landmarks have previously been aligned using Generalized Procrustes 
#'  Analysis (GPA) [e.g., with \code{\link{gpagen}}]. Linear model 
#'  fits (using the  \code{\link{lm}} function) can also be input in place of a formula.  Arguments for  \code{\link{lm}} 
#'  can also be passed on via this function.  The first formula, f1, must contain the independent variable on the left-hand
#'  side of the formula (e.g., Y ~) and either a single factor or a two factor interaction on the right-hand side.  If
#'  a single factor is provided, e.g., Y ~ A, it is assumed that groups to be described are the levels of factor A
#'  and that the data in Y comprise trajectories.  In this case, the traj.pts = NULL argument must be changed to a numeric value 
#'  to define the number of points in the trajectory.  It is also assumed that the data are structured as variables within points.  
#'  For example, y11 y21 y31 y12 y22 y32 y13 y23 y33 y14 y24 y34 would be columns of a matrix, Y, describing a 4-point trajectory
#'  in a data space defined by three variables.  This is the proper arrangment; the folowwing is an improper arrangement: 
#'  y11 y12 y13 y14 y21 y22 y23 y24 y31 y32 y33 y34, as it group points within variables.  This approach is typical when comparing
#'  motion paths (see Adams and Cerney 2007).
#'  
#'  If f1 is a two-factor factorial model, e.g., Y ~ A*B, it is assumed that the first factor defines groups, th second factor
#'  defines trajectory points, and that trajectories are to be estimated from the linear model.  In this case, the preceding example
#'  would have a Y matrix comprised only of y1, y2, and y3, but the factor B would contain levels to define the four points (se Examples).
#'  
#'  If one wishes to include other variables in the linear model, they should be indicated in the second formula, f2.
#'  This formula can be simply a right-hand formula, e.g., ~ x1 + x2 + x3 +...  Variables in this formula will typically
#'  be covariates that one wishes to include to account for extraneous sources of shape variation.  An analysis
#'  of variance (ANOVA) will be performed with type I sums of squares (SS) and a randomized residual permutation porcedure (RRPP).
#'  The variables in f2 will be added prior to the trajectory defining variables in f1.
#'
#'  Once the function has performed the analysis, a plot can be generated of the trajectories as visualized in the 
#'  original data or the space of principal components (PC1 vs. PC2). The first point in each trajectory is displayed as white, the 
#'  last point is black, and any middle points on the trajectories are in gray.  The colors of trajectories follow
#'  the order in which they are found in the dataset as a default, using R's standard color palette: black, red, green3,
#'  blue, cyan, magenta, yellow, and gray. However, one can override these colors with group.cols in plots using
#'  \code{\link{plot}}.  An additional plot argument (with default) is pca = TRUE, which uses the first two PCs
#'  for the trajectpry plot.  Changing this argument to FALSE means that the first two variables will be plotted, which is
#'  not ideal unless there is actually only two variables.  
#'  
#'  The function, \code{\link{summary}} can be used to provide an ANOVA summary plus pairwise statistics of a
#'  sn object of class "trajectory.analysis".  The argument, angle.type = c("r", "rad", "deg") can be used to
#'  toggle between vector correlations, vector angles in radians, or vector angles in degrees, respectively.
#'  
#' \subsection{Notes for geomorph 3.0}{ 
#' Previous versions of geomorph had two separate analytical approaches based on whether trajectories were estimated or
#' provided (as might be the case with motion trajectories; see Adams and Cerney 2007).  Starting with geomorph 3.0, 
#' commensurate analytical approaches are used.  This involves converting 1 x vp vectors for trajectpries, were p is the 
#' number of trajectory points and v is the number of variables in the data space, into p x v matrices, analagous to the
#' procedure for estimating trajectories.  Thus, rather than providing separate ANOVAs for size, oreintation, and shape
#' of trajectories, a general ANOVA is provided with pairwise statistics for the same attribute differences.  This change
#' does not compromise any interpretations made with previous versions of geomoprh, but enhances inferential capacity
#' by providing pairwise statistics and P-values.
#' }
#'
#' @param f1 A formula for the linear model, for trajectories (e.g., Y ~ A or Y ~ A * B)
#' @param f2 A formula for additional covariates  (e.g.,  ~ x1 + x2)
#' @param iter Number of iterations for significance testing
#' @param traj.pts An optional value specifying the number of points in each trajectory (if f1 contains a single factor)
#' @param data A data frame for the function environment, see \code{\link{geomorph.data.frame}} 
#' @param ... Arguments passed on to procD.fit (typically associated with the lm function)
#' @export
#' @keywords analysis
#' @author Dean Adams and Michael Collyer
#' @return If "estimate.traj=TRUE", the function returns a list with the following components: 
#'   \item{aov.table}{Procrustes ANOVA table.}
#'   \item{means}{The observed least squares means based on the linear model.}
#'   \item{random.means}{A list of matrices of means calculated in the random permutations.}
#'   \item{random.trajectories}{A list of all random means reconfigured as trajectories. The observed
#'   case is the first random set.}
#'   \item{path.distances}{The path distances of each observed trajectory.}
#'   \item{magnitude.diff}{A matrix of the absolute differences in pairwise path distances.}
#'   \item{trajectory.cor}{A matrix of pairwise vector correlations among trajectory PCs}
#'   \item{trajectory.angle.rad}{The vector correlations transformed into angles, in radians.}
#'   \item{trajectory.angle.deg}{The vector correlations transformed into angles, in degrees.}
#'   \item{trajectory.shape.dist}{A matrix of pairwise Procrustes distances among trajectory shapes.}
#'   \item{P.magnitude.diff}{P-values corresponding to trajectory magnitude differences.}
#'   \item{P.angle.diff}{P-values corresponding to angular differences in trajectories.}
#'   \item{P.shape.diff}{P-values corresponding to trajectory shape differences.}
#'   \item{Z.magnitude.diff}{Effect size of observed trajectory magnitude differences.}
#'   \item{Z.angle.diff}{Effect size of observed angular differences in trajectories.}
#'   \item{Z.shape.diff}{Effect size of observed trajectory shape differences.}
#'   \item{call}{The matched call.}
#'   \item{permutations}{The numer of random permutations used in the RRPP applied to the ANOVA 
#'   and trajectory statistics.}
#'   \item{trajectory.type}{A value of 1 is trajectories were provided or 2 if they were estimated.}
#' @references Collyer, M.L., and D.C. Adams. 2013. Phenotypic trajectory analysis: Comparison of 
#'  shape change patterns in evolution and ecology. Hystrix. 24:75-83.
#' @references Adams, D. C. 2010. Parallel evolution of character displacement driven by competitive 
#'   selection in terrestrial salamanders. BMC Evol. Biol. 10:1-10.
#' @references Adams, D. C., and M. M. Cerney. 2007. Quantifying biomechanical motion using Procrustes 
#'   motion analysis. J. Biomech. 40:437-444.
#' @references Adams, D. C., and M. L. Collyer. 2007. The analysis of character divergence along environmental 
#'   gradients and other covariates. Evolution 61:510-515.
#' @references Adams, D. C., and M. L. Collyer. 2009. A general framework for the analysis of phenotypic 
#'   trajectories in evolutionary studies. Evolution 63:1143-1154.
#' @references Collyer, M. L., and D. C. Adams. 2007. Analysis of two-state multivariate phenotypic change 
#'   in ecological studies. Ecology 88:683-692.
#' @examples
#' #1: Estimate trajectories from LS means in 2-factor model
#' data(plethodon) 
#' Y.gpa <- gpagen(plethodon$land)   
#' gdf <- geomorph.data.frame(Y.gpa, species = plethodon$species, site = plethodon$site)
#'
#' TA <- trajectory.analysis(coords ~ species*site, data=gdf, iter=499)
#' summary(TA, angle.type = "deg")
#' plot(TA, pca = TRUE)
#' 
#' # Retaining random values (first sets are always observed)
#' tra <- trajectory.analysis(Y.gpa~plethodon$species*plethodon$site,iter=15, verbose = TRUE)
#' tra$anova.table
#' tra$Random.values
#'
#' #2: Compare motion trajectories
#' data(motionpaths) 
#'
#' #Motion paths represented by 5 time points per motion 
#'
#' trajectory.analysis(motionpaths$trajectories~motionpaths$groups,
#' estimate.traj=FALSE, traj.pts=5,iter=15)
#' 
#' trajectory.analysis(motionpaths$trajectories~motionpaths$groups,
#' estimate.traj=FALSE, traj.pts=5,iter=15, verbose=TRUE)
trajectory.analysis <- function(f1, f2=NULL, iter=999, traj.pts = NULL, data = NULL){
  pfit1 <- procD.fit(f1, data=data, pca=FALSE)
  Terms <- pfit1$Terms
  dat <- pfit1$data
  datClasses <- sapply(model.frame(Terms, data=dat), function(x) data.class(x))
  if(any(datClasses != "factor")) stop("Only factors can be used on right hand side of first formula")
  if(length(datClasses) > 2 | length(datClasses) < 1) stop("Only one or two factors can be uesed in the first formula")
  if(length(datClasses) == 2 & ncol(attr(Terms, "factors")) != 3) stop("Two factors provided but no interaction is indicated in forst formula")
  if(length(datClasses) == 1 & is.null(traj.pts)) stop("If data are trajectories, the number of trajectory points must be defined")
  Y <- as.matrix(pfit1$Y)
  data <- geomorph.data.frame(data, Y=Y)
  f1 <- update(f1, Y~.)
  if(!is.null(f2)) {
    f2 <- update(f2, Y~.)
    full.terms <- c(attr(terms(f2), "term.labels"), attr(terms(f1), "term.labels"))
    ff <- as.formula(paste("Y~", paste(full.terms, collapse = "+"), sep=""))
  } else {
    full.terms <- pfit1$term.labels
    ff <- f1
  }
  if(length(datClasses) == 1 & is.null(f2)) {
    fr <- as.formula(Y~1)
  } else {
    red.terms <- full.terms[-length(full.terms)]
    fr <- as.formula(paste("Y~", paste(red.terms, collapse = "+"), sep=""))
  }
  pda <- procD.lm(ff, data=data, iter=iter, RRPP = TRUE)
  if(length(datClasses) == 1) pta <- traj.by.groups(ff, fr, traj.pts, data=data, iter=iter) else
    pta <- traj.w.int(ff, fr, data=data, iter=iter)
  if(length(datClasses) == 1) gp.names <- levels(pfit1$data[[length(pfit1$data)]]) else
    gp.names <- levels(pfit1$data[[length(pfit1$data)-1]]) 
  PD <- pta$PD[[1]]; names(PD) <- gp.names
  MD <- pta$MD[[1]]; Tcor <- pta$Tcor[[1]]; Tang <- pta$Tang[[1]]; SD <- pta$SD[[1]]
  rownames(MD) <- rownames(Tcor) <- rownames(Tang) <- rownames(SD) <- colnames(MD) <-
    colnames(Tcor) <- colnames(Tang) <- colnames(SD) <- gp.names
  P.MD= pta$P.MD
  P.angle = pta$P.angle
  P.SD = pta$P.SD
  Z.MD= pta$Z.MD
  Z.angle = pta$Z.angle
  Z.SD = pta$Z.SD
  rownames(P.MD) <- rownames(P.angle) <- rownames(P.SD) <- colnames(P.MD) <-
    colnames(P.angle) <- colnames(P.SD) <- gp.names
  out <- list(aov.table = pda$aov.table, 
              means = pta$means[[1]], 
              random.means = pta$means,
              random.trajectories = pta$trajectories,
              path.distances = PD, 
              magnitude.diff = MD,
              trajectory.cor = Tcor,
              trajectory.angle.rad = Tang,
              trajectory.angle.deg = Tang[[1]]*180/pi,
              trajectory.shape.dist = SD,
              P.magnitude.diff = P.MD,
              P.angle = P.angle,
              P.shape.diff = P.SD,
              Z.magnitude.diff = Z.MD,
              Z.angle = Z.angle,
              Z.shape.diff = Z.SD,              
              call = match.call(),
              permutations = iter+1,
              trajectory.type = length(datClasses)
              )
  class(out) <- "trajectory.analysis"
  out
}