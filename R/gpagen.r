#' Generalized Procrustes analysis of points, curves, and surfaces
#'
#' A general function to perform Procrustes analysis of two- or three-dimensional landmark data that 
#'  can include both fixed landmarks and sliding semilandmarks
#'
#' The function performs a Generalized Procrustes Analysis (GPA) on two-dimensional or three-dimensional
#'  landmark coordinates. The analysis can be performed on fixed landmark points, semilandmarks on 
#'  curves, semilandmarks on surfaces, or any combination. To include semilandmarks on curves, one 
#'  must specify a matrix defining which landmarks are to be treated as semilandmarks using the "curves=" 
#'  option. Likewise, to include semilandmarks 
#'  on surfaces, one must specify a vector listing which landmarks are to be treated as surface semilandmarks 
#'  using the "surfaces=" option. The "ProcD=TRUE" option will slide the semilandmarks along their tangent 
#'  directions using the Procrustes distance criterion, while "ProcD=FALSE" will slide the semilandmarks 
#'  based on minimizing bending energy. The aligned Procrustes residuals can be projected into tangent 
#'  space using the "Proj=TRUE" option. NOTE: Large datasets may exceed the memory limitations of R. 
#'
#'  Generalized Procrustes Analysis (GPA: Gower 1975, Rohlf and Slice 1990) is the primary means by which 
#'   shape variables are obtained from landmark data (for a general overview of geometric morphometrics see 
#'   Bookstein 1991, Rohlf and Marcus 1993, Adams et al. 2004, Zelditch et al. 2012, Mitteroecker and 
#'   Gunz 2009, Adams et al. 2013). GPA translates all specimens to the origin, scales them to unit-centroid 
#'   size, and optimally rotates them (using a least-squares criterion) until the coordinates of corresponding
#'   points align as closely as possible. The resulting aligned Procrustes coordinates represent the shape 
#'   of each specimen, and are found in a curved space related to Kendall's shape space (Kendall 1984). 
#'   Typically, these are projected into a linear tangent space yielding Kendall's tangent space coordinates 
#'   (Dryden and Mardia 1993, Rohlf 1999), which are used for subsequent multivariate analyses. Additionally, 
#'   any semilandmarks on curves and are slid along their tangent directions or tangent planes during the 
#'   superimposition (see Bookstein 1997; Gunz et al. 2005). Presently, two implementations are possible: 
#'   1) the locations of semilandmarks can be optimized by minimizing the bending energy between the 
#'   reference and target specimen (Bookstein 1997), or by minimizing the Procrustes distance between the two 
#'   (Rohlf 2010).  Note that specimens are NOT automatically reflected to improve the GPA-alignment.
#'
#' @param A An array (p x k x n) containing landmark coordinates for a set of specimens
#' @param Proj A logical value indicating whether or not the aligned Procrustes residuals should be projected 
#'   into tangent space 
#' @param ProcD A logical value indicating whether or not Procrustes distance should be used as the criterion
#'   for optimizing the positions of semilandmarks
#' @param PrinAxes A logical value indicating whether or not to align the shape data by principal axes 
#' @param curves An optional matrix defining which landmarks should be treated as semilandmarks on boundary 
#'   curves, and which landmarks specify the tangent directions for their sliding
#' @param surfaces An optional vector defining which landmarks should be treated as semilandmarks on surfaces
#' @param ShowPlot A logical value indicating whether or not a plot of Procrustes residuals should be displayed (calls \code{\link{plotAllSpecimens}})
#' @param ... Options to be passed to \code{\link{plotAllSpecimens}}
#' @keywords analysis
#' @export
#' @author Dean Adams
#' @return Function returns a list with the following components: 
#'   \item{coords}{A (p x k x n) array of aligned Procrustes coordinates, where p is the number of landmark 
#'     points, k is the number of landmark dimensions (2 or 3), and n is the number of specimens. The third 
#'     dimension of this array contains names for each specimen if specified in the original input array}
#'   \item{Csize}{A vector of centroid sizes for each specimen, containing the names for each specimen if 
#'     specified in the original input array}
#' @references  Adams, D. C., F. J. Rohlf, and D. E. Slice. 2004. Geometric morphometrics: ten years of 
#'    progress following the 'revolution'. It. J. Zool. 71:5-16.
#' @references Adams, D. C., F. J. Rohlf, and D. E. Slice. 2013. A field comes of age: Geometric 
#'   morphometrics in the 21st century. Hystrix.24:7-14.
#' @references Bookstein, F. L. 1991. Morphometric tools for landmark data: Geometry and Biology. 
#'  Cambridge Univ. Press, New York.
#' @references Bookstein, F. L. 1997. Landmark methods for forms without landmarks: morphometrics of 
#'   group differences in outline shape.  1:225-243.
#' @references Dryden, I. L., and K. V. Mardia. 1993. Multivariate shape analysis. Sankhya 55:460-480.
#' @references Gower, J. C. 1975. Generalized Procrustes analysis. Psychometrika 40:33-51.
#' @references Gunz, P., P. Mitteroecker, and F. L. Bookstein. 2005. semilandmarks in three dimensions. 
#'   Pp. 73-98 in D. E. Slice, ed. Modern morphometrics in physical anthropology. Klewer Academic/Plenum, New York.
#' @references Kendall, D. G. 1984. Shape-manifolds, Procrustean metrics and complex projective spaces. 
#'   Bulletin of the London Mathematical Society 16:81-121.
#' @references Mitteroecker, P., and P. Gunz. 2009. Advances in geometric morphometrics. Evol. Biol. 36:235-247.
#' @references Rohlf, F. J., and D. E. Slice. 1990. Extensions of the Procrustes method for the optimal 
#'   superimposition of landmarks. Syst. Zool. 39:40-59.
#' @references Rohlf, F. J., and L. F. Marcus. 1993. A revolution in morphometrics. Trends Ecol. Evol. 8:129-132.
#' @references Rohlf, F. J. 1999. Shape statistics: Procrustes superimpositions and tangent spaces. 
#'   Journal of Classification 16:197-223.
#' @references Rohlf, F. J. 2010. tpsRelw: Relative warps analysis. Version 1.49. Department of Ecology and 
#'   Evolution, State University of New York at Stony Brook, Stony Brook, NY.
#' @references Zelditch, M. L., D. L. Swiderski, H. D. Sheets, and W. L. Fink. 2012. Geometric morphometrics 
#'   for biologists: a primer. 2nd edition. Elsevier/Academic Press, Amsterdam.
#' @examples
#' #Example 1: fixed points only
#' data(plethodon) 
#' gpagen(plethodon$land,PrinAxes=FALSE)
#' gpagen(plethodon$land,pointscale=2)
#' 
#' #Example 2: points and semilandmarks on curves
#' data(hummingbirds)
#'
#' #Using Procrustes Distance for sliding
#' gpagen(hummingbirds$land,curves=hummingbirds$curvepts)   
#' 
#' #Using bending energy for sliding
#' gpagen(hummingbirds$land,curves=hummingbirds$curvepts,ProcD=FALSE)   
#'
#' #Example 3: points, curves and surfaces
#' data(scallops)
#' #Using Procrustes Distance for sliding
#' gpagen(A=scallops$coorddata, curves=scallops$curvslide, surfaces=scallops$surfslide) 
#' @useDynLib geomorph
gpagen<-function(A, Proj=TRUE,ProcD=TRUE,PrinAxes=TRUE,ShowPlot=TRUE,curves=NULL,surfaces=NULL,
                 ...){
  if (length(dim(A))!=3){
    stop("Data matrix not a 3D array (see 'arrayspecs').")  }
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first(see 'estimate.missing').")  }
  n<-dim(A)[3];   k<-dim(A)[2];  p<-dim(A)[1]
  nsurf<-0; ncurve<-0
  slided<-ifelse(ProcD==T,1,0)
  if(!is.null(curves)){ncurve<-nrow(curves) }
  if(!is.null(surfaces)){nsurf<-nrow(surfaces) }
  specs.size<-NULL
  for (i in 1:n)
  {specs.size[i]<-csize(A[,,i])[[1]]}
  if (is.null(curves) && is.null(surfaces)){  
    temp<-.C("DoGPA", as.integer(p),as.integer(k),as.integer(n),as.double(A),double(p*k*n),PACKAGE = "geomorph")[[5]]
    temp<-arrayspecs(matrix(temp,ncol=k,byrow=T),p,k)
  }
  else{ 
    temp<-.C("DoGPA1", as.integer(p),as.integer(k),as.integer(n),as.double(A),double(p*k*n),PACKAGE = "geomorph")[[5]]
    temp<-arrayspecs(matrix(temp,ncol=k,byrow=T),p,k)
    ref.gpa<-mshape(temp)	
    new.gpa<-.C("DoSlide", as.integer(slided), as.integer(p),as.integer(k),as.integer(n),
                as.double(temp),as.double(ref.gpa), double(p*k*n),as.integer(ncurve),as.integer(curves),as.integer(nsurf),
                as.integer(surfaces),PACKAGE = "geomorph" )[[7]]
    new.gpa<-arrayspecs(matrix(new.gpa,ncol=k,byrow=T),p,k)
    temp<-.C("DoGPA1", as.integer(p),as.integer(k),as.integer(n),as.double(new.gpa),double(p*k*n),PACKAGE = "geomorph")[[5]]
    temp<-arrayspecs(matrix(temp,ncol=k,byrow=T),p,k)
  }
  if(Proj==TRUE){temp<-orp(temp)  }
  dimnames(temp)[[3]]<-dimnames(A)[[3]]
  names(specs.size)<-dimnames(A)[[3]]
  if(PrinAxes==TRUE){
    ref<-mshape(temp); rot <- prcomp(ref)$rotation
    for(i in 1:k) if(sign(rot[i,i])!=1) rot[1:k,i] = -rot[1:k,i]
    for(i in 1:dim(temp)[[3]]){temp[,,i]<-temp[,,i] %*% rot }
  }
  if(ShowPlot==TRUE){ plotAllSpecimens(temp,...)}
  return(list(coords=temp,Csize=specs.size))
}
