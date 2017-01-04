#' Calculate linear distances between landmarks
#' 
#' A function to calculate linear distances between a set of landmark coordinates (interlandmark distances)
#' 
#' Function takes a 3D array of landmark coordinates from a set of specimens and a matrix 
#' of addresses for the start and end landmarks defining linear measurements and then
#' calculates the interlandmark distances. The function returns a matrix of linear distances for all specimens.
#' If the 'dists' matrix has row names defining the name of the linear measurements, the returned matrix will use
#' these for column names (see example).
#' 
#' @param A A 3D array (p x k x n) containing landmark coordinates for a set of specimens
#' @param dists A matrix (m x 2) of landmark addresses for the start and end landmarks defining m linear measurements
#' @export
#' @keywords utilities
#' @author Emma Sherratt
#' @return Function returns a matrix (n x m) of m linear distances for n specimens
#' @examples  
#' data(plethodon)
#' # Make a matrix defining three interlandmark distances 
#' dists <- matrix(c(8,9,6,12,4,2), ncol=2, byrow=T, dimnames = list(c("eyeW", "headL", "mouthL"),c("start", "end")))
#' # where 8-9 is eye width; 6-12 is head length; 4-2 is mouth length
#' A <- plethodon$land
#' lineardists <- interlmkdist(A, dists)

interlmkdist <- function(A, dists){
  if(!is.array(A)) {
    stop("Data matrix not a 3D array (see 'arrayspecs').") }
  if(dim(dists)[2] != 2) {
    stop("Only two landmarks (start and end) should be defined for interlandmark distances.") }
  lindist <- matrix(NA,ncol=nrow(dists), nrow=dim(A)[3])
  if(!is.null(rownames(dists))) colnames(lindist) <- rownames(dists)   
  if(!is.null(dimnames(A)[[3]])) rownames(lindist) <- dimnames(A)[[3]] 
  for(i in 1:nrow(dists)){lindist[,i] <- apply(A[dists[i,],,], 3, dist)}
  return(lindist)
}
