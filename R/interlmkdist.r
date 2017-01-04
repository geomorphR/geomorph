#' Calculate linear distances between landmarks
#' 
#' A function to calculate linear distances between a set of landmark coordinates (interlandmark distances)
#' 
#' Function takes a 3D array of landmark coordinates from a set of specimens and the addresses for the start 
#' and end landmarks defining linear measurements and then
#' calculates the interlandmark distances. The function returns a matrix of linear distances for all specimens.
#' If the 'dists' matrix has row or column names defining the name of the linear measurements, the returned matrix will use
#' these for column names (see example). If only two interlandmark distances, 'dists' input must be m x 2.
#' 
#' @param A A 3D array (p x k x n) containing landmark coordinates for a set of specimens
#' @param dists A matrix or dataframe of landmark addresses for the start and end landmarks defining m linear measurements (can be either 2-x-m or m-x-2)
#' @export
#' @keywords utilities
#' @author Emma Sherratt
#' @return Function returns a matrix (n x m) of m linear distances for n specimens
#' @examples  
#' data(plethodon)
#' # Make a matrix defining three interlandmark distances 
#' dists <- matrix(c(8,9,6,12,4,2), ncol=2, byrow=T, 
#' dimnames = list(c("eyeW", "headL", "mouthL"),c("start", "end")))
#' # where 8-9 is eye width; 6-12 is head length; 4-2 is mouth length
#' # or alternatively
#' dists <- data.frame(eyeW = c(8,9), headL = c(6,12), mouthL = c(4,2), 
#' row.names = c("start", "end")) 
#' A <- plethodon$land
#' lineardists <- interlmkdist(A, dists)

interlmkdist <- function(A, dists){
  if(!is.array(A)) {
    stop("Data matrix not a 3D array (see 'arrayspecs').") }
  dist <- as.matrix(dists)
  if(ncol(dists) != 2 && nrow(dists) == 2) dists <- t(dists)
  if(ncol(dists) == 2 && nrow(dists) == 2) {
    cat("'dists' input assumed to be m x 2","\n") }
  lindist <- matrix(NA,ncol=nrow(dists), nrow=dim(A)[3])
  if(!is.null(rownames(dists))) colnames(lindist) <- rownames(dists)   
  if(!is.null(dimnames(A)[[3]])) rownames(lindist) <- dimnames(A)[[3]] 
  for(i in 1:nrow(dists)){lindist[,i] <- apply(A[dists[i,],,], 3, dist)}
  return(lindist)
}
