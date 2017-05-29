#' Calculate linear distances between landmarks
#' 
#' A function to calculate linear distances between a set of landmark coordinates (interlandmark distances)
#' 
#' Function takes a 3D array of landmark coordinates from a set of specimens and the addresses for the start 
#' and end landmarks defining linear measurements and then
#' calculates the interlandmark distances. The function returns a matrix of linear distances for all specimens.
#' If the 'lmks' matrix has row or column names defining the name of the linear measurements, the returned matrix will use
#' these for column names (see example). If only two interlandmark distances, 'lmks' input must be m x 2.
#' 
#' @param A A 3D array (p x k x n) containing landmark coordinates for a set of specimens
#' @param lmks A matrix or dataframe of landmark addresses for the start and end landmarks defining m linear measurements (can be either 2-x-m or m-x-2).
#' Either the rows or the columns should have names 'start' and 'end' to define landmarks.  
#' @export
#' @keywords utilities
#' @author Emma Sherratt & Michael Collyer
#' @return Function returns a matrix (n x m) of m linear distances for n specimens
#' @examples  
#' data(plethodon)
#' # Make a matrix defining three interlandmark distances 
#' lmks <- matrix(c(8,9,6,12,4,2), ncol=2, byrow=TRUE, 
#' dimnames = list(c("eyeW", "headL", "mouthL"),c("start", "end")))
#' # where 8-9 is eye width; 6-12 is head length; 4-2 is mouth length
#' # or alternatively
#' lmks <- data.frame(eyeW = c(8,9), headL = c(6,12), mouthL = c(4,2), 
#' row.names = c("start", "end")) 
#' A <- plethodon$land
#' lineardists <- interlmkdist(A, lmks)

interlmkdist <- function(A, lmks){
  if(!is.array(A)) {
    stop("Data matrix not a 3D array (see 'arrayspecs').") }
  lmks <- as.matrix(lmks)
  if(ncol(lmks) != 2 && nrow(lmks) != 2) 
    stop("Only one start and one end point are required to calculate distances")
  row.match <- match(c("start", "end"), rownames(lmks))
  col.match <- match(c("start", "end"), colnames(lmks))
  if(any(is.na(col.match)) && all(!is.na(row.match))) {
    lmks <- t(lmks)
    lmks <- lmks[,order(colnames(lmks), decreasing = TRUE)]
  }
  if(any(is.na(col.match)) && any(is.na(row.match))) {
    if(ncol(lmks) != 2 && nrow(lmks) == 2) {
      lmks <- t(lmks)
      cat("\nNo 'start' and 'end' points were defined.",
          "\nIt is assumed that matrix rows are appropriately ordered.\n","\n") 
    }
    if(ncol(lmks) == 2 && nrow(lmks) == 2) {
      if(is.null(rownames(lmks)) && is.null(colnames(lmks)))
        cat("\nNo 'start' and 'end' points were defined.",
            "\nNo names for distances were provided.",
            "\n'lmks' input assumed to be m x 2","\n\n") 
      if(!is.null(rownames(lmks)) && is.null(colnames(lmks))) {
        cat("\nNo 'start' and 'end' points were defined.",
            "\nIt is assumed that matrix columns are appropriately ordered.","\n\n")
      }
      if(is.null(rownames(lmks)) && !is.null(colnames(lmks))) {
        lmks <- t(lmks)
        cat("\nNo 'start' and 'end' points were defined.",
            "\nIt is assumed that matrix rows are appropriately ordered.","\n\n")
      } 
    }
  }
  lindist <- matrix(NA,ncol=nrow(lmks), nrow=dim(A)[3])
  if(!is.null(rownames(lmks))) colnames(lindist) <- rownames(lmks)   
  if(!is.null(dimnames(A)[[3]])) rownames(lindist) <- dimnames(A)[[3]] 
  for(i in 1:nrow(lmks)){lindist[,i] <- apply(A[lmks[i,],,], 3, dist)}
  return(lindist)
}
