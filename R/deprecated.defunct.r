
#' Deprecated functions in geomorph
#'
#' The following function has been deprecated in geomorph
#'
#' This function has been deprecated. Use \code{\link{procD.lm}} instead, along with \code{\link{anova.lm.rrpp}} and \code{\link{pairwise}}.  
#' See examples in \code{\link{procD.lm}} 
#' and/or the vignette, 'Addition by subtraction: improving geomorph capabilities with fewer functions.'
#' 
#' @export
advanced.procD.lm <- function(){
  .Deprecated("procD.lm")
}

#' Deprecated functions in geomorph
#'
#' The following function has been deprecated in geomorph
#'
#' This function has been deprecated. Use procD.lm, #' advanced.procD.lm now deprecated: use \code{\link{procD.lm}} instead, along with \code{\link{anova.lm.rrpp}} and \code{\link{pairwise}}.  
#' See examples in \code{\link{procD.lm}} 
#' and/or the vignette, 'Addition by subtraction: improving geomorph capabilities with fewer functions.'
#' 
#' @export
procD.allometry <- function(){
  .Deprecated("procD.lm")
}

#' Deprecated functions in geomorph
#'
#' The following function has been deprecated in geomorph
#'
#' This function has been deprecated. Use procD.lm, followed by \code{\link{anova.lm.rrpp}}, instead, with error adjustment.  
#' See examples in \code{\link{procD.lm}}
## and/or the vignette, 'Addition by subtraction: improving geomorph capabilities with fewer functions.'
#' 
#' @export
nested.update <- function(){
  .Deprecated("anova.lm.rrpp")
}




