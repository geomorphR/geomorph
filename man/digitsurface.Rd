% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/digitsurface.r
\name{digitsurface}
\alias{digitsurface}
\title{Digitize 3D fixed landmarks and surface semilandmarks}
\usage{
digitsurface(spec, fixed, ptsize = 1, center = TRUE)
}
\arguments{
\item{spec}{An object of class shape3d/mesh3d, or matrix of 3D vertex coordinates}

\item{fixed}{Either a numeric value designating the number of fixed landmarks to be selected by \code{\link{digit.fixed}}, or a matrix of 3D coordinates collected previously}

\item{ptsize}{Size of mesh points (vertices), e.g. 0.1 for dense meshes, 3 for sparse meshes}

\item{center}{Should the object 'spec' be centered prior to digitizing?}
}
\value{
Function returns (if assigned to an object) and writes to the working directory an NTS
 file, containing the landmark coordinates. The file name corresponds to the name of the specimen.
}
\description{
An interactive function to digitize three-dimensional (3D) landmarks on a surface lacking known landmarks.
Input for the function is either a matrix of vertex coordinates defining a 3D surface object
or a mesh3d object as obtained from \code{\link{read.ply}}.
}
\details{
Function for digitizing fixed 3D landmarks and placing surface sliding semilandmarks using a previously created 
template. Following the selection of fixed points (see digitizing below), the function finds surface semilandmarks 
following the algorithm outlined in Gunz et al. (2005) and Mitteroecker and Gunz (2009). digitsurface finds the 
same number of surface semilandmarks as the template (created by \code{\link{buildtemplate}}) by downsampling the 
scanned mesh, after registering the template with the current specimen via GPA. A nearest neighbor algorithm is 
used to match template surface semilandmarks to mesh points of the current specimen. To use function digitsurface, 
the template must be constructed first, and 'template.txt' be in the working directory. Because template matching 
is based on the correspondence of fixed landmark points in the template and the specimen, a minimum of four fixed 
landmarks must be used. 


For details on the full procedure for digitizing fixed 3D landmarks and surface
sliding semilandmarks, see the relevant vignette by running \code{vignette("geomorph.digitize3D")}.
 
 NOTE: Function centers the mesh before digitizing by default (center = TRUE). If one chooses not to center,
 specimen may be difficult to manipulate in rgl window.

\subsection{Digitizing}{
Digitizing of fixed landmarks is interactive. Once a point is selected, the user is asked if the system should 
keep or discard the selection (y/n). If "y", the user is asked to continue to select the next landmark. If "n" 
the removes the last chosen landmark, and the user is asked to select it again. This can be repeated until the 
user is comfortable with the landmark chosen. 

To digitize with a standard 3-button (PC):
\enumerate{
 \item the RIGHT mouse button (primary) to select points to be digitized,
 \item the LEFT mouse button (secondary) is used to rotate mesh, 
 \item the mouse SCROLLER (third/middle) is used to zoom in and out.
}
NOTE: Digitizing functions on MACINTOSH computers using a standard 3-button mice works as specified. Macs using platform 
specific single button mice, XQuartz must be configured: go to Preferences > Input > tick "Emulate three button mouse":
\enumerate{
 \item press button to rotate 3D mesh,
 \item press button while pressing COMMAND key to select vertex to be used as a landmark,
 \item press button while pressing OPTION key to adjust mesh perspective.
 \item the mouse SCROLLER or trackpad two finger scroll is used to zoom in an out.
 }

NOTE: there is no pan (translate) functionality in rgl library for all platforms at this time.
}

\subsection{AUTO mode}{ 
The function as described above (for interactive mode) calls \code{\link{digit.fixed}}, prompting the user to select fixed landmarks
in the rgl window. However if the user has digitized these fixed landmark elsewhere (e.g., in other software), then the input for
parameter 'fixed' can be a p-x-k matrix of 3D coordinates. In this case, the function the function will automatically use these
landmarks and fit the template of sliding semilandmarks.
}
}
\references{
Gunz P, Mitteroecker P, & Bookstein FJ (2005) Semilandmarks in Three Dimensions. Modern Morphometrics in Physical Anthropology, ed Slice DE (Springer-Verlag, New York), pp 73-98.

Mitteroecker P & Gunz P (2009) Advances in Geometric Morphometrics. Evolutionary Biology 36(2):235-247.
}
\seealso{
\code{\link{buildtemplate}}

\code{\link{read.ply}}

\code{\link{digit.fixed}}

\code{\link[rgl]{rgl-package}} (used in 3D plotting)
}
\author{
Erik Otarola-Castillo & Emma Sherratt
}
\keyword{digitizing}
