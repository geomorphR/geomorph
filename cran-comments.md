## Resubmission
This is a patch release, 4.0.8.  It also provides several new functions.

## Test environments
* local OS X install, R 4.4.0
* win-builder (devel and release)
* R-hub (all platforms)

## R CMD check results
There were no ERRORs or WARNINGs. 

## R-hub check results
Windows: OK
Ubuntu Linux: OK 

Fedora Linux R-Devel (clang, gfortram) had RGL and X11 display issues:

1: In rgl.init(initValue, onlyNULL) : RGL: unable to open X11 display
2: 'rgl.init' failed 

This caused geomorph to not install on this platform. This seems to be internal to R-hub and should have no bearing on the `geomorph 4.0.8` package.

## Downstream dependencies
I checked 9 reverse dependencies:
OK: 9 
BROKEN: 0