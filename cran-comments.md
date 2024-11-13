## Resubmission
This is a patch release, 4.0.9.  It also provides several new functions.

## Test environments
* local OS X install, R 4.4.1
* win-builder (devel and release)
* R-hub (all platforms)

## R CMD check results
There were no ERRORs or WARNINGs. 

## R-hub check results
Windows: OK
macos: OK 

macos-arm64:
linux-R-devel:

1: Warning: 'rgl.init' failed, running with 'rgl.useNULL = TRUE'.

This caused geomorph to not install on Linux-R-devel and MacOS-Arm64 platforms. This seems to be internal to R-hub and 
should have no bearing on the `geomorph 4.0.9` package.

## Downstream dependencies
I checked 11 reverse dependencies:
OK: 11 
BROKEN: 0