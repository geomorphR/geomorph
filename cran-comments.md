## Resubmission
This is a patch release, 4.0.5.  It provides several new functions.

## Test environments
* local OS X install, R 4.0.5
* win-builder (devel and release)
* R-hub (all platforms)

## R CMD check results
There were no ERRORs or WARNINGs. 

## R-hub check results
Windows: OK
Ubuntu Linux: OK 
There was one error on Fedora Linux R-Devel

Error: No such container: geomorph_4.0.5.tar.gz-7b56b62b55bd422a8e26f32026bb9d97-3

This seems to be internal to R-hub and should have no bearing on the `geomorph 4.0.5` package.

## Downstream dependencies
I checked 11 reverse dependencies:
OK: 11
BROKEN: 0