## Resubmission
This is a patch release, 4.0.5.  It provides several new functions.

## Test environments
* local OS X install, R 4.0.5
* win-builder (devel and release)
* R-hub (all platforms)

## R CMD check results
There were no ERRORs or WARNINGs. 

## R-hub check results
There was one error on Fedora Linux R-Devel

Error : Bioconductor does not yet build and check packages for R version 4.2; see
  https://bioconductor.org/install

This error appears to be insurmountable until `Bioconductor 3.16` is released.  This should have no bearing on the `geomorph 4.0.4` package.

## Downstream dependencies
I checked 11 reverse dependencies:
OK: 11
BROKEN: 0