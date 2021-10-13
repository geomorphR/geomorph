## Resubmission
This is a patch release, 4.0.1, that accompanies a minor update to the `RRPP` to version 1.1.0. The package, `RRPP`, is developed by the same authors as `geomorph`. `RRPP 1.1.0` and `geomorph 4.0.1` are to be concurrently updated and released.

## Test environments
* local OS X install, R 4.0.3
* win-builder (devel and release)
* R-hub (all platforms)

## R CMD check results
There were no ERRORs or WARNINGs. 

## R-hub check results
There was a persistent error resulting from an apparent R version 4.2 installation:

Error : Bioconductor does not yet build and check packages for R version 4.2; see
  https://bioconductor.org/install
  
This error appears to be insurmountable until `Bioconductor 3.14` is released.  This should have no bearing on the `geomorph 4.0.1` package.

## Downstream dependencies
I have also run R CMD check on downstream dependencies of geomorph. All packages that I could install passed. 
