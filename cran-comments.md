## Resubmission
This is a patch release, 4.0.2.

## Test environments
* local OS X install, R 4.0.3
* win-builder (devel and release)
* R-hub (all platforms)

## R CMD check results
There were no ERRORs or WARNINGs. 

## R-hub check results
There was one error on Fedora Linux R-Devel

"Failed with error:  ‘there is no package called ‘shiny’’"

Shiny is not a dependency of 'geomorph' or any of its dependencies. This should have no bearing on the `geomorph 4.0.2` package.

## Downstream dependencies
I checked 10 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * I saw 0 new problems
 * I failed to check 0 packages
