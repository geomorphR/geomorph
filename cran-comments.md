## Resubmission
This is a patch release, 4.1.2 for minor bug fixes.

## Test environments
* local OS X install, R 4.6.1
* win-builder (devel and release)
* R-hub (all platforms)

## R CMD check results
There were no ERRORs or WARNINGs. 

## R-hub check results
Windows: OK
macos-arm64: OK
linux-R-devel: OK

## Downstream dependencies
I checked 7 reverse dependencies:
OK: 7 (NOTE: 'patternize' throws a note, which is not induced by geomorph. 
  The same issue is found on its current CRAN checks: 
  https://cran.r-project.org/web/checks/check_results_patternize.html.
  'patternize' was last updated in 2023.)
BROKEN: 0 