## Resubmission
This is a minor update, version 3.3.0, that fixes a few bugs in version 3.2.1, and adds a new function. 

## Test environments
* local OS X install, R 3.6.0
* win-builder (devel and release)
* R-hub (all platforms)

## R CMD check results
There were no ERRORs or WARNINGs. 

## Downstream dependencies
I have also run R CMD check on downstream dependencies of geomorph. All packages that I could install passed except:

* Phylocurve: This is due to a deprecated function in geomorph. The maintainer of phylocurve has been made aware of the issue, I have provided the maintainer with code that resolves this, and have given the maintainer time to make their changes. Additionally, since geomorph 3.2.1 (Jan 27, 2020) the man page for the deprecated function stated: "Notes for geomorph 3.2.1.0900 and subsequent versions. This function is deprecated and will soon be removed. The function 'gm.prcomp' can be used instead to generate phylomorphospace plots."
