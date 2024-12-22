# CHANGES IN GEOMORPH VERSION 4.0.10

### NEW FEATURES

* Added new K-components output and plots to `physignal.eigen`.

### BUG FIXES

* Updated handling of p > n conditions in `physignal.eigen`.

# CHANGES IN GEOMORPH VERSION 4.0.9

### NEW FEATURES

* New function 'extended.pgls' for evaluating within-species trends in a phylogenetic context
* New dataset 'pupfish.ws'
* New function `physignal.eigen` for evaluating the dimensions of phylogenetic signal in multivariate data
* Added an option to `gm.prcomp` to allow the phylogeny to be scaled to unit length
* Added `groups.first` argument to `gm.measurement.error`, to optionally include groups in ANOVA table.

### BUG FIXES

* Fixed missing argument in `perm.index` within `morphol.disparity`.
* Ensured proper size-scaling of predicted shapes in `gm.prcomp`.
* Fixed the missing re-order of data for `anc.BM`.

# CHANGES IN GEOMORPH VERSION 4.0.8

### NEW FEATURES

* Improved `anc.BM`

### OTHER CHANGES

* Internal changes for compliance with RRPP and Matrix changes

### Bug fixes
* Fixed issue with `integration.Vrel` and phylogeny checks
* Fixed issue with unsigned.AI in `bilat.symmetry` and removed signed.AI
* Fixed issue with `readland.tps` and negNA = TRUE
* Fixed issue with `plot.gm.prcomp` `phylo.par` default

# CHANGES IN GEOMORPH VERSION 4.0.7

### NEW FEATURES

* New function, `gm.measurement.error`

### Bug fixes
* Fixed issue with `estimate.missing` for shapes files and TPS option
* Fixed names issue with `bilat.symmetry`
* Fixed problems with ellipsis use in `readmulti.tps`
* Fixed issue with `lambda.opt` internal function, for use in `physignal.z`, which precluded verbose results.

### OTHER CHANGES

* Reversed order of asymmetry components for internal consistency
* Updated internal code and Rd file formats to match new CRAN policies
* Added a feature to print a table to identify inconsistent landmarking in both `readland.tps` and `readmulti.tps`, rather than just return an error.

# CHANGES IN GEOMORPH VERSION 4.0.6

### Bug fixes
* Fixed issue with `readland.tps` when `warnmsg = FALSE`
* Fixed issue with `bilat.symmetry` (updated nested model computations)
* Fixed issue with `readland.tps` when curves exist but are not to be read
* Fixed issue with bootstrap CI for `modularity.test` and `phylo.modularity`
* Fixed a bug in `logLikh` internal function that required `matrix` class objects for calculation.
* Fixed a typo bug in `procD.lm` that caused lost results.
* Fixed names being dropped from shape components in `bilat.symmetry`
* Fixed `readland.fcsv` to be consistent with SlicerMorph updates.

### OTHER CHANGES
* Additional error checking in `modularity.test` and `phylo.modularity` `rgl`
* Added better warning for `estimate.missing` for failing with `method = "Reg"` to estimate landmarks.
* Added `shape.predictor` examples to PLS functions.

# CHANGES IN GEOMORPH VERSION 4.0.5

### NEW FEATURES

* New function, `physignal.z`
* New function, `compare.physignal.z`

### OTHER CHANGES
* Updated `plotspec` and `picknplot.shape` functions to accommodate changes in `rgl`

# CHANGES IN GEOMORPH VERSION 4.0.4

### NEW FEATURES

* New function, `integration.Vrel`
* New function, `compare.ZVrel`
* Added a `lambda` argument to `procD.pgls`.

# CHANGES IN GEOMORPH VERSION 4.0.3

### OTHER CHANGES
* Updated `Cov.proj` function to work with `RRPP`

# CHANGES IN GEOMORPH VERSION 4.0.2

### BUG FIXES
* Fixed issue in `compare.CR` when CR.null = TRUE.
* Added catch to `anc.BM` for singleton nodes

# CHANGES IN GEOMORPH VERSION 4.0.1

### BUG FIXES
* Fixed issue in `compare.CR` when CR.null = TRUE.
* Added catch to `anc.BM` for singleton nodes

# CHANGES IN GEOMORPH VERSION 4.0.1

### NEW FEATURES
* Added a transform argument to `morphol.disparity` for PGLS models, consistent with `gm.prcomp` approaches.
* New function, `module.eigen`, plus S3 generic associated functions.
* New function, `na.omit.geomorph.data.frame`.

### BUG FIXES
* Fixed some `as.matrix` names dropping in support code.
* Fixed pairwise Z-scores in `compare.pls` and `compare.CR` to use box-cox transformed values.
* Added Pcov as needed output `phylo.integration` and `phylo.modularity`.
* Fixed small issue with label plotting in `gridPar.r`.
* Fixed bug in `readmulti.tps` to allow for any kind of tps input and correctly combine them into an array
* Fixed typo bug in `integration.test` (only affected separate matrices).
* Fixed bug in `estimate.missing` (only affecting method = "Reg").

# CHANGES IN GEOMORPH VERSION 4.0.0 

### NEW FEATURES
* Included calculation of individual asymmetry indices to `bilat.symmetry` output
* Updated `bilat.symmetry` to use arguments available with updates to `RRPP::lm.rrpp`, taking advantage of parallel processing, C++, and turbo-charged (coefficient-suppressed) calculations to maximize speed.
* More options for parallel processing.
*`make_ggplot2` function for converting `geomorph` plots to `ggplot` objects, which can be amended.

### BUG FIXES
* Fixed bug in `plotAllometry` where `crossprod` was used instead of `tcrossprod`
* Fixed typo causing error in `phylo.integration`
* Fixed issue in `integration.test` for permutations

# CHANGES IN GEOMORPH VERSION 3.3.2 (Minor Release)

### NEW FEATURES
* Added Box-Cox transformation to `effect.size`
* Added an argument to flip axes in `plot.gm.prcomp`
* Added pairwise r, Z, and P-values to `integration.test` and `phylo.integration`

### BUG FIXES
* Fixed bug in `compare.CR` so object labels are retained
* Fixed bug with `plot.bilat.symmetry` that missed changed object output in `bilat.symmetry`.
* Fixed issue with permutations in multi-module configurations for `integration.test` and `phylo.integration`
* Fixed error in `plot.gm.prcomp` 3D tree plotting with time.plot = TRUE
* Fixed error in `compare.evol.rates` for univariate data and permutations
* Fixed plotting parameter arguments for phylo.par in `plot.gm.prcomp`.
* Made it possible to assign `summary.gm.prcomp` as an object
* Fixed plotting issue in `warpRefMesh`
* Fixed single variable issue in `compare.evol.rates` related to R 4.0.2 changes
* Added factor labels to output object CR.mat in `modularity.test` & `phylo.modularity`
* Fixed bug in computations of `phylo.modularity`
* Fixed order consistency for pairwise calculations in `phylo.modularity` 
* Allows single landmarks as slider LM in gpagen
* Fixed labels in `plotRefToTarget`, allowing 3d TPS warp grids to show labels and for txt.pos to be passed to 3d plots		
* Fixed bug in txt.adj argument for `gridPar`

### OTHER CHANGES
* Updated `readmulti.nts` function turning it into a wrapper of `readland.nts` thus allowing also multiple dta files to be read in and compiled


# CHANGES IN GEOMORPH VERSION 3.3.1 (Patch Release)

### BUG FIXES
* Fixed non-centering issue in `geomorph:::pls`

# CHANGES IN GEOMORPH VERSION 3.3.0 (Minor Release)

### NEW FEATURES
* A vignette has been added (geomorphPCA) to aid the transition of users to the gm.prcomp family of functions for exploring and visualizing shape space.
* Updates to `gm.prcomp` to include PaCA.
* Added the possibility of a 3D PCA plot with a phylogeny and time on the z-axis to `plot.gm.prcomp`.
* Added the possibility of controlling plot3d definitions in `plotspec`.
* New dataset of dorsal views of lizard heads added: `lizards`.
* `mshape` now has options (`na.action` argument) for missing values.
* New function, `readland.fcsv` for importing landmark data from SlicerMorph `.fcsv` files.
* `combine.subsets` updated to include weighting options for relative centroid size.

### BUG FIXES
* Fixed row order issue for BM simulation in `compare.evol.rates` and `compare.multi.evol.rates`
* Fixed bug in FA shape components output in `bilat.symmetry`
* Fixed pGPA issue to use first specimen rather than mean in first iteration

### OTHER CHANGES
* `plotGMPhyloMorphospace` and `plotTangentSpace` now deprecated
* Removed DA.mns and FA.mns from `bilat.symmetry` (use DA.component and FA.component)
* Option for handling NAs added to mshape and plot.mshape
* Function plotTangentSpace has been deprecated 

# CHANGES IN GEOMORPH VERSION 3.2.1.0900 (Minor Release)

### NEW FEATURES
* Added the possibility of plotting an outline with method = "vector" in plotRefToTarget

### BUG FIXES
* Fixed small bug in plot.mshape par definitions
* Amended specimen name extraction from .tiff files when using imageID in readland.tps

### OTHER CHANGES
* Function plotGMPhyloMorphoSpace has been deprecated 

# CHANGES IN GEOMORPH VERSION 3.2.1

### NEW FEATURES
* New function added: `readmulti.tps` 

### BUG FIXES
* Fixed small bug in readland.tps to accurately replace negative values with NA when only one specimen has NAs and allow for NAs in the original file

# CHANGES IN GEOMORPH VERSION 3.2.0 (Minor Release)

### NEW FEATURES
* Effect sizes provided in the output of `compare.evol.rates`, `compare.multi.evol.rates`,`integration.test`,
`modularity.test`, `phylo.integration`, `phylo.modularity`, `physignal`, `two.b.pls`

### OTHER CHANGES
* Adjustment in 'class' statements to be compliant with new handling of objects in R 4.0

### BUG FIXES
* Fixed small bug in plotAllometry that returned spurious RegScore values.
* Adjusted calculations of pairwise effect size in compare.pls and compare.CR to mirror publication equations
* Fixed issue with 3D plotting in plotOutliers with inspect.outliers = T

# CHANGES IN GEOMORPH VERSION 3.1.3 (Patch Release)

### NEW FEATURES
* New function added: `compare.CR` 
* New function added: `readmulti.tps` 

### OTHER CHANGES
* Adjustment to getSurfPCs to allow for more robust surface sliding
* Added a one-tail/two-tailed option to compare.pls

### BUG FIXES
* Fixed `mshape` inability to distinguish between one p x k specimen and a n x pk matrix.
* Fixed issue with 3D TPS grids in plotRefToTarget
* Fixed issue with digit.fixed, digitsurface, plotspec, warpRefMesh, build.template, for null ply $material, due to rgl updates
* Updated example in read.ply to not give bad advice for ply $material
* Fixed issue with 3D plotting in plotOutliers with inspect.outliers = T


# CHANGES IN GEOMORPH VERSION 3.1.2 (Patch Release)

### NEW FEATURES
* Three vignettes have been added: one assists the user in implementing analyses previously included in the deprecated functions `procD.allometry`, `advanced.procD.lm` and `nested.update` through `procD.lm` and other `RRPP` tools; the second provides a comprehensive list of main geomorph functions and their use; and the third is a step-by-step demonstration of 3D digitizing procedures
* `readland.tps` now includes an argument (negNA = FALSE) to define whether negative landmark coordinates should be transformed to NAs (by setting the argument to TRUE)      

### OTHER CHANGES

### BUG FIXES 
* Fixed bug with `combine.subsets` for gpagen objects
* Tweaked `procD.lm` to work better with missing data frames
* Fixed p-value miscalculation in `bilat.symmetry`
* Added meshColor specification in plotting of 3d meshes (to match RGL change)
* Fixed text plot for 3D (to match RGL changes)
* Added missing gls component in `procD.lm` output
* Fixed missing SS.types for pgls

------

# CHANGES IN GEOMORPH VERSION 3.1.1 (Patch Release)

### NEW FEATURES

### OTHER CHANGES

### BUG FIXES
* Fixed single landmark issue in `readland.shapes`
* Updated simulator for `compar.evol.rates` to utilize all simulated datasets
* Updated rotation in `fixed.angle`
* Removed code causing warning in `GMfromShapes1`

------

# CHANGES IN GEOMORPH VERSION 3.1.0 (Minor Release)

### NEW FEATURES
* Permutation computations now performed in `RRPP` 
* New functions added: `plotAllometry`, `gm.prcomp`, `picknplot.shape`, `shape.hulls`
* Deprecated: `advanced.procD.lm`, `nested.update`, `procD.allometry` (features of these functions now found using `procD.lm`)
* Generalized plotting options
* Reduced dependencies on `geiger` and `Matrix`
* Added mean-centering to `compare.evol.rates` for permutations
* Updated `fixed.angle` to allow centroids from multiple points to be used as angle points.
* Updated `coords.subset` to avoid issue with arrayspecs and naming subsetted data.
* Changed the default for sliding to ProcD=FALSE in `gpagen`
* Added Procrustes distance matrix to `gpagen` output
* Moved `trajectory.analysis` to `RRPP` package and updated its arguments
* Added partial disparity option to `morphol.disparity`

### OTHER CHANGES

### BUG FIXES
* Fixed ratio call for permutations in `compare.evol.rates`
* Fixed bug in `phylo.modularity` call when groups > 2
* Fixed scaling default in `combine.subsets`
* Updated `trajectory.analysis` to properly center mean values for PC plotting (now in RRPP)
* Added a trap for large distance matrices to `gpagen`, to avoid errors
* Fixed bug to allow for a single specimen with missing data in `readland.tps`
* Fixed bug for HOS result mismatching criterion for `procD.allometry`
  
------
    
# CHANGES IN GEOMORPH VERSION 3.0.7 (Patch Release)

### NEW FEATURES
* Computations for `advanced.procD.lm` now performed using RRPP package
* Modified output in `bilat.symmetry` so output shapes retain order in 'ind' factor
* Modified `readland.tps` to identify missing data only when all lm coordinate are <0 and interactively prompt the user to confirm if they are to be treated as NAs  

### OTHER CHANGES

### BUG FIXED

------

# CHANGES IN GEOMORPH VERSION 3.0.6 (Patch Release)

### NEW FEATURES
* `plotOutliers` now allows to plot outlier configurations in order to compare their shape with the consensus     
* `plot.mshape` function added to plot the consensus configuration with numbers and links
* New function `gm.prcom`p implements raw and weighted PCA and allows S3 generic plotting of its output
* New function `readland.shapes` allows reading a shapes file produced by StereoMorph, including landmark data and (potentially multiple) curves, and sampling semilandmarks from these curves.
* New function `picknplot.shape` allows interactively picking points in geomorph scatterplots to visualize shape variation across morphospace.

### OTHER CHANGES
* Modified previous addition to `readland.tps` to identify missing data only when all lm dimensions are <0 and interactively prompt the user to confirm if they are to be treated as NAs.  
* Added ellipsis to `bilat.symmetry` and `geomorphShapes` option to allow more flexible `gpagen` options within
* Changed `procD.pgls` to use residuals from GLS model for permutations
* Added effect type option to `advanced.procD.lm`
* Added an additional line of summary to ANOVA tables to signal how effect sizes are calculated
* Adjusted ANOVA tables to show from which distributions P-values are estimated
* Modified `plot.procD.allometry` to allow direct control of all plotting arguments by the user
* Modified `plotRefToTarget` and tps to use plot.xy instead of plot, and thus avoid conflicts with `plot.mshape`
* Added catch for missing data to `readland.tps`, to identify negative values and recode them as NAs

### BUG FIXES
* Fixed multi-group output in `phylo.modularity`
* Fixed issue with `compare.evol.rates` permutation method for univariate trait
* Fixed issues with `estimate.missing` regression approach
* Fixed issues with `phylo.integration` where Y-dataset contained a single variable
* Fixed issue with `integration.test` when 3+ partitions with non-contiguous variables
* Fixed effect size calculation error in `procD.lm` for single factor OLS models. 
* Fixed ellipses options in `procD.fit`
    
------

# CHANGES IN GEOMORPH VERSION 3.0.5  (Patch Release)

### NEW FEATURES

### OTHER CHANGES
* Simplified options available in `procD.allometry`

### BUG FIXES
* Fixed CAC calculations in `procD.allometry`
* Fixed dimnames error in `readland.nts` when file contained a single specimen's data
* Fixed PGLS weights for factors in `morphol.disparity` when phylogeny utilized
* Fixed `readland.tps` to be general to whitespace delimitation

------

# CHANGES IN GEOMORPH VERSION 3.0.4

### NEW FEATURES
* New function: `interlmkdist` to calculate linear distances between landmarks (interlandmark distances)
* PGLS option for `advanced.procD.lm`
* An option to use types I, II, or III sums of squares and cross-products for ANOVA analyses (e.g., `procD.lm`, `procD.pgls`)
* An option to choose between SS, F-values, or Cohen's f-squared to calculate effect sizes for ANOVA analyses
* Effect sizes are now centered on mean values from distributions of random statistics.  Statistics are also log-transformed in certain cases to normalize distributions
* An option to use only an intercept as a model for `procD.lm` and `procD.pgls`, with limited output
* An option to use permutation tests when comparing net evolutionary rates for `compare.evol.rates`

### OTHER CHANGES
* Added ellipsis programming to `plotTangentSpace` to pass arguments to `prcomp`.  Also, adjusted default tolerance to remove redundant PC dimensions.
* Enhanced landmark labeling in `arrayspecs`
* Created transferability of dimnames between `arrayspecs` and `two.d.array`
* Updated Rd files for PLS functions to better cross-reference the similar functions
* Added data frame output to `procD.lm`
* Change to `plotGMPhyloMorphoSpace` to center data by phylogenetic mean
* `read.morphologika` now supports lists of file names to return a single data object
* Updated Rd files for `readland.nts` and `readmulti.nts` to avoid ambiguity and resolve misinformation in previous versions
* Updated `readland.nts` to accept specimen labels with spaces in name

### BUG FIXES
* Fixed problem with incorporating lm arguments in `procD.fit`
* Fixed univariate PLS to allow negative correlations
* Fixed permutation issues with `advanced.procD.lm`
* Fixed links not plotting in `plotAllSpecimens` if no colour specified
* Added missing random.r output to `integration.test` and `phylo.integration`

------

# CHANGES IN GEOMORPH VERSION 3.0.3 (Patch Release)
  
### NEW FEATURES
* New functions: `compare.pls`, `coords.subset`, `shape.predictor`
* `droplevels.geomorph.data.frame` added to support code

### OTHER CHANGES
* Added a sensor to support code for `procD.lm` and its allies to 
    choose the computationally fastest algorithms based on design matrix complexity
    and data dimensionality
* Updated `procD.fit` to remove unused levels from factors
* Updated all functions using `geomorph.data.frame` to drop unused factor levels
* Updated `mshape` to be used on lists, arrays, or matrices

### BUG FIXES  
* Fixed `single.factor` function to properly maintain factor levels when combining factors
* Fixed plotting issues with `phylo.integration`

------

# CHANGES IN GEOMORPH VERSION 3.0.2 (Patch Release)

### NEW FEATURES
* Added return of min/max shape matrices to `plot.procD.allometry` and `plot.pls`
* Added option to save specimen ID to ID= line in `writeland.tps`

### OTHER CHANGES
* Optimized analytical functions for faster computations
* Added force match of specimen order between blocks in `two.b.pls`
* Updated A, Ahat, ref, p and k outputs when 2D matrix is used in `procD.allometry`
* Added option to plot group labels in `plot.procD.allometry` (method="PredLine")
* Added legend option to `plotTangentSpace` when groups are specified
* Removed 'verbose' from `plotTangentSpace`; function now returns PC scores automatically when assigned to object
* Updated `plotTangentSpace` to return min and max shapes for all PC axes in $pc.shapes
    
### BUG FIXES
* Fixed label output issue with `bilat.symmetry` shape components
* Fixed error plotting TPS grids when groups were included in `plot.procD.allometry`
* Fixed error in `digit.curves` where open outline would be treated as closed if starting point was the end point
* Fixed reading of ape::chronos phylogenies in 'phylo' functions
* Fixed error in reading names in `readland.tps` for some tps files

------
    
# CHANGES IN GEOMORPH VERSION 3.0.1 (Patch Release)


### NEW FEATURES

### OTHER CHANGES
* Removed dependency on `phytools`
* Small addition to outputs in `procD.pgls`
* Added reading polygons in `read.morphologika`
* Plots of 3D shapes now plotted in pairs in single rgl window for: `plotTangentSpace`, `plot.bilat.symmetry`, `plot.procD.allometry`, `plot.pls`
* Functions performing warping of 3D meshes now present progress bar
* Added option pt.col to `plot.procD.allometry` to designate plotting colors

### BUG FIXES
* Fixed error in output order in `plotOutliers`
* Fixed error in calculation slope distances in `advanced.procD.lm`  
* Fixed 2D data input error for `phylo.modularity`
* Fixed centering error for CAC scores in `procD.Allometry`
* Fixed deformation plot issue with `phylo.integration` (via `plot.pls`)
* Fixed two small issues with `plot.pls`, regarding best fit line and matrix reduction.
* Fixed three small issues with `gpagen` source files: indexing errors, arbitrary PC rotations for surface points, maximum iteration disparity.
* Fixed `digitize2D` scaling issue when different scales used in each image
* Fixed `phylo.integration` error when 3+ partitions examined
* Fixed typo in `modularity.test` for CI intervals into matrix input
* Updated permutation procedure in `phylo.integration` so that prob(A,B|phy)~prob(B,A|phy)
* Fixed typo in `plotAllSpecimens` where links were not being plotted

------

# CHANGES IN GEOMORPH VERSION 3.0.0 (Major Release)

### NEW FEATURES
* New functions: `modularity.test`, `integration.test`, `phylo.modularity`, 
        `phylo.integration`, `procD.allometry`, `nested.update`, `geomorph.data.frame`
* Seed option added for most analytical functions
* Major overhaul of underlying support code for analytics
* Output from analytical functions provided in lists
* Added 'summary' option to most analytical functions
* Plotting for most analytical functions as 'plot(res)' if res is output from function
* Diagnostic plots for anova/regression added
* Nested ANOVA via `nested.update` of `procD.lm` objects

### OTHER CHANGES
* Deprecated: `compare.modular.partitions`, `morphol.integr`, `phylo.pls`, and `plotAllometry` 
* Removed internal C-code for `gpagen`
* Removed automatic plots for most analytical functions
* Removed method="" parameter from `physignal`. Only Kmult used
* Removed verbose option for most functions: all output provided in lists
* Added option to pre-multiply coordinates by scale in `digitize2d`
* Added plotting options to `plotAllSpecimens`

### BUG FIXES
* Fixed `morphol.disparity$Prob.Disp` displaying NAs

------

# CHANGES IN GEOMORPH VERSION 2.1.7-1 (Patch Release)


### NEW FEATURES

### OTHER CHANGES

### BUG FIXES
* Fixed small change in C-code for `gpagen`

------

# CHANGES IN GEOMORPH VERSION 2.1.7 (Patch Release)

### NEW FEATURES
* New function `compare.multi.evol.rates` for comparing rates of evolution among traits
* `plotGMPhyloMorphoSpace` now plots 3D phylomorphospaces and chronophylomorphospaces
* Added ShowPlot option to `bilat.symmetry`, `compare.modular.partitions`, `globalIntegration`, `morphol.integr`, and `phylo.pls`

### OTHER CHANGES
* Updated Imports packages as per new CRAN policies
* Small change to C-code for `gpagen` (added additional checks on alignment)
* Package ape now full Import for geomorph (no longer ImportsFrom)
* Enhanced input flexibility for `advanced.ProcD.lm`: for single-factor analyses and matrix/variable input
* Enhanced `readland.nts` flexibility with specimen labels; now supports spaces in labels
* Generalized `read.ply` to allow reading meshes with many properties
* Generalized input for `physignal` and `compare.evol.rates`: univariate data accepted as named vector
* Enhanced input for `define.links`: read and append links to existing links matrix

### BUG FIXES
* Corrected error in `warpRefMesh` where normals were incorrectly assigned to new mesh3d object
* Corrected error in `readland.tps` which read in a file containing a single specimen returned a 2D matrix rather than 3D array; this fixes the issue with `digitize2d` not working for a single file
* Corrected error in `morphol.integr` where warpgrids = F did not work for 3D datasets
* Corrected issue with `globalIntegration` for use with 3D data
* Corrected bug in `plotRefToTarget` method = "TPS" where the wrong options from gridPars were being passed

------

# CHANGES IN GEOMORPH VERSION 2.1.6 (Patch Release)

### NEW FEATURES
* New globalIntegration` function for evaluating integration vs. self-similarity of shape variation 
* Coordinates returned by digitize2d` are now unscaled, and SCALE= returns the scale

### OTHER CHANGES
* Phylogenetic simulation procedure in `compare.evol.rates` generalized to use a single evolutionary rate matrix
* `pairwiseD.test` and `pairwise.slope`.test now defunct
* `define.sliders.2d` and `define.sliders.3d` now defunct
* `read.morphologika` can read files with missing data
* Ability to not show plot added to `physignal` and `compare.evol.rates`

### BUG FIXES
* Corrected error in defining starting point in `digit.curves`

------

# CHANGES IN GEOMORPH VERSION 2.1.5 (Patch Release)


### NEW FEATURES
* ability to include pre-digitized landmarks added to `build.template` and `digitsurface`
* new `gridPar` is a new function to customize plots of `plotRefToTarget`
* new `digit.curves` is a new function to calculate equidistant semilandmarks along 2D and 3D curves
* `define.sliders` is new interactive function for defining sliding semilandmarks for 2D and 3D curves, plus an automatic mode when given a sequence of semilandmarks along a curve

### OTHER CHANGES
* `pairwiseD.test` and `pairwise.slope.test` deprecated 
* 'read' functions now allow both tab and space delimited files
* `define.sliders.2d` and `define.sliders.3d` deprecated (replaced by define.sliders)

### BUG FIXES
* Corrected an error in `plotAllometry` where verbose=T did not return

------

# CHANGES IN GEOMORPH VERSION 2.1.4 (Patch Release)


### NEW FEATURES

### OTHER CHANGES
* `warpRefMesh` generalized - now takes a mesh3d object (i.e. made from `read.ply`) rather than calling `read.ply` directly
* `read.morphologika` now reads [groups] option and adds these data to the $labels matrix
* `plotOutliers` now has option groups to plot outliers by levels(groups) using group means
* `morphol.disparity` help file updated to correctly indicate that group shape residuals, rather than shape values, themselves, are randomized in the permutation procedure
* Internal changes to support functions for compatibility with R 3.1.3 
* Generalized plot inputs in `gpagen`

### BUG FIXES
* Corrected error `readland.tps`
* Corrected errors `trajectory.analysis`
* Corrected an issue with `gpagen` that flipped principal axes
* Fixed error in `read.morphologika` with reading [wireframe] in some morphologika files

------

# CHANGES IN GEOMORPH VERSION 2.1.3 (Patch Release)


### NEW FEATURES
* new `plotOutliers` function to identify potential outliers
* new `define.links` function for enhanced plotting of shapes

### OTHER CHANGES
* Additional input options added to `pairwiseD.test` and `pairwise.slope.test` 
* Ability to accommodate singular phylogenetic covariance matrices in: `physignal`, `compare.evol.rates`, `procD.pgls` and `phylo.pls`
* Enhanced digitizing capability in: `build.template`, `define.modules`, `define.sliders.3d`, `digit.fixed`, `digitsurface`, and `editTemplate` )
* `plotAllometry` input can be 2D matrix or 3D array
* `read.ply` reads normals for enhanced downstream digitizing from ply files
* `readland.tps`  reads curves from tps files and convert them to landmarks (semilandmarks)
* `plotTangentSpace` has enhanced plotting flexibility with labels and colors

### BUG FIXES
* Corrected error printing output of ANOVA table of `bilat.symmetry`
* Removed redundant permutation loop in `phylo.pls` 

------

# CHANGES IN GEOMORPH VERSION 2.1.2 (Patch Release)

### NEW FEATURES
* New function `advanced.procD.lm` for statistically comparing two or explanatory models
* Added warping of outline in `plotRefToTarget`
* New function to warp a specimen outline to the reference: `warpRefOutline`
* Added ShowPlot option to `two.b.pls
* Added RRPP option to `procD.pgls`
* New dataset pupfish

### OTHER CHANGES
* Enhanced underlying code in `procD.lm`, `procD.pgls` `pairwiseD.test` and `pairwise.slope.test`
* Removed k>3 restriction in `arrayspecs`
* Added is.numeric check to `phylo.pls`
* Removed arrows from plots when groups included in `plotAllometry`

### BUG FIXES
* Fixed output of `procD.lm` when verbose=TRUE
* Fixed reflections of aligned coordinate axes when PrinAxes=TRUE in `gpagen`

------

# CHANGES IN GEOMORPH VERSION 2.1.1 (Patch Release)

### NEW FEATURES
* Specimens rotated to their principal axes in `gpagen` with option to disable

### OTHER CHANGES
* Underlying code for ancestral state estimation changed to use `fastAnc` (phytools)

### BUG FIXES
* Fixed concatenated SSCP matrix issue in `procD.lm`, `pairwise.D.test` and `pairwise.slope.test`
* Corrected issue reading specimen names in `readland.tps`
* Corrected color options for plotting groups in `plotTangentSpace`
* Corrected test for slope:group interaction in `plotAllometry`

------

# CHANGES IN GEOMORPH VERSION 2.1 (Minor Release)

### NEW FEATURES
* `procD.pgls` added to assess high-dimensional ANOVA and regression models in a phylogenetic context
* `pairwise.slope.test` added to compare slopes of regression lines
* Residual randomization options added to `procD.lm` and `pairwise.d.test`
* Enhanced capabilities of `digitize2d`. Function now reads multiple images and outputs TPS file, can be used with missing data, and digitizing session can be restarted where previous session stopped

### OTHER CHANGES
* Ability to plot specimen labels added to `two.b.pls`, `morphol.integr`, and `phylo.pls`
* Slight ANOVA table output adjustment in `bilat.symmetry`
* Vector of labels can be added for plotting in `plotAllometry` and `plotTangentSpace`
* Labels for ancestral states added to `plotGMPhyloMorphoSpace`

### BUG FIXES
* Fixed scale issue in `digitize2d`

------

# CHANGES IN GEOMORPH VERSION 2.0.1 (Patch Release)

### NEW FEATURES

### OTHER CHANGES

### BUG FIXES
* Small change to C-code for `gpagen`

------

# CHANGES IN GEOMORPH VERSION 2.0 (Major Release)

### NEW FEATURES
* New function `phylo.pls` for assessing the multivariate association between two blocks of variables in a phylogenetic context
* New function `two.b.pls` for assessing the multivariate association between two blocks of variables
* New function `morphol.disparity` to compare Procrustes variance disparity among groups
* F-ratios and R-squared values added to output of `ProcD.lm`
* 3D Visualizations now include "surface" option to view shape deformation as warped mesh3d surfaces in the following: `plotRefToTarget`, `plotTangentSpace`, `plotAllometry`, and `bilat.symmetry`
* New I/O functions: `warpRefMesh` to create a mesh3d surface that represents the mean shape, `findMeanSpec` to assist in choosing a template ply file for use with `warpRefMesh` that identifies specimen closest to the mean shape, and `defline.modules` to interactively assign landmarks to modular partitions [currently 2D only]
* Generalized data input to allow 3D array or 2D matrix of data added to the following analysis functions: `compare.evol.rates`, `phylo.pls`, `morphol.integr`, `two.b.pls`, `physignal`, and `plotGMPhyloMorphoSpace`
* Ability to input univariate data added to `compare.evol.rates` and `physignal`
* Verbose output = T/F added to the following functions: `bilat.symmetry`, `phylo.pls`, `two.b.pls`, `morphol.integr`, `plotAllometry`, `plotTangentSpace`, `physignal`

### OTHER CHANGES
* Added calculation of pairwise Pvalues, and the option to assess a single group in  function `compare.evol.rates` 
* Additional graphical output added to `morphol.integr`
* Missing data handling altered (now NA is used)
* byLand option in `arrayspecs` has been removed
* Residual shape component (RSC) plot added and scores returned for `plotAllometry` (method = "CAC")
* Procrustes ANOVA added to `plotAllometry` output
* Centering = T/F option added to following 3D digitizing functions: `build.Template`, `digit.fixed`, `digitsurface`, and `plotSpec`
* `read.vrml` now defunct
* Major re-organization of underlying R code structure and format
* Optimized code to improve speed and performance in following functions: `arrayspecs`, `readland.tps`, `readland.nts`, `readmulti.nts`, `two.d.array`, `plotTangentSpace`, `trajectory.analysis`, `bilat.symmetry`, `gpagen` 
* Simplified plotting options in `bilat.symmetry`

### BUG FIXES
* Corrected small coding error in `digitize2D` and updated flexibility of the function

------

# CHANGES IN GEOMORPH VERSION 1.1-6 (Patch Release)

### NEW FEATURES

### OTHER CHANGES
* Minor I/O enhancements in `readmulti.nts`
* Minor I/O enhancements in `define.sliders.3d` to allow sliders to be in any order
* Simplified `pPsup` (original code from J. Claude) to not include size re-scaling by beta (underlying function used in `trajectory.analysis` only)
* Added name.check for groups in `compare.evol.rates`

### BUG FIXES

------

# CHANGES IN GEOMORPH VERSION 1.1-5 (Patch Release)

### NEW FEATURES
* New function `compare.evol.rates` for comparing multivariate evolutionary rates on phylogenies
* `define.sliders.2d` and `define.sliders.3d` replace `curves2d` and `digit.curves` 
* Option allowing specimens to be colored by group added to `plotTangentSpace` and `PlotAllometry` 

### OTHER CHANGES
* Simplified options in `morphol.integr`
* `curves2d` and `digit.curves` deprecated 

### BUG FIXES
* Corrected parameter estimates when groups specified for Regression Score option in `plotAllometry`

------

# CHANGES IN GEOMORPH VERSION 1.1-4 (Patch Release)

### NEW FEATURES
* Enhanced plotting of ply files in `read.ply`
* `digitsurface`, `buildtemplate`, `plotspec`, `digitfixed`, and `digitcurves` now support ply file input

### OTHER CHANGES
* Minor changes to plot window options in plotting functions
* Change to magnification factor usage in `plotRefToTartget`

### BUG FIXES

------

# CHANGES IN GEOMORPH VERSION 1.1-3 (Patch Release)

### NEW FEATURES

### OTHER CHANGES
* Improved NAMESPACE file and package usage
* Generalized `read.vrml` code for additional file formats

### BUG FIXES

------

# CHANGES IN GEOMORPH VERSION 1.1-2 (Patch Release)

### NEW FEATURES
* `pairwiseD.test` function added
* `bilat.symmetry` output includes symmetric and asymmetric shape components

### OTHER CHANGES
* Adjusted plotting routines in `morphol.integr` to be compatible with new CRAN guidelines
* Adjusted plotting routines in `bilat.symmetry` to be compatible with new CRAN guidelines
* Alternative ancestral state reconstruction and tests for bifurcating tree implemented in `physignal`
* Alternative ancestral state reconstruction and tests for bifurcating tree implemented in `PlotGMPhyloMorphoSpace`

### BUG FIXES
* Corrected `readland.tps` to allow for non-numeric ID and reading a single specimen per file
* Corrected landmark plotting issue and added greater directory flexibility in `curves2D` 
* Added greater directory flexibility and fixed header output in `digitize2D` 
* Added greater flexibility in reading distinct file formats in `read.morphologika`
* Corrected angle calculations in `fixed.angle`
* Corrected plotting of deformations grids in `plotTangentSpace`
* Corrected ancestral state output in `physignal` and `PlotGMPhyloMorphoSpace`

------

# CHANGES IN GEOMORPH VERSION 1.1-1 (Patch Release)

### NEW FEATURES

### OTHER CHANGES
* Removed dependency in `physignal` on `getAncStates` from `geiger`, which is no longer supported
* Removed dependency in `plotGMPhyloMorphoSpace` on `getAncStates` from `geiger`, which is no longer supported

### BUG FIXES

------

# CHANGES IN GEOMORPH VERSION 1.1-0 (Patch Release)


### NEW FEATURES
* `bilat.symmetry` function added
* `writeland.tps` function added
* `fixed.angle` function added
* `compare.modular.partitions` generalized to allow 2 or more partitions
* `morphol.integr` generalized to allow 2 or more partitions
* `trajectory.analysis` re-written to accept formulas, allowing greater flexibility for motion analysis
* PLS scores added to output `morphol.integr` 
* Ancestral states added to output `physignal`
* Centroid size and allometry scores added as output in `plotAllometry`
* PC scores added to output of `plotTangentSpace`
* Option added to select PC axes for plot in `plotTangentSpace`
* Option added to include specimen numbers to `plotTangentSpace`

### OTHER CHANGES
* `read.morphologika` accepts greater variety of input file formats

### BUG FIXES
* `buildtemplate` positional error in plot between template and scan corrected
* `digit.curves` error with passing objects to internal function corrected
* `gpgen` occasional reflection issue corrected

* Added a `NEWS.md` file to track changes to the package.
