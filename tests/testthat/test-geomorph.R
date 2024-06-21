### -----------------------------------------------------------------------
###
### All parallel processing options must be run outside of CRAN check
###
### All tests that require dependencies should be run outside of CRAN check
###
### -----------------------------------------------------------------------

### arrayspecs --------------------------------------------------------------

test_that("arrayspecs.works", {
  x <- matrix(rnorm(18), nrow = 3)
  succeed(arrayspecs(x, 3, 2)) 
})

### bilat.symmetry --------------------------------------------------------------

test_that("bilat.symmetry1.works", {
  data(mosquito)
  gdf <- geomorph.data.frame(wingshape = mosquito$wingshape, 
         ind = mosquito$ind, side = mosquito$side,
         replicate = mosquito$replicate)
  succeed(mosquito.sym <- bilat.symmetry(A = wingshape, ind = ind, side = side,
         replicate = replicate, object.sym = FALSE, RRPP = TRUE, 
         iter = 3, data = gdf))
})

## THIS ONE DOESN'T WORK.  AND IMPLIES THE HELP EXAMPLE REQUIRES TWEAKING FOR PROPER NOTATION
test_that("bilat.symmetry2.works", {
  data(mosquito)
  Y.gpa <- gpagen(mosquito$wingshape)
  succeed(mosquito.sym2 <- bilat.symmetry(A = Y.gpa, ind = mosquito$ind, 
    side = mosquito$side, replicate = mosquito$replicate, 
    object.sym = FALSE, RRPP = TRUE, iter = 3))
})

test_that("bilat.symmetry3.works", {
  data(lizards)
  gdf <- geomorph.data.frame(shape = lizards$coords, 
    ind = lizards$ind, 
    replicate = lizards$rep)
  succeed(liz.sym <- bilat.symmetry(A = shape, ind = ind, rep = rep, 
    object.sym = TRUE, land.pairs = lizards$lm.pairs, 
    iter = 3, data = gdf, RRPP = TRUE))
})


