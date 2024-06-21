### -----------------------------------------------------------------------
###
### All parallel processing options must be run outside of CRAN check
###
### All tests that require dependencies should be run outside of CRAN check
###
### -----------------------------------------------------------------------

##NOT EXAMINED
  # all interactive functions (define.links, etc.)
  # make.ggplot
  # read and write data functions

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
  succeed(summary(mosquito.sym))
})

test_that("bilat.symmetry2.works", {
  data(mosquito)
  Y.gpa <- gpagen(mosquito$wingshape)
  succeed(mosquito.sym2 <- bilat.symmetry(A = Y.gpa, ind = mosquito$ind, 
    side = mosquito$side, replicate = mosquito$replicate, 
    object.sym = FALSE, RRPP = TRUE, iter = 3))
  succeed(summary(mosquito.sym2))
})

test_that("bilat.symmetry3.works", {
  data(lizards)
  gdf <- geomorph.data.frame(shape = lizards$coords, 
    ind = lizards$ind, 
    replicate = lizards$rep)
  succeed(liz.sym <- bilat.symmetry(A = shape, ind = ind, rep = rep, 
    object.sym = TRUE, land.pairs = lizards$lm.pairs, 
    iter = 3, data = gdf, RRPP = TRUE))
  succeed(summary(liz.sym))
})

test_that("bilat.symmetry4.works", {
  data(scallops)
  gdf <- geomorph.data.frame(shape = scallops$coorddata, ind = scallops$ind)
  succeed(scallop.sym <- bilat.symmetry(A = shape, ind = ind, 
    object.sym = TRUE, curves= scallops$curvslide, surfaces = scallops$surfslide,
    land.pairs=scallops$land.pairs, iter = 3, data = gdf, RRPP = TRUE))
  succeed(summary(scallop.sym))
})

### combine.subsets --------------------------------------------------------------

test_that("combine.subsets1.works", {
  data(larvalMorph) 
  head.gpa <- gpagen(larvalMorph$headcoords, curves = larvalMorph$head.sliders)
  tail.gpa <- gpagen(larvalMorph$tailcoords, curves = larvalMorph$tail.sliders)
  succeed(all.lm <- combine.subsets(head = larvalMorph$headcoords,
    tail = larvalMorph$tailcoords, gpa = FALSE, CS.sets = NULL))
  succeed(plotAllSpecimens(all.lm$coords))
})

test_that("combine.subsets2.works", {
  data(larvalMorph) 
  head.gpa <- gpagen(larvalMorph$headcoords, curves = larvalMorph$head.sliders)
  tail.gpa <- gpagen(larvalMorph$tailcoords, curves = larvalMorph$tail.sliders)
  succeed(comb.lm <- combine.subsets(head = head.gpa, tail = tail.gpa, gpa = TRUE))
  succeed(comb.lm)
})

test_that("combine.subsets3.works", {
  data(larvalMorph) 
  head.gpa <- gpagen(larvalMorph$headcoords, curves = larvalMorph$head.sliders)
  tail.gpa <- gpagen(larvalMorph$tailcoords, curves = larvalMorph$tail.sliders)
  succeed(comb.lm <- combine.subsets(head = head.gpa$coords, tail = tail.gpa$coords, 
          gpa = FALSE, CS.sets = NULL))
  succeed(summary(comb.lm))
})

test_that("combine.subsets4.works", {
  data(larvalMorph) 
  head.gpa <- gpagen(larvalMorph$headcoords, curves = larvalMorph$head.sliders)
  tail.gpa <- gpagen(larvalMorph$tailcoords, curves = larvalMorph$tail.sliders)
  succeed(comb.lm <- combine.subsets(head = head.gpa, tail = tail.gpa, 
          gpa = TRUE, norm.CS = TRUE))
  succeed(summary(comb.lm))
})

test_that("combine.subsets5.works", {
  data(larvalMorph) 
  head.gpa <- gpagen(larvalMorph$headcoords, curves = larvalMorph$head.sliders)
  tail.gpa <- gpagen(larvalMorph$tailcoords, curves = larvalMorph$tail.sliders)
  succeed(comb.lm <- combine.subsets(head = head.gpa, 
         tail = tail.gpa, gpa = TRUE, norm.CS = FALSE, weights = c(0.3, 0.7)))
  succeed(summary(comb.lm))
})

### Compare.CR --------------------------------------------------------------

test_that("compareCR1.works", {
  data(pupfish) 
  Y.gpa<-gpagen(pupfish$coords)   
  land.gps<-rep('a',56); land.gps[39:48]<-'b'
  group <- factor(paste(pupfish$Pop, pupfish$Sex, sep = "."))
  coords.gp <- coords.subset(Y.gpa$coords, group)
  succeed(modul.tests <- Map(function(x) modularity.test(x, land.gps,iter=3, 
           print.progress = FALSE), coords.gp)) 
  succeed(group.Z <- compare.CR(modul.tests, CR.null = FALSE))
  succeed(summary(group.Z))
})

test_that("compareCR2.works", {
  data(pupfish) 
  Y.gpa<-gpagen(pupfish$coords)   
  land.gps<-rep('a',56); land.gps[39:48]<-'b'
  group <- factor(paste(pupfish$Pop, pupfish$Sex, sep = "."))
  coords.gp <- coords.subset(Y.gpa$coords, group)
  modul.tests <- Map(function(x) modularity.test(x, land.gps,iter=3, 
            print.progress = FALSE), coords.gp)
  land.gps3 <- rep('a',56); land.gps3[39:48]<-'b'
  land.gps3[c(6:9,28:38)] <- 'c' 
  land.gps4 <- rep('a',56); land.gps4[39:48]<-'b'
  land.gps4[c(6:9,28:38)] <- 'c'; land.gps4[c(10,49:56)] <- 'd'  
  m3.test <- modularity.test(coords.gp$Marsh.F,land.gps3, iter = 3, print.progress = FALSE)
  m4.test <- modularity.test(coords.gp$Marsh.F,land.gps4, iter = 3, print.progress = FALSE)
  succeed(model.Z <- compare.CR(modul.tests$Marsh.F,m3.test,m4.test, CR.null = TRUE))
  succeed(summary(model.Z))
})

### Compare.evol.rates --------------------------------------------------------------

test_that("compare.evol.rates1.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land)    
  gp.end <- factor(c(0,0,1,0,0,1,1,0,0))  
  names(gp.end) <- plethspecies$phy$tip
  succeed(ER<-compare.evol.rates(A = Y.gpa$coords, phy = plethspecies$phy,
          method = "simulation", iter = 3, gp = gp.end))
  succeed(summary(ER))
})

### Compare.multi.evol.rates --------------------------------------------------------------

test_that("compare.evol.rates1.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land)   
  land.gp <- c("A","A","A","A","A","B","B","B","B","B","B")
  succeed(EMR <- compare.multi.evol.rates(A = Y.gpa$coords, gp = land.gp, 
           Subset = TRUE, phy = plethspecies$phy))
  succeed(summary(EMR))
})

### Compare.physignal --------------------------------------------------------------

test_that("compare.physignal1.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land)
  jaw <- 1:5
  cranium <- 6:11
  PS.jaw <- physignal.z(A = Y.gpa$coords[jaw,,], phy = plethspecies$phy, 
                        lambda = "front", PAC.no = 7, iter=3)
  PS.cranium <- physignal.z(A = Y.gpa$coords[cranium,,], phy = plethspecies$phy, 
                            lambda = "front", PAC.no = 7, iter=3)
  PS.list <-list(PS.jaw, PS.cranium)
  names(PS.list) <- c("jaw", "cranium")
  
  succeed(PS.Z <- compare.physignal.z(PS.list))
  succeed(summary(PS.Z))
})

### Compare.pls --------------------------------------------------------------

test_that("compare.pls1.works", {
  data(pupfish) 
  group <- factor(paste(pupfish$Pop, pupfish$Sex, sep = "."))
  tail.LM <- c(1:3, 5:9, 18:38)
  head.LM <- (1:56)[-tail.LM]
  tail.coords <- pupfish$coords[tail.LM,,]
  head.coords <- pupfish$coords[head.LM,,]
  tail.coords.gp <- coords.subset(tail.coords, group)
  head.coords.gp <- coords.subset(head.coords, group)
  integ.tests <- Map(function(x,y) integration.test(x, y, iter=3, 
       print.progress = FALSE), head.coords.gp, tail.coords.gp)
  succeed(group.Z <- compare.pls(integ.tests))
  succeed(summary(group.Z))
})

### Compare.ZVrel --------------------------------------------------------------

test_that("compare.ZVrel1.works", {
  data("plethodon")
  Y.gpa <- gpagen(plethodon$land)
  coords.gp <- coords.subset(Y.gpa$coords, plethodon$species)
  Vrel.gp <- Map(function(x) integration.Vrel(x), coords.gp) 
  succeed(out <- compare.ZVrel(Vrel.gp$Jord, Vrel.gp$Teyah))
  succeed(summary(out))
})

### Coords.subset --------------------------------------------------------------

test_that("coords.subset1.works", {
  data(pupfish) 
  group <- factor(paste(pupfish$Pop, pupfish$Sex))
  succeed(new.coords <- coords.subset(A = pupfish$coords, group = group))
})

### estimate.missing --------------------------------------------------------------

test_that("estimate.missing1.works", {
  data(plethodon)
  plethland <- plethodon$land
     plethland[3,,2] <- plethland[8,,2] <- NA  #create missing landmarks
     plethland[3,,5] <- plethland[8,,5] <- plethland[9,,5] <- NA  
     plethland[3,,10] <- NA  
  succeed(estimate.missing(plethland, method = "TPS"))
})

test_that("estimate.missing2.works", {
  data(plethodon)
  plethland <- plethodon$land
  plethland[3,,2] <- plethland[8,,2] <- NA  #create missing landmarks
  plethland[3,,5] <- plethland[8,,5] <- plethland[9,,5] <- NA  
  plethland[3,,10] <- NA  
  succeed(estimate.missing(plethland, method = "Reg"))
})

### fixed.angle --------------------------------------------------------------

test_that("fixed.angle1.works", {
  data(plethspecies) 
  newLM1 <- fixed.angle(plethspecies$land,
    art.pt = 1, angle.pts.1 = 5, 
    angle.pts.2 = 6, rot.pts = c(2,3,4,5))
  succeed(Y.gpa1 <- gpagen(newLM1))
  succeed(plot(Y.gpa1, mean = FALSE))
})

test_that("fixed.angle2.works", {
  data(plethspecies) 
  newLM2 <- fixed.angle(plethspecies$land, art.pt = 1, 
    angle.pts.1 = c(1, 6:11), 
    angle.pts.2 = 2:5, 
    rot.pts = NULL, angle = 20, degrees = TRUE) 
  succeed(Y.gpa2 <- gpagen(newLM2))
  succeed(plot(Y.gpa2, mean = FALSE))
})

### geomorph.data.frame --------------------------------------------------------------

test_that("fixed.angle1.works", {
  data(plethodon) 
  Y.gpa <- gpagen(plethodon$land, PrinAxes = FALSE)
  succeed(gdf <- geomorph.data.frame(Y.gpa))
  succeed(attributes(gdf))
})

test_that("fixed.angle2.works", {
  data(plethodon) 
  Y.gpa <- gpagen(plethodon$land, PrinAxes = FALSE)
  succeed(gdf <- geomorph.data.frame(Y.gpa, species = plethodon$species, 
      site = plethodon$site))
  succeed(attributes(gdf))
})

### global.integration --------------------------------------------------------------

test_that("global.integration1.works", {
  data(plethodon) 
  Y.gpa <- gpagen(plethodon$land)
  succeed(res <- globalIntegration(Y.gpa$coords))
  succeed(summary(res))
})

### gm.measurement.error --------------------------------------------------------------

test_that("gm.measurement.error1.works", {
  data(fishy)
  fishy$coordsarray <- arrayspecs(fishy$coords, p = 11, k = 2)
  rep1 <- matrix(fishy$coords[1,], 11, 2, byrow = TRUE)
  rep2 <- matrix(fishy$coords[61,], 11, 2, byrow = TRUE)
  succeed(ME1 <- gm.measurement.error(coords = "coordsarray",
    subjects = "subj", replicates = "reps", data = fishy, iter = 3))
  succeed(anova(ME1))
  succeed(ICCstats(ME1, subjects = "Subjects", with_in = "Systematic ME"))
  succeed(plot(ME1))
})

test_that("gm.measurement.error2.works", {
  data(fishy)
  fishy$coordsarray <- arrayspecs(fishy$coords, p = 11, k = 2)
  rep1 <- matrix(fishy$coords[1,], 11, 2, byrow = TRUE)
  rep2 <- matrix(fishy$coords[61,], 11, 2, byrow = TRUE)
  succeed(ME2 <- gm.measurement.error(coords = "coordsarray", subjects = "subj", 
    replicates = "reps", groups = "groups", data = fishy, iter = 3))
  succeed(anova(ME2))
  succeed(ICCstats(ME2, subjects = "Subjects", 
      with_in = "Systematic ME", groups = "groups"))
  succeed(P <- plot(ME2))
  succeed(focusMEonSubjects(P, subjects = 18:20, shadow = TRUE))
  succeed(int.var <- interSubVar(ME2, type = "var"))
  succeed(plot(int.var))
})

### gm.prcomp --------------------------------------------------------------

test_that("gm.prcomp1.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land)
  succeed(PCA <- gm.prcomp(Y.gpa$coords))
  succeed(summary(PCA))
  succeed(plot(PCA, main = "PCA"))
  succeed(plot(PCA, main = "PCA", flip = 1))
  succeed(plot(PCA, main = "PCA", axis1 = 3, axis2 = 4))
})

test_that("gm.prcomp2.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land)
  succeed(PCA.w.phylo <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy))
  succeed(summary(PCA.w.phylo))
  succeed(plot(PCA.w.phylo, phylo = TRUE, main = "PCA.w.phylo"))
})

test_that("gm.prcomp3.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land)
  succeed(phylo.PCA <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, GLS = TRUE))
  succeed(summary(phylo.PCA))
  succeed(plot(phylo.PCA, phylo = TRUE, main = "phylo PCA"))
})

test_that("gm.prcomp4.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land)
  succeed(phylo.tPCA <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, 
    GLS = TRUE, transform = TRUE))
  succeed(summary(phylo.tPCA))
  succeed(plot(phylo.tPCA, phylo = TRUE, main = "phylo PCA"))
})

test_that("gm.prcomp5.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land)
  succeed(PaCA.ols <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, 
      align.to.phy = TRUE))
  succeed(summary(PaCA.ols))
  succeed(plot(PaCA.ols, phylo = TRUE, main = "PaCA using OLS"))
})

test_that("gm.prcomp6.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land)
  succeed(PaCA.gls <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, 
    align.to.phy = TRUE, GLS = TRUE))
  succeed(summary(PaCA.gls))
  succeed(plot(PaCA.gls, phylo = TRUE, main = "PaCA using GLS"))
})

test_that("gm.prcomp7.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land)
  succeed(PaCA.gls <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, 
    align.to.phy = TRUE, GLS = TRUE, transform = TRUE))
  succeed(summary(PaCA.gls))
  succeed(plot(PaCA.gls, phylo = TRUE, 
    main = "PaCA using GLS and transformed projection"))
})

test_that("gm.prcomp8.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land)
  PaCA.ols <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, 
    align.to.phy = TRUE)
  gps <- as.factor(c(rep("gp1", 5), rep("gp2", 4)))
  succeed(plot(PaCA.ols, pch=22, cex = 1.5, bg = gps, phylo = TRUE)) 
  succeed(text(par()$usr[1], 0.1*par()$usr[3], labels = "PC1 - 45.64%", 
    pos = 4, font = 2))
  succeed(text(0, 0.95*par()$usr[4], labels = "PC2 - 18.80%", pos = 4, font = 2))
  succeed(legend("topleft", pch=22, pt.bg = unique(gps), legend = levels(gps)))
})

### gpagen --------------------------------------------------------------

test_that("gpagen1.works", {
  data(plethodon) 
  succeed(Y.gpa <- gpagen(plethodon$land, PrinAxes = FALSE))
  succeed(summary(Y.gpa))
  succeed(plot(Y.gpa))
})

test_that("gpagen2.works", {
  data(hummingbirds)
  succeed(Y.gpa <- gpagen(hummingbirds$land, curves = hummingbirds$curvepts,
      ProcD = FALSE))   
  succeed(summary(Y.gpa))
  succeed(plot(Y.gpa))
})

test_that("gpagen3.works", {
  data(hummingbirds)
  succeed(Y.gpa <- gpagen(hummingbirds$land, curves = hummingbirds$curvepts,
                          ProcD = TRUE))   
  succeed(summary(Y.gpa))
  succeed(plot(Y.gpa))
})

test_that("gpagen4.works", {
  data(scallops)
  succeed(Y.gpa <- gpagen(A = scallops$coorddata, curves = scallops$curvslide, 
    surfaces = scallops$surfslide))
  succeed(summary(Y.gpa))
})

### integration.test --------------------------------------------------------------

test_that("integration.test1.works", {
  data(plethodon) 
  Y.gpa <- gpagen(plethodon$land)
  land.gps <- c("A","A","A","A","A","B","B","B","B","B","B","B") 
  succeed(IT <- integration.test(Y.gpa$coords, partition.gp = land.gps))
  succeed(summary(IT))
  succeed(P <- plot(IT))
  succeed(IT$left.pls.vectors)
})

### integration.Vrel --------------------------------------------------------------

test_that("integration.Vrel1.works", {
  data(plethodon) 
  Y.gpa <- gpagen(plethodon$land)     
  succeed(res <- integration.Vrel(Y.gpa$coords))
  succeed(print(res))
})


### interlmkdist --------------------------------------------------------------

test_that("interlmkdist1.works", {
  data(plethodon) 
  lmks <- data.frame(eyeW = c(8,9), headL = c(6,12), mouthL = c(4,2), 
  row.names = c("start", "end")) 
  A <- plethodon$land
  succeed(lineardists <- interlmkdist(A, lmks))
  succeed(lineardists)
})

### modularity.test --------------------------------------------------------------

test_that("modularity.test1.works", {
  data(pupfish) 
  Y.gpa <- gpagen(pupfish$coords, print.progress = FALSE)
  land.gps <- rep('a',56); land.gps[39:48] <- 'b'
  succeed(MT <- modularity.test(Y.gpa$coords, land.gps, CI = FALSE, iter = 3))
  succeed(summary(MT))
  succeed(plot(MT))
})

### morphol.disparity --------------------------------------------------------------

test_that("morphol.disparity1.works", {
  data(plethodon)
  Y.gpa <- gpagen(plethodon$land, print.progress = FALSE)
  gdf <- geomorph.data.frame(Y.gpa, species = plethodon$species, 
    site = plethodon$site)
  succeed(morphol.disparity(coords ~ 1, groups = NULL, data = gdf, 
    print.progress = FALSE, iter = 3))
  succeed(morphol.disparity(coords ~ Csize, groups= NULL, data = gdf, 
    print.progress = FALSE, iter = 3))
  succeed(morphol.disparity(coords ~ 1, groups= ~ species * site, data = gdf, 
    print.progress = FALSE, iter = 3))
  succeed(morphol.disparity(coords ~ 1, groups= ~ species * site, partial = TRUE, 
    data = gdf, print.progress = FALSE, iter = 3))
  succeed(morphol.disparity(coords ~ species * site, groups= ~species * site, 
    data = gdf, print.progress = FALSE, iter = 3))
  succeed(morphol.disparity(coords ~ Csize + species * site, groups= ~ species, 
    data = gdf, print.progress = FALSE, iter = 3))
})

test_that("morphol.disparity2.works", {
  data(plethodon)
  Y.gpa <- gpagen(plethodon$land, print.progress = FALSE)
  gdf <- geomorph.data.frame(Y.gpa, species = plethodon$species, 
     site = plethodon$site)
  succeed(MD <- morphol.disparity(coords ~ Csize + species * site, groups= ~ species, 
    data = gdf, print.progress = FALSE, iter = 3))
  succeed(MD$Procrustes.var)
})

test_that("morphol.disparity3.works", {
  data(plethspecies)
  Y.gpa <- gpagen(plethspecies$land)
  gp.end <- factor(c(0,0,1,0,0,1,1,0,0))
  names(gp.end) <- plethspecies$phy$tip
  gdf <- geomorph.data.frame(Y.gpa, phy = plethspecies$phy, 
    gp.end = gp.end)
  pleth.ols <- procD.lm(coords ~ Csize + gp.end, 
    data = gdf, iter = 3)
  pleth.pgls <- procD.pgls(coords ~ Csize + gp.end, phy = phy, 
    data = gdf, iter = 3)
  succeed(morphol.disparity(f1 = pleth.ols, groups = ~ gp.end, data = gdf, 
    print.progress = FALSE))
  succeed(morphol.disparity(f1 = pleth.pgls, groups = ~ gp.end, 
    transform = FALSE, data = gdf, print.progress = FALSE))
  succeed(morphol.disparity(f1 = pleth.pgls, groups = ~ gp.end,
    transform = TRUE, data = gdf, print.progress = FALSE))
  succeed(PW <- pairwise(pleth.ols, groups = gp.end))
  succeed(summary(PW, test.type = 'var'))
  succeed(PW2 <- pairwise(pleth.pgls, groups = gp.end))
  succeed(summary(PW2, test.type = 'var'))
})

### phylo.integration --------------------------------------------------------------

test_that("phylo.integration1.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land)
  land.gps <- c("A","A","A","A","A","B","B","B","B","B","B") 
  succeed(IT <- phylo.integration(Y.gpa$coords, partition.gp = land.gps,
    phy = plethspecies$phy, iter = 3))
  succeed(summary(IT))
  succeed(P <- plot(IT))
})

### phylo.moduarity --------------------------------------------------------------

test_that("phylo.modularity1.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land)
  land.gps <- c("A","A","A","A","A","B","B","B","B","B","B") 
  succeed(MT <- phylo.modularity(Y.gpa$coords, partition.gp = land.gps, 
    phy = plethspecies$phy, CI = FALSE, iter = 3))
  succeed(summary(MT))
  succeed(plot(MT)) 
})