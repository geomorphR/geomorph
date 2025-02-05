### -----------------------------------------------------------------------
###
### All parallel processing options must be run outside of CRAN check
###
### All tests that require dependencies should be run outside of CRAN check
###
### -----------------------------------------------------------------------

##NOT EXAMINED
  # all interactive functions (define.links, etc.)
  # make.ggplot; no 3D plots (RGL)
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
         iter = 3, data = gdf, print.progress = F))
  succeed(summary(mosquito.sym))
})

test_that("bilat.symmetry2.works", {
  data(mosquito)
  Y.gpa <- gpagen(mosquito$wingshape, print.progress = F)
  succeed(mosquito.sym2 <- bilat.symmetry(A = Y.gpa, ind = mosquito$ind, 
    side = mosquito$side, replicate = mosquito$replicate, 
    object.sym = FALSE, RRPP = TRUE, iter = 3, print.progress = F))
  succeed(summary(mosquito.sym2))
})

test_that("bilat.symmetry3.works", {
  data(lizards)
  gdf <- geomorph.data.frame(shape = lizards$coords, 
    ind = lizards$ind, 
    replicate = lizards$rep)
  succeed(liz.sym <- bilat.symmetry(A = shape, ind = ind, rep = rep, 
    object.sym = TRUE, land.pairs = lizards$lm.pairs, 
    iter = 3, data = gdf, RRPP = TRUE, print.progress = F))
  succeed(summary(liz.sym))
})

test_that("bilat.symmetry4.works", {
  data(scallops)
  gdf <- geomorph.data.frame(shape = scallops$coorddata, ind = scallops$ind)
  succeed(scallop.sym <- bilat.symmetry(A = shape, ind = ind, 
    object.sym = TRUE, curves= scallops$curvslide, surfaces = scallops$surfslide,
    land.pairs=scallops$land.pairs, iter = 3, data = gdf, RRPP = TRUE, 
         print.progress = F))
  succeed(summary(scallop.sym))
})

### combine.subsets --------------------------------------------------------------

test_that("combine.subsets1.works", {
  data(larvalMorph) 
  head.gpa <- gpagen(larvalMorph$headcoords, curves = larvalMorph$head.sliders, 
                     print.progress = F)
  tail.gpa <- gpagen(larvalMorph$tailcoords, curves = larvalMorph$tail.sliders, 
                     print.progress = F)
  succeed(all.lm <- combine.subsets(head = larvalMorph$headcoords,
    tail = larvalMorph$tailcoords, gpa = FALSE, CS.sets = NULL))
  succeed(plotAllSpecimens(all.lm$coords))
})

test_that("combine.subsets2.works", {
  data(larvalMorph) 
  head.gpa <- gpagen(larvalMorph$headcoords, curves = larvalMorph$head.sliders, 
                     print.progress = F)
  tail.gpa <- gpagen(larvalMorph$tailcoords, curves = larvalMorph$tail.sliders, 
                     print.progress = F)
  succeed(comb.lm <- combine.subsets(head = head.gpa, tail = tail.gpa, gpa = TRUE))
  succeed(comb.lm)
})

test_that("combine.subsets3.works", {
  data(larvalMorph) 
  head.gpa <- gpagen(larvalMorph$headcoords, curves = larvalMorph$head.sliders, 
                     print.progress = F)
  tail.gpa <- gpagen(larvalMorph$tailcoords, curves = larvalMorph$tail.sliders, 
                     print.progress = F)
  succeed(comb.lm <- combine.subsets(head = head.gpa$coords, tail = tail.gpa$coords, 
          gpa = FALSE, CS.sets = NULL))
  succeed(summary(comb.lm))
})

test_that("combine.subsets4.works", {
  data(larvalMorph) 
  head.gpa <- gpagen(larvalMorph$headcoords, curves = larvalMorph$head.sliders, 
                     print.progress = F)
  tail.gpa <- gpagen(larvalMorph$tailcoords, curves = larvalMorph$tail.sliders, 
                     print.progress = F)
  succeed(comb.lm <- combine.subsets(head = head.gpa, tail = tail.gpa, 
          gpa = TRUE, norm.CS = TRUE))
  succeed(summary(comb.lm))
})

test_that("combine.subsets5.works", {
  data(larvalMorph) 
  head.gpa <- gpagen(larvalMorph$headcoords, curves = larvalMorph$head.sliders, 
                     print.progress = F)
  tail.gpa <- gpagen(larvalMorph$tailcoords, curves = larvalMorph$tail.sliders, 
                     print.progress = F)
  succeed(comb.lm <- combine.subsets(head = head.gpa, 
         tail = tail.gpa, gpa = TRUE, norm.CS = FALSE, weights = c(0.3, 0.7)))
  succeed(summary(comb.lm))
})

### Compare.CR --------------------------------------------------------------

test_that("compareCR1.works", {
  data(pupfish) 
  Y.gpa<-gpagen(pupfish$coords, print.progress = F)   
  land.gps<-rep('a',56); land.gps[39:48]<-'b'
  group <- factor(paste(pupfish$Pop, pupfish$Sex, sep = "."))
  coords.gp <- coords.subset(Y.gpa$coords, group)
  succeed(modul.tests <- Map(function(x) modularity.test(x, land.gps,iter=3, 
               opt.rot = FALSE, print.progress = FALSE), coords.gp)) 
  succeed(group.Z <- compare.CR(modul.tests, CR.null = FALSE))
  succeed(summary(group.Z))
})

test_that("compareCR2.works", {
  data(pupfish) 
  Y.gpa<-gpagen(pupfish$coords, print.progress = F)   
  land.gps<-rep('a',56); land.gps[39:48]<-'b'
  group <- factor(paste(pupfish$Pop, pupfish$Sex, sep = "."))
  coords.gp <- coords.subset(Y.gpa$coords, group)
  modul.tests <- Map(function(x) modularity.test(x, land.gps,iter=3, 
            opt.rot = FALSE, print.progress = FALSE), coords.gp)
  land.gps3 <- rep('a',56); land.gps3[39:48]<-'b'
  land.gps3[c(6:9,28:38)] <- 'c' 
  land.gps4 <- rep('a',56); land.gps4[39:48]<-'b'
  land.gps4[c(6:9,28:38)] <- 'c'; land.gps4[c(10,49:56)] <- 'd'  
  m3.test <- modularity.test(coords.gp$Marsh.F,land.gps3, iter = 3, opt.rot = FALSE, 
                             print.progress = FALSE)
  m4.test <- modularity.test(coords.gp$Marsh.F,land.gps4, iter = 3, opt.rot = FALSE, 
                             print.progress = FALSE)
  succeed(model.Z <- compare.CR(modul.tests$Marsh.F,m3.test,m4.test, CR.null = TRUE))
  succeed(summary(model.Z))
})

### Compare.evol.rates --------------------------------------------------------------

test_that("compare.evol.rates1.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land, print.progress = F)    
  gp.end <- factor(c(0,0,1,0,0,1,1,0,0))  
  names(gp.end) <- plethspecies$phy$tip
  succeed(ER<-compare.evol.rates(A = Y.gpa$coords, phy = plethspecies$phy,
          method = "simulation", iter = 3, gp = gp.end, print.progress = F))
  succeed(summary(ER))
})

### Compare.multi.evol.rates --------------------------------------------------------------

test_that("compare.evol.rates1.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land, print.progress = F)   
  land.gp <- c("A","A","A","A","A","B","B","B","B","B","B")
  succeed(EMR <- compare.multi.evol.rates(A = Y.gpa$coords, gp = land.gp, 
           Subset = TRUE, phy = plethspecies$phy, iter = 3, print.progress = F))
  succeed(summary(EMR))
})

### Compare.physignal --------------------------------------------------------------

test_that("compare.physignal1.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land, print.progress = F)
  jaw <- 1:5
  cranium <- 6:11
  PS.jaw <- physignal.z(A = Y.gpa$coords[jaw,,], phy = plethspecies$phy, 
                        lambda = "front", PAC.no = 7, iter=3, print.progress = F)
  PS.cranium <- physignal.z(A = Y.gpa$coords[cranium,,], phy = plethspecies$phy, 
                            lambda = "front", PAC.no = 7, iter=3, print.progress = F)
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
  Y.gpa <- gpagen(plethodon$land, print.progress = F)
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
  succeed(Y.gpa1 <- gpagen(newLM1, print.progress = F))
  succeed(plot(Y.gpa1, mean = FALSE))
})

test_that("fixed.angle2.works", {
  data(plethspecies) 
  newLM2 <- fixed.angle(plethspecies$land, art.pt = 1, 
    angle.pts.1 = c(1, 6:11), 
    angle.pts.2 = 2:5, 
    rot.pts = NULL, angle = 20, degrees = TRUE) 
  succeed(Y.gpa2 <- gpagen(newLM2, print.progress = F))
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
  Y.gpa <- gpagen(plethodon$land, PrinAxes = FALSE, print.progress = F)
  succeed(gdf <- geomorph.data.frame(Y.gpa, species = plethodon$species, 
      site = plethodon$site))
  succeed(attributes(gdf))
})

### global.integration --------------------------------------------------------------

test_that("global.integration1.works", {
  data(plethodon) 
  Y.gpa <- gpagen(plethodon$land, print.progress = F)
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
    subjects = "subj", replicates = "reps", data = fishy, turbo = FALSE, 
    iter = 3, print.progress = F))
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
    replicates = "reps", groups = "groups", data = fishy, turbo = FALSE, 
    iter = 3, print.progress = F))
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
  Y.gpa <- gpagen(plethspecies$land, print.progress = F)
  succeed(PCA <- gm.prcomp(Y.gpa$coords))
  succeed(summary(PCA))
  succeed(plot(PCA, main = "PCA"))
  succeed(plot(PCA, main = "PCA", flip = 1))
  succeed(plot(PCA, main = "PCA", axis1 = 3, axis2 = 4))
})

test_that("gm.prcomp2.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land, print.progress = F)
  succeed(PCA.w.phylo <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy))
  succeed(summary(PCA.w.phylo))
  succeed(plot(PCA.w.phylo, phylo = TRUE, main = "PCA.w.phylo"))
})

test_that("gm.prcomp3.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land, print.progress = F)
  succeed(phylo.PCA <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, GLS = TRUE))
  succeed(summary(phylo.PCA))
  succeed(plot(phylo.PCA, phylo = TRUE, main = "phylo PCA"))
})

test_that("gm.prcomp4.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land, print.progress = F)
  succeed(phylo.tPCA <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, 
    GLS = TRUE, transform = TRUE))
  succeed(summary(phylo.tPCA))
  succeed(plot(phylo.tPCA, phylo = TRUE, main = "phylo PCA"))
})

test_that("gm.prcomp5.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land, print.progress = F)
  succeed(PaCA.ols <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, 
      align.to.phy = TRUE))
  succeed(summary(PaCA.ols))
  succeed(plot(PaCA.ols, phylo = TRUE, main = "PaCA using OLS"))
})

test_that("gm.prcomp6.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land, print.progress = F)
  succeed(PaCA.gls <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, 
    align.to.phy = TRUE, GLS = TRUE))
  succeed(summary(PaCA.gls))
  succeed(plot(PaCA.gls, phylo = TRUE, main = "PaCA using GLS"))
})

test_that("gm.prcomp7.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land, print.progress = F)
  succeed(PaCA.gls <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, 
    align.to.phy = TRUE, GLS = TRUE, transform = TRUE))
  succeed(summary(PaCA.gls))
  succeed(plot(PaCA.gls, phylo = TRUE, 
    main = "PaCA using GLS and transformed projection"))
})

test_that("gm.prcomp8.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land, print.progress = F)
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
  succeed(Y.gpa <- gpagen(plethodon$land, PrinAxes = FALSE, print.progress = F))
  succeed(summary(Y.gpa))
  succeed(plot(Y.gpa))
})

test_that("gpagen2.works", {
  data(hummingbirds)
  succeed(Y.gpa <- gpagen(hummingbirds$land, curves = hummingbirds$curvepts,
      ProcD = FALSE, print.progress = F))   
  succeed(summary(Y.gpa))
  succeed(plot(Y.gpa))
})

test_that("gpagen3.works", {
  data(hummingbirds)
  succeed(Y.gpa <- gpagen(hummingbirds$land, curves = hummingbirds$curvepts,
                          ProcD = TRUE, print.progress = F))   
  succeed(summary(Y.gpa))
  succeed(plot(Y.gpa))
})

test_that("gpagen4.works", {
  data(scallops)
  succeed(Y.gpa <- gpagen(A = scallops$coorddata, curves = scallops$curvslide, 
    surfaces = scallops$surfslide, print.progress = F))
  succeed(summary(Y.gpa))
})

### integration.test --------------------------------------------------------------

test_that("integration.test1.works", {
  data(plethodon) 
  Y.gpa <- gpagen(plethodon$land, print.progress = F)
  land.gps <- c("A","A","A","A","A","B","B","B","B","B","B","B") 
  succeed(IT <- integration.test(Y.gpa$coords, partition.gp = land.gps, 
                                 iter = 3, print.progress = F))
  succeed(summary(IT))
  succeed(P <- plot(IT))
  succeed(IT$left.pls.vectors)
})

### integration.Vrel --------------------------------------------------------------

test_that("integration.Vrel1.works", {
  data(plethodon) 
  Y.gpa <- gpagen(plethodon$land, print.progress = F)     
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
  succeed(MT <- modularity.test(Y.gpa$coords, land.gps, opt.rot = FALSE, 
                                CI = FALSE, iter = 3, print.progress = F))
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
  Y.gpa <- gpagen(plethspecies$land, print.progress = F)
  gp.end <- factor(c(0,0,1,0,0,1,1,0,0))
  names(gp.end) <- plethspecies$phy$tip
  gdf <- geomorph.data.frame(Y.gpa, phy = plethspecies$phy, 
    gp.end = gp.end)
  pleth.ols <- procD.lm(coords ~ Csize + gp.end, 
    data = gdf, iter = 3, print.progress = F)
  pleth.pgls <- procD.pgls(coords ~ Csize + gp.end, phy = phy, 
    data = gdf, iter = 3, print.progress = F)
  succeed(morphol.disparity(f1 = pleth.ols, groups = ~ gp.end, data = gdf, 
    print.progress = FALSE, iter = 3))
  succeed(morphol.disparity(f1 = pleth.pgls, groups = ~ gp.end, 
    transform = FALSE, data = gdf, print.progress = FALSE, iter = 3))
  succeed(morphol.disparity(f1 = pleth.pgls, groups = ~ gp.end,
    transform = TRUE, data = gdf, print.progress = FALSE, iter = 3))
  succeed(PW <- pairwise(pleth.ols, groups = gp.end, print.progress = F))
  succeed(summary(PW, test.type = 'var'))
  succeed(PW2 <- pairwise(pleth.pgls, groups = gp.end, print.progress = F))
  succeed(summary(PW2, test.type = 'var'))
})

### phylo.integration --------------------------------------------------------------

test_that("phylo.integration1.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land, print.progress = F)
  land.gps <- c("A","A","A","A","A","B","B","B","B","B","B") 
  succeed(IT <- phylo.integration(Y.gpa$coords, partition.gp = land.gps,
    phy = plethspecies$phy, iter = 3, print.progress = F))
  succeed(summary(IT))
  succeed(P <- plot(IT))
})

### phylo.moduarity --------------------------------------------------------------

test_that("phylo.modularity1.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land, print.progress = F)
  land.gps <- c("A","A","A","A","A","B","B","B","B","B","B") 
  succeed(MT <- phylo.modularity(Y.gpa$coords, partition.gp = land.gps, 
    phy = plethspecies$phy, CI = FALSE, iter = 3, print.progress = F))
  succeed(summary(MT))
  succeed(plot(MT)) 
})

### physignal --------------------------------------------------------------

test_that("physignal1.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land, print.progress = F)
  succeed(PS.shape <- physignal(A = Y.gpa$coords, phy = plethspecies$phy,
              iter = 3, print.progress = F))
  succeed(summary(PS.shape))
  succeed((PS.shape))
  succeed(plot(PS.shape$PACA, phylo = TRUE))
  succeed(PS.shape$K.by.p) 
})

### physignal.z --------------------------------------------------------------

test_that("physignal.z1.works", {
  skip_on_cran()
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land, print.progress = F)   
  succeed(PS.shape <- physignal.z(A = Y.gpa$coords, phy = plethspecies$phy, 
       lambda = "front", PAC.no = 7, iter = 3, print.progress = F))
  succeed(summary(PS.shape))
  succeed(plot(PS.shape))
  succeed(plot(PS.shape$PACA, phylo = TRUE))
})

### plotAllometry --------------------------------------------------------------

test_that("plotAllemtry1.works", {
  data(plethodon) 
  Y.gpa <- gpagen(plethodon$land, print.progress = FALSE) 
  gdf <- geomorph.data.frame(Y.gpa, site = plethodon$site, 
    species = plethodon$species) 
  fit <- procD.lm(coords ~ log(Csize), data = gdf, 
    iter = 3, print.progress = FALSE)
  succeed(plotAllometry(fit, size = gdf$Csize, logsz = TRUE, 
    method = "PredLine", pch = 19))
  logSize <- log(gdf$Csize)
  succeed(plot(fit, type = "regression", reg.type = "PredLine", 
    predictor = logSize, pch = 19))
})

test_that("plotAllemtry2.works", {
  data(plethodon) 
  Y.gpa <- gpagen(plethodon$land, print.progress = FALSE) 
  gdf <- geomorph.data.frame(Y.gpa, site = plethodon$site, 
                             species = plethodon$species) 
  logSize <- log(gdf$Csize)
  fit <- procD.lm(coords ~ log(Csize), data = gdf, 
                  iter = 3, print.progress = FALSE)
  succeed(plotAllometry(fit, size = gdf$Csize, logsz = TRUE, 
    method = "RegScore", pch = 19))
  succeed(plot(fit, type = "regression", reg.type = "RegScore", 
    predictor = logSize, pch = 19))
})

test_that("plotAllemtry3.works", {
  data(plethodon) 
  Y.gpa <- gpagen(plethodon$land, print.progress = FALSE) 
  gdf <- geomorph.data.frame(Y.gpa, site = plethodon$site, 
                             species = plethodon$species) 
  logSize <- log(gdf$Csize)
  fit <- procD.lm(coords ~ log(Csize), data = gdf, 
                  iter = 3, print.progress = FALSE)
  succeed(plotAllometry(fit, size = gdf$Csize, logsz = TRUE, 
    method = "CAC", pch = 19))
  succeed(PLS <- two.b.pls(log(gdf$Csize), gdf$coords, iter = 3, 
                           print.progress = FALSE))
  succeed(plot(PLS))
})

test_that("plotAllemtry4.works", {
  data(plethodon) 
  Y.gpa <- gpagen(plethodon$land, print.progress = FALSE) 
  gdf <- geomorph.data.frame(Y.gpa, site = plethodon$site, 
                             species = plethodon$species) 
  logSize <- log(gdf$Csize)
  fit <- procD.lm(coords ~ Csize * species * site, data = gdf, 
    iter = 3, print.progress = FALSE)
  succeed(plotAllometry(fit, size = gdf$Csize, logsz = TRUE, method = "CAC", 
    pch = 19))
  succeed(plotAllometry(fit, size = gdf$Csize, logsz = TRUE, method = "PredLine", 
    pch = 19, col = as.numeric(interaction(gdf$species, gdf$site))))
  succeed(plotAllometry(fit, size = gdf$Csize, logsz = TRUE, method = "RegScore", 
    pch = 19, col = as.numeric(interaction(gdf$species, gdf$site))))
  succeed(pc.plot <- plotAllometry(fit, size = gdf$Csize, logsz = TRUE, 
    method = "size.shape", pch = 19, 
    col = as.numeric(interaction(gdf$species, gdf$site))))
  succeed(summary(pc.plot$size.shape.PCA))
})

test_that("plotAllemtry5.works", {
  data(plethodon) 
  Y.gpa <- gpagen(plethodon$land, print.progress = FALSE) 
  gdf <- geomorph.data.frame(Y.gpa, site = plethodon$site, 
                             species = plethodon$species) 
  fit3 <- procD.lm(coords ~ species, data = gdf, 
    iter = 3, print.progress = FALSE)
  succeed(plotAllometry(fit3, size = gdf$Csize, logsz = TRUE, method = "RegScore", 
    pch = 19, col = as.numeric(gdf$species)))
})

### plotAllSpecimens --------------------------------------------------------------

test_that("plotAllSpecimens1.works", {
  data(plethodon) 
  Y.gpa <- gpagen(plethodon$land, print.progress = F)
  succeed(plotAllSpecimens(Y.gpa$coords, links = plethodon$links))
})

### plotOutliers --------------------------------------------------------------

test_that("plotOutliers1.works", {
  data(plethodon)
  newland <- plethodon$land
  newland[c(1,8),,2] <- newland[c(8,1),,2]
  newland[c(3,11),,26] <- newland[c(11,3),,2]
  Y<- gpagen(newland, print.progress = F) 
  succeed(out <- plotOutliers(Y$coords))
  succeed(plotOutliers(Y$coords, inspect.outliers = TRUE))
  succeed(plotOutliers(Y$coords, groups = plethodon$species, 
    inspect.outliers = TRUE))
})

### plotRefToTarget --------------------------------------------------------------

test_that("plotRefToTarget1.works", {
  data(plethodon) 
  Y.gpa <- gpagen(plethodon$land, print.progress = F)
  ref <- mshape(Y.gpa$coords)
  succeed(plotRefToTarget(ref, Y.gpa$coords[,,39]))
  succeed(plotRefToTarget(ref, Y.gpa$coords[,,39], mag = 2, 
    outline = plethodon$outline)) 
})

test_that("plotRefToTarget2.works", {
  data(plethodon) 
  Y.gpa <- gpagen(plethodon$land, print.progress = F)
  ref <- mshape(Y.gpa$coords)
  succeed(plotRefToTarget(ref, Y.gpa$coords[,,39], method = "vector", mag = 3))
  succeed(plotRefToTarget(ref, Y.gpa$coords[,,39], method = "points", 
    outline = plethodon$outline))
  succeed(plotRefToTarget(ref, Y.gpa$coords[,,39], method = "vector", 
    outline = plethodon$outline, mag = 2.5))
  succeed(plotRefToTarget(ref, Y.gpa$coords[,,39], 
    gridPars = gridPar(pt.bg = "green", pt.size = 1),
    method = "vector", mag = 3))
})

### procD.lm --------------------------------------------------------------

test_that("procD.lm.example1.works", {
  data(plethodon) 
  Y.gpa <- gpagen(plethodon$land, print.progress = F)
  gdf <- geomorph.data.frame(Y.gpa, 
                             site = plethodon$site, 
                             species = plethodon$species) 
  succeed(fit1 <- procD.lm(coords ~ species * site, 
                           data = gdf, iter = 3, verbose = TRUE,
                           RRPP = FALSE, print.progress = FALSE))
  succeed(summary(fit1))
})

test_that("procD.lm.example2.works", {
  data(plethodon) 
  Y.gpa <- gpagen(plethodon$land, print.progress = F)
  gdf <- geomorph.data.frame(Y.gpa, 
                             site = plethodon$site, 
                             species = plethodon$species) 
  succeed(fit2 <- procD.lm(coords ~ species * site, 
                           data = gdf, iter = 3, verbose = TRUE,
                           RRPP = TRUE, print.progress = FALSE))
  succeed(summary(fit2))
  succeed(coef(fit2))
  succeed(coef(fit2, test = TRUE))
  succeed(anova(fit2)) 
  succeed(anova(fit2, effect.type = "Rsq"))
  succeed(anova(fit2, error = c("species:site", "species:site", "Residuals")))
})

test_that("procD.lm.example3.works", {
  data(plethodon) 
  Y.gpa <- gpagen(plethodon$land, print.progress = F)
  gdf <- geomorph.data.frame(Y.gpa, 
                             site = plethodon$site, 
                             species = plethodon$species) 
  succeed(fit.null <- procD.lm(coords ~ log(Csize) + species + site, data = gdf, 
                               iter = 3, print.progress = FALSE, verbose = TRUE))
  succeed(fit.full <- procD.lm(coords ~ log(Csize) + species * site, data = gdf, 
                               iter = 3, print.progress = FALSE, verbose = TRUE))
  succeed(anova(fit.null, fit.full, print.progress = FALSE))
  succeed(gp <-  interaction(gdf$species, gdf$site))
  succeed(PW <- pairwise(fit.full, groups = gp, covariate = NULL))
  succeed(summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE))
  succeed(summary(PW, test.type = "dist", confidence = 0.95, stat.table = FALSE))
  succeed(summary(PW, test.type = "var", confidence = 0.95, stat.table = TRUE))
  succeed(summary(PW, test.type = "var", confidence = 0.95, stat.table = FALSE))
  succeed(morphol.disparity(fit.full, groups = gp, iter = 3))
})

test_that("procD.lm.example4.works", {  
  data(ratland)
  rat.gpa<-gpagen(ratland, print.progress = F)         
  gdf <- geomorph.data.frame(rat.gpa)  
  succeed(fit <- procD.lm(coords ~ Csize, data = gdf, iter = 3, verbose = TRUE,
                          RRPP = TRUE, print.progress = FALSE))
  succeed(summary(fit))
  succeed(plot(fit, type = "diagnostics"))
  succeed(plot(fit, type = "diagnostics", outliers = TRUE))
  succeed(plot(fit, type = "PC", pch = 19, col = "blue"))
  succeed(plot(fit, type = "regression", 
               predictor = gdf$Csize, reg.type = "RegScore", 
               pch = 19, col = "green"))
  succeed(rat.plot <- plot(fit, type = "regression", 
                           predictor = gdf$Csize, reg.type = "RegScore", 
                           pch = 21, bg = "yellow"))
  succeed(preds <- shape.predictor(fit$GM$fitted, x = rat.plot$RegScore, 
                                   predmin = min(rat.plot$RegScore), 
                                   predmax = max(rat.plot$RegScore)))
  M <- rat.gpa$consensus
  succeed(plotRefToTarget(M, preds$predmin, mag=2))
  succeed(plotRefToTarget(M, preds$predmax, mag=2))
})

test_that("procD.lm.example5.works", {  
  data("larvalMorph")
  Y.gpa <- gpagen(larvalMorph$tailcoords, 
                  curves = larvalMorph$tail.sliders,
                  ProcD = TRUE, print.progress = FALSE)
  gdf <- geomorph.data.frame(Y.gpa, treatment = larvalMorph$treatment, 
                             family = larvalMorph$family)
  succeed(fit <- procD.lm(coords ~ treatment/family, data = gdf, verbose = TRUE,
                          print.progress = FALSE, iter = 3))
  succeed(anova(fit))
  succeed(anova(fit, error = c("treatment:family", "Residuals")))
})


### procD.pgls --------------------------------------------------------------

test_that("procD.pgls.example1.works", {
  data(plethspecies)
  Y.gpa <- gpagen(plethspecies$land, print.progress = F)
  gdf <- geomorph.data.frame(Y.gpa, phy = plethspecies$phy)
  succeed(pleth.pgls <- procD.pgls(coords ~ Csize, phy = phy, 
        iter = 3, data = gdf, print.progress = F))
  succeed(anova(pleth.pgls))
  succeed(summary(pleth.pgls))
})

test_that("procD.pgls.example2.works", {
  data(plethspecies)
  Y.gpa <- gpagen(plethspecies$land, print.progress = F)
  gdf <- geomorph.data.frame(Y.gpa, phy = plethspecies$phy)
  pleth.pgls <- procD.pgls(coords ~ Csize, phy = phy, 
                                   iter = 3, data = gdf, print.progress = F)
  succeed(predict(pleth.pgls))
  succeed(plot(pleth.pgls, type = "regression", reg.type = "RegScore", 
    predictor = gdf$Csize))
  succeed(pleth.pgls$LM$Pcov)
  succeed(pleth.pgls2 <- procD.pgls(coords ~ Csize, phy = phy, lambda = 0.5, 
   iter = 3, data = gdf))
  succeed(anova(pleth.pgls))
  succeed(anova(pleth.pgls2))
})

### rotate.coords --------------------------------------------------------------

test_that("rotate.coords1.works", {
  data(plethodon)
  Y.gpa <- gpagen(plethodon$land, print.progress = F)
  succeed(Y.gpa2 <- rotate.coords(Y.gpa, "flipX"))
  succeed(plot(Y.gpa2))
  succeed(Y.gpa3 <- rotate.coords(Y.gpa2, "rotateCC"))
  succeed(plot(Y.gpa3))
  spec1 <- Y.gpa$coords[,,1]
  succeed(spec1 <- rotate.coords(spec1, "flipY"))
})

### shape.predictor --------------------------------------------------------------

test_that("shape.predictor1.works", {
  data("plethodon")
  Y.gpa <- gpagen(plethodon$land, print.progress = F)      
  succeed(preds <- shape.predictor(Y.gpa$coords, x = NULL, Intercept = FALSE, 
    pred1 = -0.1, pred2 = 0.1))
  M <- mshape(Y.gpa$coords)
  succeed(plotRefToTarget(M, preds$pred1))
  succeed(plotRefToTarget(M, preds$pred2))
})

test_that("shape.predictor2.works", {
  data("plethodon")
  Y.gpa <- gpagen(plethodon$land, print.progress = F)     
  M <- mshape(Y.gpa$coords)
  PCA <- gm.prcomp(Y.gpa$coords)
  PC <- PCA$x[,1]
  succeed(preds <- shape.predictor(Y.gpa$coords, x = PC, Intercept = FALSE, 
    pred1 = min(PC), pred2 = max(PC)))
  succeed(plotRefToTarget(M, preds$pred1))
  succeed(plotRefToTarget(M, preds$pred2))
  PC <- PCA$x[,1:2]
  succeed(preds <- shape.predictor(Y.gpa$coords, x = PC, Intercept = FALSE, 
    pred1 = c(0.045,-0.02),   pred2 = c(-0.025,0.06), pred3 = c(-0.06,-0.04))) 
  succeed(plotRefToTarget(M, preds$pred1))
  succeed(plotRefToTarget(M, preds$pred2))
  succeed(plotRefToTarget(M, preds$pred3))
})

test_that("shape.predictor3.works", {
  data("plethodon")
  Y.gpa <- gpagen(plethodon$land, print.progress = F)   
  M <- mshape(Y.gpa$coords)
  succeed(preds <- shape.predictor(Y.gpa$coords, x = log(Y.gpa$Csize), 
    Intercept = TRUE,  predmin = min(log(Y.gpa$Csize)), 
    predmax = max(log(Y.gpa$Csize)))) 
  succeed(plotRefToTarget(M, preds$predmin, mag = 3))
  succeed(plotRefToTarget(M, preds$predmax, mag = 3))
})

test_that("shape.predictor4.works", {
  data("plethodon")
  Y.gpa <- gpagen(plethodon$land, print.progress = F) 
  M <- mshape(Y.gpa$coords)
  gdf <- geomorph.data.frame(Y.gpa)
  plethAllometry <- procD.lm(coords ~ log(Csize), iter = 3, data = gdf, 
                             print.progress = F)
  allom.plot <- plot(plethAllometry, type = "regression", 
    predictor = log(gdf$Csize), reg.type ="RegScore") 
  succeed(preds <- shape.predictor(plethAllometry$GM$fitted, 
      x = allom.plot$RegScore, Intercept = FALSE, 
      predmin = min(allom.plot$RegScore), 
      predmax = max(allom.plot$RegScore))) 
  succeed(plotRefToTarget(M, preds$predmin, mag = 3))
  succeed(plotRefToTarget(M, preds$predmax, mag = 3))
  succeed(preds <- shape.predictor(plethAllometry$GM$fitted, 
    x = allom.plot$PredLine, Intercept = FALSE, 
    predmin = min(allom.plot$PredLine), 
    predmax = max(allom.plot$PredLine))) 
  succeed(plotRefToTarget(M, preds$predmin, mag = 3))
  succeed(plotRefToTarget(M, preds$predmax, mag = 3))
})

test_that("shape.predictor5.works", {
  data("plethodon")
  Y.gpa <- gpagen(plethodon$land, print.progress = F) 
  M <- mshape(Y.gpa$coords)
  gdf <- geomorph.data.frame(Y.gpa, species = plethodon$species, 
            site = plethodon$site)
  pleth <- procD.lm(coords ~ species * site, data=gdf, iter = 3, print.progress = F)
  PCA <- prcomp(pleth$fitted)
  means <- unique(round(PCA$x, 3))
  succeed(preds <- shape.predictor(arrayspecs(pleth$fitted, 12,2), x = PCA$x[,1:3],
    Intercept = FALSE, pred1 = means[1,1:3],   pred2 = means[2,1:3],
    pred3 = means[3,1:3],  pred4 = means[4,1:3]))                   
  succeed(plotRefToTarget(M, preds$pred1, mag = 2))
  succeed(plotRefToTarget(M, preds$pred2, mag = 2))
  succeed(plotRefToTarget(M, preds$pred3, mag = 2))
  succeed(plotRefToTarget(M, preds$pred4, mag = 2))
})

test_that("shape.predictor6.works", {
  data("plethodon")
  Y.gpa <- gpagen(plethodon$land, print.progress = F) 
  M <- mshape(Y.gpa$coords)
  gdf <- geomorph.data.frame(Y.gpa, species = plethodon$species, 
                             site = plethodon$site)
  pleth <- procD.lm(coords ~ species * site, data=gdf, iter = 3, print.progress = F)
  X <- pleth$X
  X <- X[,-1]
  symJord <- c(0,1,0) 
  alloJord <- c(0,0,0) 
  succeed(preds <- shape.predictor(arrayspecs(pleth$fitted, 12,2), x = X, 
    Intercept = TRUE, symJord=symJord, alloJord=alloJord))
  succeed(plotRefToTarget(M, preds$symJord, mag = 2))
  succeed(plotRefToTarget(M, preds$alloJord, mag = 2))
})

test_that("shape.predictor7.works", {
  data(plethShapeFood) 
  Y.gpa <- gpagen(plethShapeFood$land, print.progress = F)
  PLS <- two.b.pls(A1 = plethShapeFood$food, A2 = Y.gpa$coords, iter = 3, 
                   print.progress = F) 
  succeed(preds <- shape.predictor(Y.gpa$coords, plethShapeFood$food, 
    Intercept = FALSE,  method = "PLS", 
    pred1 = 2, pred2 = -4, pred3 = 2.5)) 
  M <- mshape(Y.gpa$coords)
  succeed(plotRefToTarget(M, preds$pred1, mag = 2))
  succeed(plotRefToTarget(M, preds$pred2, mag = 2))
  succeed(plotRefToTarget(M, preds$pred3, mag = 2))
})

### shapeHulls --------------------------------------------------------------

test_that("shapeHulls1.works", {
  data("pupfish")
  gdf <- geomorph.data.frame(coords = pupfish$coords, Sex = pupfish$Sex,
    Pop = pupfish$Pop)
  fit <- procD.lm(coords ~ Pop * Sex, data = gdf, 
      iter = 3, print.progress = FALSE)
  succeed(pc.plot <- plot(fit, type = "PC", pch = 19))
  succeed(shapeHulls(pc.plot))
  succeed(pc.plot <- plot(fit, type = "PC", pch = 19))
  groups <- interaction(gdf$Pop, gdf$Sex)
  succeed(shapeHulls(pc.plot, groups = groups, 
    group.cols = c("dark red", "dark red", "dark blue", "dark blue"),
    group.lwd = rep(2, 4), group.lty = c(2, 1, 2, 1)))
  succeed(pc.plot <- plot(fit, type = "PC", pch = 19))
  succeed(shapeHulls(pc.plot, groups = gdf$Sex, group.cols = c("black", "black"), 
    group.lwd = rep(2, 2), group.lty = c(2, 1)))
})

test_that("shapeHulls2.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land, print.progress = F)    #GPA-alignment
  pleth.phylo <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy)
  summary(pleth.phylo)
  succeed(pc.plot <- plot(pleth.phylo, phylo = TRUE))
  gp <- factor(c(rep(1, 5), rep(2, 4)))
  succeed(shapeHulls(pc.plot, groups = gp, group.cols = 1:2, 
    group.lwd = rep(2, 2), group.lty = c(2, 1)))
})

### two.b.pls --------------------------------------------------------------

test_that("two.b.pls1.works", {
  data(plethShapeFood) 
  Y.gpa <- gpagen(plethShapeFood$land, print.progress = F)  
  succeed(PLS <- two.b.pls(Y.gpa$coords, plethShapeFood$food, iter = 3, 
                           print.progress = F))
  succeed(summary(PLS))
  succeed(P <- plot(PLS))
})

### two.d.array --------------------------------------------------------------

test_that("two.d.array1.works", {
  data(plethodon) 
  succeed(two.d.array(plethodon$land))   
})

### physignal.eigen --------------------------------------------------------------

test_that("physignal.eigen.works", {
  data(plethspecies) 
  Y.gpa <- gpagen(plethspecies$land, print.progress = F)
  succeed(PSe.shape <- physignal.eigen(Y = Y.gpa$coords, phy = plethspecies$phy,
                                       iter = 3))
  succeed(summary(PSe.shape))
  succeed(plot(PSe.shape))
})

### extended.pgls --------------------------------------------------------------

#test_that("extended.pgls.works", {
#  data(pupfish.ws) 
#  succeed(fit <- extended.pgls(f1 = coords~Species * Sex + Population, 
#                               data = pupfish.ws, species = "Species",
#                               phy = pupfish.ws$phy, iter = 3))
#  succeed(anova(fit))
#  succeed(fit.mult <- manova.update(fit, PC.no = 40))
#  succeed(summary(fit.mult, test = "Wilks"))
#  succeed(fit2 <- extended.pgls(f1 = coords ~ Species * Sex + Population, 
#                                data = pupfish.ws, species = "Species",
#                                Cov = pupfish.ws$Cov, iter = 3))
#  succeed(anova(fit2))
#  succeed(fit2.mult <- manova.update(fit2, PC.no = 40))
#  succeed(summary(fit2.mult, test = "Wilks"))
#})
# TOO LONG RUNTIME