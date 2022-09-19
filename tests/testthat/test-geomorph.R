### -----------------------------------------------------------------------
###
### All parallel processing options must be run outside of CRAN check
###
### All tests that require dependencies should be run outside of CRAN check
###
### -----------------------------------------------------------------------


### procD.lm examples --------------------------------------------------------------

test_that("procD.lm.examples.works", {
  library(geomorph)
  data(plethodon) 
  Y.gpa <- gpagen(plethodon$land)    #GPA-alignment  
  gdf <- geomorph.data.frame(Y.gpa, 
                             site = plethodon$site, 
                             species = plethodon$species) # geomorph data frame
  
  succeed(fit1 <- procD.lm(coords ~ species * site, 
                   data = gdf, iter = 3, turbo = TRUE,
                   RRPP = FALSE, print.progress = FALSE))
  succeed(fit2 <- procD.lm(coords ~ species * site, 
                   data = gdf, iter = 3, turbo = TRUE,
                   RRPP = TRUE, print.progress = FALSE))
  
  succeed(summary(fit1))
  succeed(summary(fit2))
  succeed(coef(fit2))
  succeed(coef(fit2, test = TRUE))
  succeed(anova(fit2)) 
  succeed(anova(fit2, effect.type = "Rsq"))
  succeed(anova(fit2, error = c("species:site", "species:site", "Residuals")))
  M <- Y.gpa$consensus
  succeed(plotRefToTarget(M, fit2$GM$fitted[,,1], mag = 3))
  succeed(plotRefToTarget(M, fit2$GM$fitted[,,20], mag = 3))
  
  succeed(fit.null <- procD.lm(coords ~ log(Csize) + species + site, data = gdf, 
                       iter = 3, print.progress = FALSE, turbo = TRUE))
          succeed(fit.full <- procD.lm(coords ~ log(Csize) + species * site, data = gdf, 
                       iter = 3, print.progress = FALSE, turbo = TRUE))
  succeed(anova(fit.null, fit.full, print.progress = FALSE))
  succeed(gp <-  interaction(gdf$species, gdf$site))
  succeed(PW <- pairwise(fit.full, groups = gp, covariate = NULL))
  succeed(summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE))
  succeed(summary(PW, test.type = "dist", confidence = 0.95, stat.table = FALSE))
  succeed(summary(PW, test.type = "var", confidence = 0.95, stat.table = TRUE))
  succeed(summary(PW, test.type = "var", confidence = 0.95, stat.table = FALSE))
  succeed(morphol.disparity(fit.full, groups = gp, iter = 3))
  data(ratland)
  rat.gpa<-gpagen(ratland)         #GPA-alignment
  gdf <- geomorph.data.frame(rat.gpa) # geomorph data frame is easy 
  succeed(fit <- procD.lm(coords ~ Csize, data = gdf, iter = 3, turbo = TRUE,
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

  data("larvalMorph")
  Y.gpa <- gpagen(larvalMorph$tailcoords, 
                  curves = larvalMorph$tail.sliders,
                  ProcD = TRUE, print.progress = FALSE)
  gdf <- geomorph.data.frame(Y.gpa, treatment = larvalMorph$treatment, 
                             family = larvalMorph$family)
  
  succeed(fit <- procD.lm(coords ~ treatment/family, data = gdf, turbo = TRUE,
                  print.progress = FALSE, iter = 3))
  succeed(anova(fit))
  succeed(anova(fit, error = c("treatment:family", "Residuals")))
  
})

# Other tests that are not examples can also go here.
