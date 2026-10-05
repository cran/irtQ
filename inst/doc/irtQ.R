## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6.5,
  fig.height = 4
)

## ----metadata-----------------------------------------------------------------
library(irtQ)

flex_file <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
x <- bring.flexmirt(file = flex_file, type = "par")$Group1$full_df
x <- x[1:25, ]
head(x)

## ----shape-df-----------------------------------------------------------------
shape_df(
  par.drm = list(a = c(1.0, 1.2), b = c(-0.5, 0.5), g = c(0.2, 0.2)),
  cats = 2, model = "3PLM"
)

## ----simulate-----------------------------------------------------------------
set.seed(2026)
theta <- rnorm(500, mean = 0, sd = 1)
resp <- simdat(x = x, theta = theta, D = 1)
dim(resp)

## ----calibrate----------------------------------------------------------------
fit <- est_irt(
  data = resp, D = 1, model = "3PLM", cats = 2, item.id = x$id,
  use.gprior = TRUE, gprior = list(dist = "beta", params = c(5, 16)),
  verbose = FALSE
)
fit

## ----par-est------------------------------------------------------------------
par_est <- getirt(fit, what = "par.est")
head(par_est)

## ----score--------------------------------------------------------------------
scores <- est_score(fit, method = "ML", range = c(-4, 4))
head(scores$est.theta)
cor(scores$est.theta, theta)

## ----sx2----------------------------------------------------------------------
sx2 <- sx2_fit(fit)
head(sx2$fit_stat)

## ----irtfit-plot--------------------------------------------------------------
fit_irt <- irtfit(
  x = fit, score = scores$est.theta, group.method = "equal.width",
  n.width = 10, loc.theta = "middle"
)
plot(fit_irt, item.loc = 1, type = "both", ci.method = "wald",
     show.table = FALSE, ylim.sr.adjust = TRUE)

## ----info---------------------------------------------------------------------
theta_grid <- seq(-4, 4, 0.1)
tif <- info(x = fit, theta = theta_grid, tif = TRUE)
plot(tif)

