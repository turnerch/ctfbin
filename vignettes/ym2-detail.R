## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(bozesim)

## -----------------------------------------------------------------------------
filter.example <- synthetic_ym2_pf()

## -----------------------------------------------------------------------------
filter.example2 <- synthetic_ym2_pf(seed = 10000)

## -----------------------------------------------------------------------------
GPP(lai = 2.641, c_eff = 20, maxt = 12.83247, mint = 6.867838, Ca = 400.95, lat = 46.23, yearday = 268, rad = 11.26568)

## -----------------------------------------------------------------------------
params.scaled <- synthetic_ym2_pars()
params.unscaled <- synthetic_ym2_pars(scaled = FALSE)

## -----------------------------------------------------------------------------
filter.example <- synthetic_ym2_pf()
params.scaled <- synthetic_ym2_pars()

dalec <- synthetic_ym2(filter = filter.example, pars = params.scaled)

## -----------------------------------------------------------------------------
filter.example2 <- synthetic_ym2_pf(seed = 10000)

dalec2 <- synthetic_ym2(filter = filter.example2, pars = params.scaled)

## -----------------------------------------------------------------------------
params.unscaled <- synthetic_ym2_pars(scaled = FALSE)

dalec3 <- synthetic_ym2(filter = filter.example, pars = params.unscaled, scaled = FALSE)

## -----------------------------------------------------------------------------
dalec <- synthetic_ym2(filter = filter.example, pars = params.scaled, verbose = FALSE)

## -----------------------------------------------------------------------------
dalec4 <- synthetic_ym2(filter = filter.example, pars = params.scaled, n.part = 300, verbose = FALSE)

## -----------------------------------------------------------------------------
dalec$logLik       # marginal log-likelihood
plot(dalec$traj)   # trajectory

