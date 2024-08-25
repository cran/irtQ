# set starting value for the item parameter estimation
set_startval <- function(pars = NULL, item, use.startval = FALSE, mod, score.cat,
                         fix.a.1pl = FALSE, fix.g = FALSE, fix.a.gpcm = FALSE,
                         n.1PLM = NULL) {
  if (use.startval) {
    if (mod == "1PLM") {
      if (fix.a.1pl) {
        startval <- pars[item, 2]
      } else {
        startval <- c(pars[item[1], 1], pars[item, 2])
      }
    } else if (mod == "2PLM" | (mod == "3PLM" & fix.g)) {
      startval <- pars[item, 1:2]
    } else if (mod == "3PLM" & !fix.g) {
      startval <- pars[item, 1:3]
    } else if (mod == "GRM" | (mod == "GPCM" & !fix.a.gpcm)) {
      startval <- pars[item, 1:score.cat]
    } else if (mod == "GPCM" & fix.a.gpcm) {
      startval <- pars[item, 2:score.cat]
    }
  } else {
    if (mod == "1PLM") {
      if (fix.a.1pl) {
        startval <- 0
      } else {
        startval <- c(1, rep(0, n.1PLM))
      }
    } else if (mod == "2PLM" | (mod == "3PLM" & fix.g)) {
      startval <- c(1, 0)
    } else if (mod == "3PLM" & !fix.g) {
      startval <- c(1, 0, 0.2)
    } else if (mod == "GRM" | (mod == "GPCM" & !fix.a.gpcm)) {
      startval <- c(1, seq(-1.0, 1.0, length.out = (score.cat - 1)))
    } else if (mod == "GPCM" & fix.a.gpcm) {
      startval <- seq(-1.0, 1.0, length.out = (score.cat - 1))
    }
  }

  # return the results
  startval
}
