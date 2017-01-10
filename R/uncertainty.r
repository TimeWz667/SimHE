#' Parameter table maker
#'
#' @param pars.sim list of probabilistic distributions of parameters, double value if constant
#' @param n.sim size of parameter sets
#'
#' @return data.frame of parameter table
#' @export
#'
#' @examples
#' pars.def <- list(TrAB=0.02, TrBC='exp(0.2)', TrAC='gamma(0.01, 0.1)')
#' make.parameters(pars.def, 10)
#'
make.parameters <- function(pars.sim, n.sim) {
  pars.tab <- lapply(pars.sim, function(x) {
    if (is.double(x)) {
      return (x)
    } else {
      x <- sub('\\(', paste0('(', n.sim, ','), paste0('r', x))
      return (eval(parse(text=x)))
    }})
  data.frame(pars.tab)
}


#' Simulate model with multiple input.
#'
#' @param n.sim size of parameter sets
#' @param nfunc Constructor of model, (DetMC, StoMC, ODE, IBM)
#' @param md model definition
#' @param pars.tab list of probabilistic distributions of parameters, double value if constant
#' @param yini initial value for simulation
#' @param fr start time
#' @param to end time
#' @param dt observation interval
#'
#' @return list of outputs of each simulation
#' @export
#'
#' @examples
#'
#'
simHE.multi <- function(mfunc, md, pars.tab, yini, fr, to, dt=1, n.sim=100) {
  if (!is.data.frame(pars.tab)) {
    pars.tab <- make.parameters(pars.sim, n.sim)
  } else {
    n.sim <- nrow(pars.tab)
  }

  prog <- txtProgressBar(0, n.sim, style=3)
  mods <- apply(pars.tab, 1, function(pa) mfunc(md, pa))

  sims <- lapply(mods, function(mod) {
    setTxtProgressBar(prog, getTxtProgressBar(prog)+1)
    simHE(mod, yini, fr, to, dt)
  })
  close(prog)
  list(Simulations=sims, Parameters=pars.tab)
}
