
#' Summarise cost and qol with time preference
#'
#' @param out output table from simulation model
#' @param rho_cost interest rate for future cost
#' @param rho_qol interst rate for future qol
#'
#' @return a vector of cost and qol
#' @export
#'
#' @examples
#'
summarise.CE <- function(out, rho_cost=0.02, rho_qol=0.02) {
  obs <- out$Obs
  return (c(
    'Cost'=exp(-obs$Time*rho_cost) %*% obs$Cost,
    'QOL'=exp(-obs$Time*rho_qol) %*% obs$QOL
  ))
}


#' Simulate model with multiple input.
#'
#' @param n.sim size of parameter sets
#' @param nfunc Constructor of model, (DetMC, StoMC, ODE, IBM)
#' @param md model definition
#' @param pars.tab0 list of probabilistic distributions of parameters for standard treatment
#' @param pars.tab1 list of probabilistic distributions of parameters for new intervention
#' @param yini initial value for simulation
#' @param rho_cost interest rate for future cost
#' @param rho_qol interst rate for future qol
#' @param fr start time
#' @param to end time
#' @param dt observation interval
#'
#' @return object of cea
#' @export
#'
#' @examples
#'
#'
cea <- function(mfunc, md, pars.tab0, pars.tab1, yini, rho_cost, rho_qol, fr, to, dt=1, n.sim=100) {
  cat('Simulating standard model\n')
  outs0 <- simMultiHE(mfunc, md, pars.tab0, yini, fr, to, dt, n.sim)
  cat('Simulating intervention model\n')
  outs1 <- simMultiHE(mfunc, md, pars.tab1, yini, fr, to, dt, n.sim)

  ce0 <- sapply(outs0$Simulations, function(out) summarise.CE(out, rho_cost, rho_qol))
  ce1 <- sapply(outs1$Simulations, function(out) summarise.CE(out, rho_cost, rho_qol))

  delta <- data.frame(t(ce1-ce0))

  print(t.test(ce0['QOL',], ce1['QOL',])$p.value)
  if (t.test(delta$QOL)$p.value > 0.05) {
    cat('Zero division expected. Can not calculate ICER\n')
    icer <- NA
  } else {
    icer <- delta['Cost']/delta['QOL']
  }

  lambda <- quantile(delta$Cost/delta$QOL, 0.95)
  if (mean(delta$QOL - delta$Cost/lambda > 0) >= 0.95) {
    cat('pr[Net Health Benefit > 0] > 0.95 if willingness to pay > ', lambda, '\n')
    nhb <-  c(summary(delta$QOL - delta$Cost/lambda))
  } else {
    cat('Intervention never cost effectiveness\n')
    nhb <- 0
  }

  output <- list(
    'CE0'=data.frame(t(ce0)),
    'CE1'=data.frame(t(ce1)),
    'dC'=delta$Cost,
    'dE'=delta$QOL,
    'ICER'=icer,
    'NHB'=nhb,
    'WTP'=ifelse(is.na(nhb), 0, lambda)
  )


  class(output) <-'CEA'

  return(output)
}

