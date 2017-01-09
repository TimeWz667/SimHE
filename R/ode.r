ODE <- function(m, pars, default.cost=1, default.qol=1, default.trs=1) {

  if (!check.parameters(m, pars)) {
    stop('Parameter Checking Error')
  }

  m <- m$model.def

  n.states <- nrow(m$States)
  name.states <- rownames(m$States)

  im <- matrix(pars[m$Targets], n.states, n.states, dimnames=list(name.states, name.states))
  im[is.na(im)] <- 0
  im[is.na(im) & m$Targets != ''] <- default.trs
  im <- im - diag(rowSums(im))

  costs <- pars[m$States[,'Cost']]
  names(costs) <- name.states
  costs[is.na(costs)] <- default.cost

  qols <- pars[m$States[,'QOL']]
  names(qols) <- name.states
  qols[is.na(qols)] <- default.qol

  delta <- function(Time, state, Pars) {
    return (list(c(state %*% im)))
  }

  dy <- list(Delta=delta, Costs=costs, QOLs=qols)
  class(dy) <- 'ODE'
  dy
}


goto.ODE <- function(model, y, ti, dt=1) {
  ts <- seq(ti, ti+dt, length.out=10)
  y <- deSolve::ode(y, ts, model$Delta)[length(ts), -1]
  y
}

observe.ODE <- function(model, y, ti) {
  c(Time=ti, y, N=sum(y), Cost=sum(y*model$Costs), QOL=sum(y*model$QOLs))
}


