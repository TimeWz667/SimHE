#' Stochastic Markov Chain Model
#' create a Stochastic Markov model. Basicly, the model is defined under continuous time
#' @param dm definiation of model
#' @param pars list or vector of parameters
#' @param default.cost default value of cost per simulation interval per person
#' @param default.qol default value of qol per simulation interval per person, ranged [0, 1]
#' @param default.trs default value of intensity (transition rate)
#'
#' @return Stochastic Markov Chain Model for simulation
#' @export
#'
#' @examples
#'
StoMC <- function(m, pars, default.cost=1, default.qol=1, default.trs=1) {

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

  dy <- list(IntensityMatrix=im, Costs=costs, QOLs=qols)
  class(dy) <- 'StoMC'
  dy
}


goto.StoMC <- function(model, y, ti, dt=1) {
  pm <- Matrix::expm(model$IntensityMatrix*dt)
  sts <- names(y)
  y <- sapply(1:length(y), function(i) stats::rmultinom(1, y[i], pm[i, ]))
  y <- rowSums(y)
  names(y) <- sts
  y
}


observe.StoMC <- function(model, y, ti) {
  c(Time=ti, y, N=sum(y), Cost=sum(y*model$Costs), QOL=sum(y*model$QOLs))
}
