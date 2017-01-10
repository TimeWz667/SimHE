#' Deterministic Markov Chain Model
#' create a deterministic Markov model. Basicly, the model is defined under continuous time
#' @param dm definiation of model
#' @param pars list or vector of parameters
#' @param default.cost default value of cost per simulation interval per person
#' @param default.qol default value of qol per simulation interval per person, ranged [0, 1]
#' @param default.trs default value of intensity (transition rate)
#'
#' @return Deterministic Markov Chain Model for simulation
#' @export
#'
#' @examples
#'
DetMC <- function(dm, pars, default.cost=1, default.qol=1, default.trs=1) {

  if (!check.parameters(dm, pars)) {
    stop('Parameter Checking Error')
  }

  m <- dm$model.def

  n.states <- nrow(m$States)
  name.states <- rownames(m$States)

  im <- matrix(pars[m$Targets], n.states, n.states, dimnames=list(name.states, name.states))
  im[is.na(im)] <- 0
  im[is.na(im) & m$Targets != ''] <- default.trs
  im <- im - diag(rowSums(im))

  costs <- fill.pars.vector(m$States[, 'Cost'], pars, default.cost)

  qols <- fill.pars.vector(m$States[, 'QOL'], pars, default.qol)

  dy <- list(IntensityMatrix=im, Costs=costs, QOLs=qols)
  class(dy) <- 'DetMC'
  dy
}


goto.DetMC <- function(model, y, ti, dt=1) {
  pm <- Matrix::expm(model$IntensityMatrix*dt)
  y <- y %*%pm
  states <- as.numeric(y)
  names(states) <- colnames(y)
  states
}


observe.DetMC <- function(model, y, ti) {
  c(Time=ti, y, N=sum(y), Cost=sum(y*model$Costs), QOL=sum(y*model$QOLs))
}
