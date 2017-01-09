#' Individual-Based Model (microsimulation)
#' create an individual based model. The state
#' @param dm definiation of model
#' @param pars list or vector of parameters
#' @param default.cost default value of cost per simulation interval per person
#' @param default.qol default value of qol per simulation interval per person, ranged [0, 1]
#' @param default.trs default value of intensity (transition rate)
#' @param ini.attr function for generating attributes of an agent ()=>vector(Attr: Val)
#'
#' @return individual based model for simulation
#' @export
#'
#' @examples
#'
IBM <- function(m, pars, default.cost=1, default.qol=1, default.trs=1, ini.attr=F) {
  if (is.list(pars) & !check.parameters(m, pars)) {
    stop('Parameter Checking Error')
  }
  ## check list
  m <- m$model.def

  n.states <- nrow(m$States)
  name.states <- rownames(m$States)

  Agent <- function(id, state, ...) {
    return (c(list(ID=id, State=state, Next=state, TTE=Inf), ...))
  }

  make.agents <- function(y) {
    states <- sample(rep(names(y),y))
    if (ini.attr)
      ags <- lapply(1:sum(y), function(x) Agent(x, states[x], ini.attr()))
    else
      ags <- lapply(1:sum(y), function(x) Agent(x, states[x]))
    return (data.frame(rbindlist(ags)))
  }

  targets <- lapply(name.states, function(x) {
      tr <- m$Targets[x,]
      tr <- tr[tr != '']
      to <- names(tr)
      names(to) = tr
      to
    })
  names(targets) <- name.states
  trans <- pars[unlist(lapply(targets, function(x) names(x)))]
  trans[is.na(trans)] <- default.trs

  trans <- lapply(trans, function(tr) {
    if (is.numeric(tr)) {
      return (function(ag) rexp(1, rate=tr))
    } else if (is.function(tr)) {
      return (tr)
    } else if (is.data.frame(tr)) {
      return (function(ag) {
        atr <- tr[names(tr) != "Rate"]
        return (rexp(1, rate=tr[atr == ag[,names(atr)], "Rate"]))
      })
    }
  })

  find.next <- function(ag, ti=0) {
    ta <- targets[[ag$State]]
    if (length(ta) == 0) {
      ag$TTE <- Inf
    } else {
      tte <- sapply(names(ta), function(f) {trans[[f]](ag)})
      tr <- which.min(tte)
      ag$Next <- ta[tr]
      ag$TTE <- tte[tr] + ti
    }
    return (ag)
  }

  transit <- function(ag) {
    ag$State <- ag$Next
    return (ag)
  }

  initialise <- function(ags, ti) {
    n <- nrow(ags)
    ags <- lapply(1:n, function(i) find.next(ags[i, ], ti))
    return (data.frame(rbindlist(ags)))
  }

  costs <- pars[m$States[,'Cost']]
  names(costs) <- name.states
  costs[is.na(costs)] <- default.cost

  qols <- pars[m$States[,'QOL']]
  names(qols) <- name.states
  qols[is.na(qols)] <- default.qol

  dy <- list(make.agents=make.agents,
             initialise=initialise,
             find.next=find.next,
             transit=transit,
             Targets=targets,
             Costs=costs, QOLs=qols)
  class(dy) <- 'IBM'
  dy
}


goto.IBM <- function(model, y, ti, dt=1) {
  if (is.data.frame(y)) {
    ags <- y
  } else {
    ags <- model$make.agents(y)
    ags <- model$initialise(ags, ti)
  }

  ti.start <- ti
  ti.now <- ti.start
  ti.end <- ti+dt
  while (TRUE) {
    evt.i <- which.min(ags$TTE)
    evt.agent <- ags[evt.i,]
    evt.tte <- evt.agent$TTE
    if (evt.tte > ti.end) {
      break
    }
    ti.now <- evt.tte
    evt.agent <- model$transit(evt.agent)
    evt.agent <- model$find.next(evt.agent, ti.now)
    ags[evt.i,] <- evt.agent

  }
  return (ags)
}


observe.IBM <- function(model, y, ti) {
  if (is.data.frame(y)) {
    y <- table(factor(y$State, level=names(model$Targets)))
  }
  c(Time=ti, y, N=sum(y), Cost=sum(y*model$Costs), QOL=sum(y*model$QOLs))
}

