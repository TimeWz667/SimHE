## trimming utility functions
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
trimall <- function (x) gsub("\\s+", "", x)
gfun <- function(x) gsub( "[:space:]*\\|.*'[:space:]*]","']",x) #for dropping Q/C


wrap4diag <- function(x){
  ## get rid of comments
  lnz <- unlist(strsplit(trim(x),split='\n')) #lines
  lnz <- lnz[lnz!=""]               #ditch empty lines
  lnz <- unlist(lapply(as.list(lnz),FUN=function(x) unlist(strsplit(x,split="\\\\"))[1] )) ## strip comment
  x <- paste0(lnz,collapse="\n")
  ## top/tail and substitute
  ssb <- "digraph dot{\ngraph [layout = 'dot',rankdir=LR];\nnode[shape=record]\n"
  sse <- "\n}"
  x <- gsub("\\[","[label=",x)
  x <- paste0(ssb,x,sse)
  x
}


## parse the string and extract the data
parseData <- function(x){
  lnz <- unlist(strsplit(trim(x),split='\n')) #lines
  lnz <- lnz[lnz!=""]               #ditch empty lines
  lnz <- unlist(lapply(as.list(lnz),FUN=function(x) unlist(strsplit(x,split="\\\\"))[1] )) ## strip comment
  lnz <- trimall(lnz)
  lnz <- lnz[lnz!=""]                 #ditch empty lines
  edz <- grepl("->",lnz)
  edges <- lnz[edz]
  nodes <- lnz[!edz]
  ## ---edges---
  espl <- strsplit(edges,split='->')
  efrom <- unlist(lapply(espl,function(x)x[[1]]))  #the froms
  eto <- unlist(lapply(espl, function(x) unlist(strsplit(x[2],split='\\['))[1])) #the tos
  erate <- unlist(lapply(espl, function(x) unlist(strsplit(x[2],split="'"))[2])) #the probabilities
  ## ---nodes---
  nnames <- unlist(lapply(nodes, function(x) unlist(strsplit(x,split='\\['))[1])) #the node names
  nct <- unlist(lapply(nodes, function(x) unlist(strsplit(x,split="'"))[2])) #node data
  ndat <- strsplit(nct,split="\\|")     #node data
  lbz <- unlist(lapply(ndat,function(x)x[1]))
  cstz <- unlist(lapply(ndat,function(x)x[2]))
  qlz <- unlist(lapply(ndat,function(x)x[3]))
  ## return values
  list(edges=list(from=efrom, to=eto, rate=erate),
       nodes=list(names=nnames, labels=lbz, costs=cstz, qols=qlz))
}

wrap4dynamicmodel <- function(dat) {
  name.states <- dat$nodes$names
  n.states <- length(name.states)

  States <- matrix("", n.states, 2, dimnames=list(name.states, c("QOL", "Cost")))
  for (i in 1:n.states) {
    States[i,] <- c(dat$nodes$qols[i], dat$nodes$costs[i])
  }

  Targets <- matrix("", n.states, n.states, dimnames=list(name.states, name.states))
  for (i in 1:length(dat$edges$from)) {
    Targets[dat$edges$from[i], dat$edges$to[i]] <- dat$edges$rate[i]
  }

  list(
    States=States,
    Targets=Targets
  )
}





#' Extract information in model definition script.
#'
#' The algorithm was a modified version of dtree package (https://github.com/petedodd/dtree/)
#'
#' @param ss script
#'
#' @return list of \itemize{
#'   \item \code{dot} {string in DOT language for export and visualization}
#'   \item \code{dotc} {string in DOT language for export and visualization (cost/qol node info dropped for conceptual overview)}
#'   \item \code{model.def} {Definition of dynamic model}
#'   \item \code{(todo) parmeters}{list of parameters}
#' }
#' @export
#'
#' @examples
#'  script0 <- "
#'  \\nodes...(this line is  comment)
#'    \\ nodename ['label | cost variable | QoL variable']
#'    A[ 'description A | costA | qolA'] \\ a node
#'    B[ 'description B | costB | qolB']
#'  C[  'description C | costC | qolC']
#'  \\edges
#'  \\ from_statename -> to_statename['transition function or name for this edge']
#'  A -> B ['TrAB']
#'  A -> C ['TrAC']
#'  B -> A ['TrBC']
#'
#'  \\ (whitespace is ignored in parsing for calculation)
#'  "
#'
#'  md = ddm(script0)
#'
ddm <- function(ss) {
  ssv <- wrap4diag(ss)                #graph code
  ssc <- paste( unlist(lapply( unlist(strsplit(trim(ssv),split='\n')), gfun)), collapse='\n')
  dm <- wrap4dynamicmodel(parseData(ss))
  list(
    dot=ssv,
    dotc=ssc,
    model.def=dm
  )
}


#' Check the values of parameters
#' print the ill-defined parameters and parameters applied default values
#'
#' @param dm definition sheet of model
#' @param pars a vector or list of the model
#'
#' @return True if checking failed
#' @export
#'
#' @examples
#'
check.parameters <- function(dm, pars) {
  m <- dm$model.def

  err <- c()
  fail <- FALSE
  par.in <- c(m$States[, 'Cost'], m$States[, 'QOL'], c(m$Targets))
  par.in <- unique(par.in)
  par.in <- par.in[par.in != '']
  for (name in par.in) {
    val <- pars[name]
    if (is.na(val)) {
      err <- c(err, paste0('- ', name, ' uses default value'))
    } else if (val < 0) {
      err <- c(err, paste0('- ', name, ' is a negative value'))
      fail <- TRUE
    }
  }

  for (er in err) cat(er, '\n')

  cat('The Parameters', ifelse(fail, 'are not', 'are'),'well-placed\n')
  return (!fail)
}


#' Simulate forward one unit
#'
#' @param model dynamic model
#' @param y initial value
#' @param ti initial time
#' @param dt step size
#'
#' @return y after simulation
#' @export
#'
#' @examples
#'
goto <- function(model, y, ti, dt) {
  UseMethod("goto", model)
}


#' Make observations of dynamic model
#'
#' @param model dynamic model
#' @param y status of the model
#' @param ti current time
#'
#' @return vector of obasevations
#' @export
#'
#' @examples
#'
observe <- function(model, y, ti, ...) {
  UseMethod("observe", model)
}


update.model <- function(model, y, ti, dt, rec='m',...) {
  if (rec=='s') {
    obs <- observe(model, y, ti, ...)
    y <- goto(model, y, ti, dt)
  } else if (rec=='m') {
    y <- goto(model, y, ti, dt/2)
    obs <- observe(model, y, ti+dt/2, ...)
    y <- goto(model, y, ti+dt/2, ti+dt)
  } else {
    y <- goto(model, y, ti, dt)
    obs <- observe(model, y, ti+dt, ...)
  }
  return (list(Time=ti, Y=y, Observations=obs))
}


#' Simulate observations from dynamic model
#'
#' @param model dynamic model
#' @param yini initial values
#' @param from double, start time
#' @param to double, end time
#' @param dt double, observation interval
#' @param rec characher, timing for observation, (s: starting, m: middle, e: ending)
#'
#'
#' @return \code{Y}: last value in the end of simulation
#' @return \code{Obs}: the observed time series (data.frame)
#' @export
#'
#' @examples
#'
simHE <- function(model, yini, from, to, dt=1, rec='m', ...) {
  UseMethod('simulate', model)
}

simHE.default <- function(model, yini, from, to, dt=1, rec='m', ...) {
  times <- seq(from, to, by=dt)
  if (rec != 's') {
    observations <- observe(model, yini, from, ...)
  } else {
    observations <- c()
  }

  y <- yini
  ti.sim <- proc.time()
  for (i in 1:(length(times) - 1)) {
    fr <- times[i]
    dt <- times[i+1] - fr
    tmp <- update.model(model, y, fr, dt, rec, ...)
    y <- tmp$Y
    observations <- rbind(observations, tmp$Observations)
  }
  ti.sim <- proc.time() - ti.sim
  if (rec != 'e') observations <- rbind(observations, observe(model, y, to, ...))
  row.names(observations) <- 1:nrow(observations)
  res <- list(Y=y, Obs=data.frame(observations), ModelType=class(model), SimTi=ti.sim)
  class(res) <- c('OutputDM', paste0('Output', class(model)))
  res
}


#' @export
print.OutputDM <- function(out,...) {
  cat('Results: ')
  print(out$Obs)
  cat('\nSimulation Time: \n')
  print(out$SimTi)
}




