library(data.table)


make.abm <- function(state, transitions, futures, ini.state, ini.attr) {
  Agent <- function(id, state, ...) {
    return (c(list(ID=id, State=state, Next=state, TTE=Inf), ...))
  }

  rand.agents <- function(n) {
    ags <- lapply(1:n, function(x) Agent(x, ini.state(), ini.attr()))
    return (data.frame(rbindlist(ags)))
  }

  trans <- lapply(transitions, function(tr) {
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
    fu <- futures[[ag$State]]
    if (length(fu) == 0) {
      ag$TTE <- Inf
    } else {
      tte <- sapply(names(fu), function(f) {trans[[f]](ag)})
      tr <- which.min(tte)
      ag$Next <- fu[tr]
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

  return (list(
    rand.agents = rand.agents,
    initialise = initialise,
    find.next = find.next,
    transit = transit
  ))
}


simulate <- function (abm, n, t0, t1, dt=1) {
  with(abm, {
    ags <- rand.agents(n)
    ags <- initialise(ags, t0)

    ti.start <- t0
    ti.now <- ti.start
    out <- table(factor(ags$State, level=rownames(States)))

    while (ti.now < t1) {
      ti.end <- min(ti.start + dt, t1)
      while (TRUE) {
        evt.i <- which.min(ags$TTE)
        evt.agent <- ags[evt.i,]
        evt.tte <- evt.agent$TTE
        if (evt.tte > ti.end) {
          break
        }
        ti.now <- evt.tte
        evt.agent <- transit(evt.agent)
        evt.agent <- find.next(evt.agent, ti.now)
        ags[evt.i,] <- evt.agent
        ## make a record

      }
      out <- rbind(out, table(factor(ags$State, level=rownames(States))))
      ## make a record
      ti.start <- ti.end
      ti.now <- ti.end
    }

    return (out)
  })

}







