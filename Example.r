source("r/abm.r")

n.agents = 1000
time.start = 0
time.end = 10
time.diff = 1


States <- rbind(A = c("costA", "qolA"),
                B = c("costB", "qolB"),
                C = c("costC", "qolC"))

colnames(States) <- c("Cost", "Qol")


Transitions <- list(
  trAB = 0.05,
  trAC = data.frame(Rate=c(0.05, 0.1), Sex=c("F", "M")),
  trBC = function(ag) rgamma(1, 1, 0.1)
)

Futures <- list(
  A = c(trAB="B", trAC="C"),
  B = c(trBC="C"),
  C = c()
)


Ini.state <- function() {
  return (sample(c("A", "B", "C"), 1, prob=c(0.8, 0.2, 0)))
}

Ini.attr <- function() {
  return (list(Sex=sample(c("F", "M"), 1, prob=c(0.5, 0.5))))
}


abm <- make.abm(States, Transitions, Futures, Ini.state, Ini.attr)

simulate(abm, n.agents, time.start, time.end, time.diff)






