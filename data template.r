rm(list=ls())
library(SimHE)

script0 <- "
\\nodes...(this line is  comment)
\\ nodename ['label | cost variable | QoL variable']
A[ 'description A | costA | qolA'] \\ a node
B[ 'description B | costB | qolB']
C[  'description C | costC | qolC']
\\edges
\\ from_statename -> to_statename['transition function or name for this edge']
A -> B ['TrAB']
A -> C ['TrAC']
B -> C ['TrBC']

\\ (whitespace is ignored in parsing for calculation)
"

md = ddm(script0)

pars <- c(qolA=1, qolB=0.8, colC=0.5, costA=1, costB=2, TrAB=0.02, TrBC=0.5, TrAC=0.2)
yini <- c(A=900, B=100, C=0)


model <- DetMC(md, pars)
out <- simHE(model, yini, 0, 10, 1)


model <- StoMC(md, pars)
out <- simHE(model, yini, 0, 10, 1)


model <- ODE(md, pars)
out <- simHE(model, yini, 0, 10, 1)


pars <- list(qolA=1, qolB=0.8, colC=0.5, costA=1, costB=2, TrAB=0.02, TrBC=0.5,
          TrAC=function(ag, ...) rgamma(1, 0.01, 0.1))

model <- IBM(md, pars)
out <- simHE(model, yini, 0, 10, 1)



pars.sim <- list(qolA=1, qolB=0.8, colC=0.5, costA=1, costB=2, TrAB=0.02, TrBC='exp(0.2)',
     TrAC='gamma(0.1,0.1)')

outs <- simMultiHE(DetMC, md, pars.sim, yini, fr=0, to=10, dt=1, 100)










md = ddm(script0)
yini <- c(A=900, B=100, C=0)

pars.sim1 <- list(qolA=1, qolB=0.8, qolC=0.5, costA=1, costB=2, TrAB=0.02, TrBC='exp(0.2)',
                  TrAC='gamma(0.1,0.1)')

pars.sim2 <- list(qolA=1, qolB=0.8, qolC=0.5, costA=1, costB=2, TrAB=0.02, TrBC='exp(0.5)',
                  TrAC='gamma(0.05,0.1)')


ca <- cea(ODE, md, pars.sim2, pars.sim1, yini, 0.02, 0.02, fr=0, to=10, dt=1, 100)
