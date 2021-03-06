% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uncertainty.r
\name{simMultiHE}
\alias{simMultiHE}
\title{Simulate model with multiple input.}
\usage{
simMultiHE(mfunc, md, pars.tab, yini, fr, to, dt = 1, n.sim = 100)
}
\arguments{
\item{md}{model definition}

\item{pars.tab}{list of probabilistic distributions of parameters, double value if constant}

\item{yini}{initial value for simulation}

\item{fr}{start time}

\item{to}{end time}

\item{dt}{observation interval}

\item{n.sim}{size of parameter sets}

\item{nfunc}{Constructor of model, (DetMC, StoMC, ODE, IBM)}
}
\value{
list of outputs of each simulation
}
\description{
Simulate model with multiple input.
}
\examples{


}

