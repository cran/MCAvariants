\name{satisfaction}
\alias{satisfaction}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Patient Satisfaction }
\description{
The data set consists of 235 rows and 5 columns.
The rows represent the individuals (patients in an hospital) and the columns concern the five variables of satisfaction 
(Tangibility, Reliability, Capacity of Response, Capacity of Assurance and Empathy
}
\usage{data(satisfaction)}
\format{
  A data frame with 235 observations on the following 5 variables.
  \describe{
    \item{\code{TANG}}{a numeric vector of ordered categories ranging from 1 to 5.}
    \item{\code{REL}}{a numeric vector of ordered categories ranging from 1 to 5.}
    \item{\code{CRES}}{a numeric vector of ordered categories ranging from 1 to 5.}
    \item{\code{CASS}}{a numeric vector of ordered categories ranging from 1 to 5.}
    \item{\code{EMPAT}}{a numeric vector of ordered categories ranging from 1 to 5.}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
Della Valle E (2010) Della Valle E 2010 Analisi Multidimensionale dei Dati: la Qualit\'{a} del Lavoro nelle Cooperative Sociali. 
Unpublished thesis, Seconda Universit\'{a}  degli Studi di Napoli, Capua Italy.
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(satisfaction)
#dim(satisfaction)
#dimnames(satisfaction)
}
\keyword{datasets}
