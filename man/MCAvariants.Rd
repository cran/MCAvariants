\name{MCAvariants}
\alias{MCAvariants}
\docType{package}
\title{
Classic and Ordered Multiple Correspondence Analysis
}
\description{
It performs Classic Multiple Correspondence analysis for nominal variables
(setting \code{catype}="mca") and Ordered Multiple Correspondence analysis via orthogonal polynomials (setting \code{catype}="omca").
}
\details{
\tabular{ll}{
Package: \tab MCAvariants\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-06-27\cr
%License: \tab What license is it under?\cr
}
}
\author{
Rosaria Lombardo
Maintainer: Rosaria Lombardo <rosaria.lombardo@unina2.it>
}
%\details{
%  ~~ This function belongs to the object class called "cacorporateplus"  ~~
%}
\value{
Description of the output returned
\item{DataMatrix }{The starting table of variables in reduced (primitive) code.}
\item{rows }{The row number of the starting table.}
\item{cols }{The column number of the starting table (coincident with the variable number).}
\item{rowlabels }{The label of the row individuals.}
\item{columnlabels }{The label of the column variable categories.}
\item{Rprinccoord }{The  coordinates of individuals.}
\item{Cprinccoord }{The category variable coordinates.}
\item{inertiaXsum}{The total inertia when multiple  correspondence analysis is performed on the indicator table.}
\item{inertiaBurtsum}{The total inertia when multiple  correspondence analysis is performed on the Burt table.}
\item{inertias}{ Benzecri's Adjusted Inertia values, percentages and cumulative values.}
\item{inertiasAdjusted}{The adjusted inertia values.}
\item{catype}{The kind of multiple correspondence analysis chosen, classical or ordered, that is \code{catype} is \dQuote{mca} or \dQuote{omca}.}
\item{printdims}{The  dimension of a matrix in print. By default it is equal to 3.}
\item{comp}{The polynomial components of inertia when \code{catype} is \dQuote{omca}. The total inertia is partitioned in terms
of polynomial components.}
\item{componentpvalue1}{The p-value of the polynomial components of total inertia, when \code{catype} is \dQuote{omca}.}
\item{degreef}{The degree of freedom of polynomial components of total inertia when, \code{catype} is \dQuote{omca}.}
}

\references{
Lombardo R and Meulman JJ (2010) Journal of Classification, 27, 191-210. \cr
Beh EJ Lombardo R (2014) Correspondence Analysis, Theory, Practice and New Strategies. Wiley
}

%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
satisfaction<-structure(c(5 , 3 , 4 , 4 , 4 , 5 , 4 , 3 , 5 , 4 , 
3 , 3 , 4 , 4 , 4 , 3 , 3 , 1 , 4 , 3 , 3 , 2 , 3 , 4 , 2 , 3 , 
3 , 4 , 4 , 2 , 4 , 3 , 3 , 3 , 4 , 3 , 2 , 3 , 3 , 5 , 5 , 5 , 
4 , 2 , 4 , 5 , 3 , 3 , 4 , 3 , 4 , 2 , 3 , 4 , 4 , 5 , 1 , 3 , 
4 , 4 , 5 , 4 , 5 , 3 , 4 , 3 , 5 , 4 , 3 , 2 , 4 , 3 , 3 , 3 , 
3 , 3 , 4 , 3 , 4 , 4 , 2 , 3 , 4 , 3 , 4 , 4 , 4 , 5 , 3 , 3 , 
2 , 5 , 3 , 3 , 4 , 5 , 2 , 3 , 3 , 1 , 4 , 4 , 1 , 4 , 5 , 3 , 
3 , 3 , 2 , 4 , 1 , 1 , 1 , 4 , 4 , 4 , 4 , 1 , 3 , 5 , 1 , 3 , 
3 , 3 , 5 , 2 , 2 , 5 , 3 , 4 , 5 , 3 , 4 , 3 , 4 , 4 , 4 , 5 , 
3 , 5 , 4 , 3 , 4 , 3 , 4 , 5 , 3 , 5 , 4 , 4 , 3 , 3 , 4 , 5 , 
3 , 3 , 2 , 3 , 3 , 3 , 4 , 4 , 3 , 2 , 5 , 3 , 3 , 3 , 3 , 3 , 
2 , 2 , 3 , 5 , 5 , 3 , 3 , 3 , 3 , 4 , 2 , 2 , 3 , 1 , 4 , 4 , 
3 , 5 , 3 , 3 , 4 , 4 , 3 , 3 , 3 , 4 , 5 , 3 , 3 , 3 , 3 , 4 , 
5 , 5 , 3 , 4 , 2 , 3 , 4 , 3 , 3 , 4 , 2 , 1 , 3 , 4 , 3 , 3 , 
3 , 2 , 3 , 3 , 4 , 4 , 3 , 3 , 3 , 3 , 3 , 4 , 4 , 3 , 4 , 5 , 
4 ,5 , 3 , 5 , 5 , 3 , 5 , 4 , 4 , 3 , 3 , 4 , 3 , 
4 , 3 , 4 , 3 , 4 , 2 , 3 , 5 , 2 , 2 , 2 , 3 , 4 , 2 , 5 , 4 , 
4 , 1 , 4 , 3 , 4 , 2 , 4 , 4 , 2 , 3 , 4 , 5 , 5 , 5 , 5 , 3 , 
3 , 4 , 3 , 4 , 5 , 3 , 5 , 3 , 4 , 4 , 3 , 4 , 3 , 3 , 4 , 4 , 
5 , 4 , 5 , 2 , 5 , 2 , 5 , 4 , 4 , 5 , 4 , 2 , 5 , 3 , 4 , 3 , 
5 , 4 , 3 , 3 , 2 , 3 , 5 , 4 , 4 , 5 , 3 , 5 , 4 , 3 , 4 , 5 , 
3 , 3 , 4 , 5 , 3 , 4 , 2 , 2 , 4 , 4 , 1 , 5 , 5 , 5 , 4 , 5 , 
2 , 4 , 3 , 2 , 2 , 4 , 4 , 5 , 4 , 4 , 4 , 5 , 1 , 4 , 3 , 3 , 
5 , 3 , 4 , 5 , 4 , 5 , 5 , 4 , 3 , 3 , 5 , 4 , 5 , 4 , 3 , 4 , 
5 , 3 , 4 , 5 , 5 , 5 , 4 , 4 , 4 , 4 , 5 , 2 , 5 , 5 , 3 , 3 , 
3 , 4 , 3 , 3 , 4 , 3 , 2 , 2 , 4 , 4 , 2 , 3 , 3 , 3 , 1 , 1 , 
4 , 5 , 4 , 4 , 3 , 4 , 4 , 3 , 1 , 2 , 2 , 4 , 5 , 4 , 5 , 5 , 
4 , 2 , 4 , 4 , 3 , 3 , 4 , 4 , 4 , 4 , 5 , 5 , 5 , 5 , 5 , 4 , 
2 , 4 , 4 , 3 , 4 , 2 , 4 , 4 , 3 , 1 , 1 , 2 , 5 , 1 , 5 , 2 , 
3 , 4 , 5 , 4 , 4 , 4 , 4 , 5 , 2 , 1 , 4 , 3 , 4 , 5 , 4 , 
5 , 5 , 5 , 4 , 2 , 5 , 4 , 3 , 3 , 4 , 3 , 3 , 
4 , 5 , 4 , 3 , 2 , 3 , 3 , 4 , 2 , 2 , 1 , 3 , 4 , 2 , 4 , 
4 , 4 , 1 , 4 , 3 , 3 , 4 , 4 , 4 , 2 , 2 , 4 , 5 , 5 , 5 , 
4 , 3 , 3 , 5 , 3 , 4 , 4 , 3 , 5 , 4 , 5 , 3 , 5 , 5 , 3 , 
3 , 4 , 5 , 5 , 3 , 5 , 3 , 4 , 3 , 5 , 4 , 3 , 4 , 4 , 3 , 
4 , 2 , 3 , 3 , 5 , 4 , 3 , 3 , 2 , 3 , 4 , 3 , 5 , 5 , 4 , 
5 , 3 , 4 , 4 , 5 , 3 , 2 , 4 , 5 , 3 , 4 , 3 , 1 , 4 , 4 , 
1 , 3 , 5 , 5 , 4 , 4 , 2 , 4 , 5 , 4 , 2 , 5 , 3 , 4 , 3 , 
2 , 3 , 5 , 2 , 4 , 3 , 2 , 5 , 4 , 4 , 4 , 3 , 3 , 5 , 3 , 
3 , 3 , 5 , 3 , 4 , 5 , 2 , 5 , 5 , 4 , 4 , 5 , 4 , 4 , 3 , 
5 , 4 , 3 , 5 , 3 , 4 , 4 , 2 , 3 , 3 , 5 , 3 , 3 , 3 , 3 , 
2 , 1 , 3 , 3 , 3 , 3 , 3 , 2 , 1 , 1 , 4 , 5 , 4 , 3 , 3 , 
3 , 3 , 2 , 1 , 2 , 2 , 4 , 5 , 5 , 4 , 5 , 4 , 3 , 4 , 4 , 
3 , 4 , 5 , 3 , 4 , 1 , 4 , 4 , 5 , 3 , 4 , 5 , 3 , 4 , 4 , 
3 , 4 , 3 , 3 , 1 , 3 , 1 , 4 , 4 , 5 , 2 , 5 , 2 , 2 , 3 , 
2 , 3 , 5 , 3 , 3 , 5 , 2 , 4 , 4 , 3 , 3 , 5 , 4 ,5 , 
4 , 5 , 4 , 2 , 4 , 4 , 4 , 2 , 3 , 3 , 2 , 4 , 5 , 4 , 3 , 
2 , 3 , 4 , 5 , 2 , 2 , 1 , 3 , 4 , 3 , 3 , 3 , 4 , 2 , 4 , 
3 , 3 , 2 , 4 , 4 , 2 , 2 , 4 , 4 , 5 , 5 , 5 , 3 , 4 , 5 , 
3 , 4 , 5 , 3 , 5 , 4 , 5 , 4 , 5 , 5 , 4 , 3 , 4 , 4 , 5 , 
4 , 5 , 3 , 5 , 3 , 5 , 4 , 3 , 5 , 4 , 3 , 5 , 3 , 3 , 3 , 
5 , 3 , 4 , 3 , 2 , 3 , 4 , 3 , 5 , 5 , 5 , 5 , 3 , 3 , 4 , 
5 , 4 , 3 , 4 , 4 , 3 , 4 , 3 , 1 , 3 , 4 , 1 , 3 , 5 , 5 , 
5 , 4 , 2 , 4 , 3 , 4 , 3 , 4 , 2 , 5 , 4 , 5 , 4 , 5 , 2 , 
4 , 2 , 3 , 5 , 4 , 3 , 4 , 4 , 4 , 4 , 4 , 4 , 4 , 3 , 3 , 
5 , 4 , 3 , 5 , 5 , 3 , 4 , 2 , 5 , 5 , 3 , 5 , 4 , 3 , 4 , 
3 , 4 , 4 , 2 , 3 , 3 , 4 , 3 , 4 , 4 , 3 , 3 , 2 , 4 , 3 , 
4 , 3 , 3 , 3 , 1 , 2 , 5 , 5 , 4 , 3 , 3 , 2 , 3 , 2 , 1 , 
2 , 3 , 3 , 5 , 5 , 4 , 5 , 3 , 5 , 5 , 4 , 4 , 4 , 5 , 2 , 
4 , 3 , 3 , 3 , 5 , 4 , 4 , 5 , 4 , 5 , 4 , 4 , 4 , 3 , 4 , 
2 , 4 , 2 , 3 , 3 , 5 , 2 , 3 , 2 , 3 , 4 , 4 , 4 , 3 , 4 , 
3 , 4 , 3 , 4 , 5 , 3 , 4 , 3 , 4 , 5 , 4 , 5 , 
4 , 3 , 4 , 4 , 4 , 3 , 4 , 3 , 2 , 3 , 5 , 3 , 4 , 1 , 2 , 
3 , 4 , 2 , 1 , 1 , 1 , 4 , 3 , 3 , 3 , 5 , 1 , 4 , 3 , 3 , 
3 , 2 , 3 , 2 , 2 , 4 , 5 , 5 , 5 , 4 , 2 , 5 , 5 , 5 , 4 , 
4 , 3 , 3 , 2 , 3 , 3 , 5 , 5 , 3 , 3 , 4 , 4 , 5 , 4 , 4 , 
5 , 4 , 3 , 5 , 4 , 3 , 3 , 4 , 4 , 4 , 2 , 3 , 3 , 5 , 2 , 
4 , 2 , 1 , 3 , 4 , 3 , 4 , 5 , 4 , 5 , 3 , 3 , 3 , 5 , 2 , 
3 , 3 , 5 , 4 , 3 , 3 , 1 , 4 , 5 , 2 , 3 , 5 , 4 , 5 , 3 , 
3 , 4 , 3 , 3 , 2 , 5 , 3 , 5 , 3 , 4 , 4 , 5 , 2 , 4 , 3 , 
3 , 5 , 2 , 3 , 4 , 4 , 3 , 3 , 4 , 4 , 3 , 2 , 3 , 5 , 5 , 
1 , 5 , 5 , 3 , 4 , 2 , 4 , 5 , 3 , 5 , 4 , 3 , 5 , 3 , 4 , 
4 , 2 , 4 , 3 , 4 , 3 , 3 , 4 , 2 , 2 , 2 , 4 , 3 , 3 , 3 , 
3 , 3 , 2 , 1 , 5 , 5 , 2 , 3 , 3 , 3 , 3 , 2 , 1 , 1 , 3 , 
3 , 4 , 5 , 3 , 4 , 4 , 2 , 5 , 3 , 4 , 4 , 3 , 2 , 4 , 1 , 
4 , 4 , 3 , 4 , 4 , 5 , 2 , 3 , 4 , 4 , 4 , 3 , 3 , 1 , 3 , 
1 , 3 , 4 , 5 , 2 , 4 , 2 , 3 , 3 , 5 , 4 , 4 , 4 , 2 , 4 , 
2 , 4 , 3 , 3 , 3 , 4 , 4 ),.Dim=c(235,5),
.Dimnames=list(paste("i",1:235, sep=""), c("TANG", "REL", 
"CRES", "CASS", "EMPAT")))
MCAvariants(satisfaction, catype = "omca")
MCAvariants(satisfaction, catype = "mca")
}
\keyword{multivariate}

