\name{MCAvariants}
\alias{MCAvariants}
\docType{package}
\title{
Classic and Ordered Multiple Correspondence Analysis
}
\description{
It performs Classic Multiple Correspondence analysis for nominal variables
(setting \code{catype} = "mca") and Ordered Multiple Correspondence analysis via orthogonal polynomials (setting \code{catype}="omca").
When the categorical variables are nominal and ordinal, you can specify writing \code{FALSE} or \code{TRUE} in the input parameter \code{vordered}. 
}
\usage{
MCAvariants(Xtable, catype = "omca", np = 5, vordered=c(TRUE,TRUE,TRUE,TRUE,TRUE)) 
}
\arguments{
\item{Xtable}{
The two-way contingency table.
}
\item{catype}{
The input parameter for specifying what variant of multiple correspondence analysis is considered.  
By default,  \code{catype = "mca"}. 
}
\item{np}{The input parameter for specifying the number of categorical variables. By default, \code{np = 5}.
}
\item{vordered}{
The flag parameter for specifying what variable  is ordered, the ordered variables should be in column close each other.
By default,  all the five variables  are ordered: \code{vordered = c(TRUE,TRUE,TRUE,TRUE,TRUE)}. 
}
}
%\details{
%This function belongs to the object class called \code{mcacorporateplus}  
%}
\value{
Description of the output returned
\item{Xtable}{The starting table of variables in reduced (primitive) code.}
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
\author{
Rosaria Lombardo and Eric J Beh
}
\note{
This function recalls internally two other functions, depending on the setting of the input parameter \code{catype}, it recalls
multiple correspondence analysis or ordered multiple correspondence analysis.
It gives the output object necessary for printing and plotting the results. These two
important functions are \code{print.MCAvariants} and \code{plot.MCAvariants}. 
}
%\examples{
%data(satisfaction)
%MCAvariants(satisfaction, catype = "omca",np=5,vordered=c#%(TRUE,TRUE,TRUE,TRUE,TRUE))
%MCAvariants(satisfaction, catype = "mca",np=5)
%}
\keyword{multivariate}

