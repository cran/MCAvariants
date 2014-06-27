\name{mcabasic}
\alias{mcabasic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Classical multiple correspondence analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is used in the main function \code{MCAvariants} when 
the input parameter is \code{catype="mca"}. 
%It performs the singular value decomposition of the weighted super-indicator matrix and
%compute principal axes, coordinates, weights of rows and columns and
%total inertia.
}
\usage{
mcabasic(xo, np = 2, nmod = c(5, 5), tmod = 10, rows = 511, idr, idc, idcv)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xo}{
The starting table of variables in reduced code.
}
  \item{np}{
The column number of the starting table (coincident with the variable number).
}
  \item{nmod}{
The number of variable catgories of each variable.}
  \item{tmod}{
The total number of variable catgories.}
  \item{rows}{
The row number of the starting table (coincident with the individual number).}
  \item{idr}{
The row labels of the data table.}
  \item{idc}{
The column labels of the data table.
}
  \item{idcv}{
The labels of the categories of each variable.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Lombardo R and Meulman JJ (2010) Journal of Classification, 27, 191-210. \cr
Beh EJ Lombardo R (2014) Correspondence Analysis, Theory, Practice and New Strategies. Wiley}
\author{
Rosaria Lombardo}
\note{
This function belongs to the \code{R} object class called \code{mcabasicresults}.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
\keyword{multivariate}% __ONLY ONE__ keyword per line
