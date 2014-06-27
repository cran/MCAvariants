\name{omcabasic}
\alias{omcabasic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Ordered multiple correspondence analysis via orthogonal polynomials
}
\description{
This function is used in the main function \code{MCAvariants} when 
the input parameter is \code{catype="omca"}. It requires that all categorical variables are ordered variables.
It performs the hybrid decomposition of the weighted super-indicator matrix and
compute polynomial axes, coordinates, weights of rows and columns and
total inertia.
}
\usage{
omcabasic(xo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xo}{
The starting table of variables in reduced code.
}
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
