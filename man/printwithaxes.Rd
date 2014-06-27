\name{printwithaxes}
\alias{printwithaxes}
\title{
Secondary printing function 
}
\description{
The function is called from the main print function \code{printmcacorporateris}. 
It adds the names to objects. 
}
\usage{
printwithaxes(res, thenames)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
An \code{R} object.}
  \item{thenames}{
A character vector of up to the same length as \code{x}.}
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
It is called from \code{printmcacorporateris}.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.

%## The function is currently defined as
%function (res, thenames) 
%{
%    names(res) <- thenames
%    print(res, digits = 4)
%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
